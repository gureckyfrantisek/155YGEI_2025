%% --------------------------------------------------------------------
% Rychlá segmentace lesa z topografické mapy pomocí k-means
% Doba výpočtu: 5–20 sekund
%% --------------------------------------------------------------------

I = imread('data/TM25_sk3.jpg');
I = im2double(I);
[H, W, ~] = size(I);

%% 1) Převod na L*a*b* (lepší pro k-means)
lab = rgb2lab(I);
L = lab(:,:,1);
A = lab(:,:,2);
B = lab(:,:,3);

%% 2) Příprava dat pro kmeans – pouze subsampling
N = H * W;

sample_size = 150000;        % 150k pixelů – rychlé a stabilní
idx = randperm(N, sample_size);

features = [L(:), A(:), B(:)];
features_sample = features(idx, :);

%% 3) K-means clustering (3 clustery)
k = 3;
opts = statset('MaxIter', 150, 'Display', 'final');

[idx_sample, C] = kmeans(features_sample, k, ...
    'Distance', 'sqeuclidean', ...
    'Replicates', 3, ...
    'Options', opts);

%% 4) Přiřazení clusterů pro celý obrázek
dist = pdist2(features, C);
[~, idx_full] = min(dist, [], 2);
clusters = reshape(idx_full, H, W);

%% 5) Identifikace clusteru lesa
% Les = nejzelenější cluster → nejnižší A a nejvyšší B
cluster_stats = zeros(k,2);
for i = 1:k
    mask = (clusters == i);
    cluster_stats(i,1) = mean(A(mask));  % průměr A
    cluster_stats(i,2) = mean(B(mask));  % průměr B
end

% les = nejmenší A (více do zelena)
[~, forest_cluster] = min(cluster_stats(:,1));

%% 6) Získání lesní masky
mask_forest = (clusters == forest_cluster);

%% 7) Odstranění vrstevnic a černých cest
HSV = rgb2hsv(I);
V = HSV(:,:,3);
mask_forest(V < 0.25) = 0;  % černé linie ven

%% 8) Morfologické dočištění
mask_forest = imopen(mask_forest, strel('disk', 3));
mask_forest = imclose(mask_forest, strel('disk', 5));

%% 9) Otvory větší než 50 px zachovat
mask_forest = imfill(mask_forest, 'holes');
holes = ~mask_forest;
Lholes = bwlabel(holes);
props = regionprops(Lholes, 'Area');

for i = 1:length(props)
    if props(i).Area < 50
        mask_forest(Lholes == i) = 1;
    end
end

%% 10) Výsledek
figure; imshow(mask_forest); title('K-means maska lesa');

%% 11) Export pixelů
[y, x] = find(mask_forest==1);
forest_pixels = [x,y];
save('forest_pixels.mat','forest_pixels');
