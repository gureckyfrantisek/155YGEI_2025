%% === Načtení ===
I = imread('data/TM25_sk3.jpg');
Igray = im2double(rgb2gray(I));

%% === RYCHLÝ STD FILTR ===
window = 9;
h = ones(window) / window^2;

mu = imfilter(Igray, h, 'replicate');
mu2 = imfilter(Igray.^2, h, 'replicate');
T = sqrt(mu2 - mu.^2);      % místní směrodatná odchylka
T = mat2gray(T);

%% === Zesílení textury (nutné pro váš typ mapy) ===
T = T .^ 2.5;
T = imgaussfilt(T, 1) - imgaussfilt(T, 6);
T = mat2gray(T);

%% === Segmentace K-means ===
X = T(:);
idx = kmeans(X,2,'Replicates',3);

mask = reshape(idx,size(T));

m1 = mean(T(mask==1));
m2 = mean(T(mask==2));

if m1 < m2
    forestMask = mask==1;
else
    forestMask = mask==2;
end

%% Čištění
forestMask = imopen(forestMask, strel('disk', 2));
forestMask = imclose(forestMask, strel('disk', 4));
forestMask = bwareaopen(forestMask, 50);

%% Výstup souřadnic
[y,x] = find(forestMask);
forestPixels = [x,y];

%% Zobrazení
figure;
subplot(131); imshow(I); title('Mapa');
subplot(132); imshow(T,[]); title('Rychlá STD textura');
subplot(133); imshow(forestMask); title('Segmentovaný les');
