% forest_segmentation_pipeline.m
% Vstup: rastrový obraz topografické mapy (TM25)
% Výstup: binary mask (forest_mask), pixel coordinates (rows,cols) a polygony pro export

clc; clear; close all;

%% --- PARAMETERS (ladit podle mapy) ----------------
imgPath = 'data/TM25_sk3.jpg';
kClusters = 5;         % počet clusterů v k-means (experimentovat 4-6)
gaborMinWL = 4/sqrt(2);
gaussSigma = 2;
minHolePxToKeep = 50;  % díry větší než toto zachovat (nevyplňovat)
minObjectArea = 20;    % malé artefakty odstranit
prusek_max_width = 6;  % max šířka (px) pro průsek (odhadem) 
prusek_min_length = 30; % minimální délka (px) pro identifikaci průseku

%% --- READ & PREP -----------------------------------
I = imread(imgPath);
if size(I,3) == 3
    lab = rgb2lab(I);
else
    lab = cat(3, rgb2lab(repmat(I,1,1,3)));
end
L = lab(:,:,1);
a = lab(:,:,2);
b = lab(:,:,3);

% standardize / convert to single
a_s = im2single(a);
b_s = im2single(b);

%% --- SMOOTH & TEXTURE FEATURES ---------------------
% Gauss smoothing to remove tiny speckle
a_blur = imgaussfilt(a_s,gaussSigma);

% Gabor bank - create wavelengths & orientations (adapted z vašeho skriptu)
[numRows,numCols,~] = size(I);
wavelengthMin = gaborMinWL;
wavelengthMax = hypot(numRows,numCols);
n = floor(log2(wavelengthMax/wavelengthMin));
wavelength = 2.^(0:(max(1,n-2))) * wavelengthMin;
orientations = 0:30:150;  % více orientací než 45deg krok

g = gabor(wavelength, orientations);
gabormag = imgaborfilt(a_blur, g);   % gabormag is height x width x numFilters
% aggregate Gabor responses to a single texture map (max energy)
gaborEnergy = max(gabormag,[],3);
gaborEnergy = im2single(mat2gray(gaborEnergy)); % normalize

% local STD (texture cue)
localStd = stdfilt(a_blur, ones(5));

%% --- BUILD FEATURE VECTOR & PCA (optional) ----------
% Use channels: a_blur, b_s, gaborEnergy, localStd
F = cat(3, mat2gray(a_blur), mat2gray(b_s), gaborEnergy, mat2gray(localStd));
[nr,nc,nf] = size(F);
X = reshape(F, nr*nc, nf);

% optional: PCA to 3 dims for stability
[coeff,score,~] = pca(X);
Xred = score(:,1:3);

%% --- K-MEANS SEGMENTATION --------------------------
[Lseg,centers] = imsegkmeans(reshape(Xred, nr, nc, 3), kClusters, NumAttempts=10);
% Lseg is label image (1..kClusters)

figure; imshow(label2rgb(Lseg)); title('Raw k-means segmentation');

%% --- SELECT CLUSTER(S) CORRESPONDING TO FOREST -----
% Heuristika: cluster s tmavšími a/b hodnotami + vyšší texturou
% Compute mean a,b,gaborEnergy per cluster:
maskForestGuess = false(nr,nc);
stats = zeros(kClusters,4);
for k=1:kClusters
    idx = (Lseg==k);
    stats(k,1) = mean(a(idx),'all','omitnan');
    stats(k,2) = mean(b(idx),'all','omitnan');
    stats(k,3) = mean(gaborEnergy(idx),'all','omitnan');
    stats(k,4) = mean(localStd(idx),'all','omitnan');
end
% Výběr clusteru: nízký L (tmavé) => v lab to může být nízké L,
% ale jsme použili spíše a/b a texturu. Níže jednoduchá heuristika:
[~, forestCluster] = min(stats(:,1)); % nejtmavší v L canal approximation
% pokud to selže, můžeš ručně vybrat:
% forestCluster = 3;

maskForestGuess = (Lseg == forestCluster);

figure; imshow(maskForestGuess); title('Initial forest guess');

%% --- REMOVE LINEAR FEATURES (vrstevnice) -----------
% Extract likely line-map by thresholding strong Gabor responses
lineMap = gaborEnergy > graythresh(gaborEnergy) * 0.9; % práh (ladit)
% Morphologically thin; then dilate slightly to cover line width
lineMapThin = bwmorph(lineMap,'thin',Inf);
lineMapDil = imdilate(lineMapThin, strel('disk',1));

% Remove lines from forest mask (maskForestGuess)
mask_no_lines = maskForestGuess & ~lineMapDil;

figure; imshowpair(maskForestGuess, mask_no_lines,'montage');
title('Before / After removing linear features');

%% --- FILL SMALL HOLES (<= 50 px) & preserve bigger holes ----
% Find holes in current mask (holes = background components inside object)
% Approach: find connected components of inverted mask that's inside bounding box
holesMask = ~mask_no_lines;
CC = bwconncomp(holesMask,8);
S = regionprops(CC, 'Area', 'PixelIdxList', 'BoundingBox');

% We want holes that are fully enclosed by the forest object. A reliable way:
% For each connected component in holesMask check if its pixels touch image border.
holeIsBoundary = false(CC.NumObjects,1);
for i=1:CC.NumObjects
    [r,c] = ind2sub([nr,nc], CC.PixelIdxList{i});
    if any(r==1 | r==nr | c==1 | c==nc)
        holeIsBoundary(i) = true;
    end
end

% Consider only true holes (not background touching border) and with size <= threshold -> fill them
maskFilled = mask_no_lines;
for i=1:CC.NumObjects
    if ~holeIsBoundary(i)
        area_i = S(i).Area;
        if area_i <= minHolePxToKeep
            maskFilled(S(i).PixelIdxList) = true; % fill small hole
        end
    end
end

figure; imshowpair(mask_no_lines, maskFilled,'montage');
title('Before / After filling small holes <= 50 px');

%% --- REMOVE SMALL OBJECTS (artefakty) ----------------
maskClean = bwareaopen(maskFilled, minObjectArea);

%% --- IDENTIFY & REMOVE "LESNI PRUSEKY" (úzké průseky) -------
% Goal: v masce najít objekty uvnitř lesa, které jsou dlouhé a úzké
% We'll invert maskClean to find background components within the forest and test geometry.
forestCC = bwconncomp(maskClean,8);
forestProps = regionprops(forestCC, 'Area', 'PixelIdxList', 'BoundingBox', 'Solidity', 'Eccentricity');

% Create a mask that removes elongated clearings (pruseky) from forest:
maskFinal = maskClean;

% Strategy: find background components inside forest (holes) with small width and large length:
invMask = ~maskClean;
CCholes = bwconncomp(invMask,8);
holeProps = regionprops(CCholes, 'Area','BoundingBox','MajorAxisLength','MinorAxisLength','Eccentricity');

for i=1:CCholes.NumObjects
    a_major = holeProps(i).MajorAxisLength;
    a_minor = holeProps(i).MinorAxisLength;
    ecc = holeProps(i).Eccentricity;
    area = holeProps(i).Area;
    % condition: elongated (minor small), length long, area not too big
    if (a_minor <= prusek_max_width) && (a_major >= prusek_min_length) && (ecc > 0.9)
        % treat as prusek -> remove from forest: ensure we subtract these holes from forest mask
        maskFinal(CCholes.PixelIdxList{i}) = false;  % ensure hole remains (i.e., not forest)
        % Also optionally dilate hole slightly to remove fringes
    end
end

figure; imshowpair(maskClean, maskFinal,'montage');
title('Before / After removing pruseky');

%% --- OPTIONAL: Clean small islands & smooth boundaries ------------
maskFinal = imopen(maskFinal, strel('disk',2));
maskFinal = bwareaopen(maskFinal, minObjectArea);

figure; imshow(maskFinal); title('Final forest mask');

%% --- OUTPUT: pixel coordinates of forest mask --------------------
[rows, cols] = find(maskFinal);  % row = Y, col = X in image coords
pixel_coords = [cols, rows];     % [x,y] pairs - common order (col, row)

fprintf('Detected forest pixels: %d\n', size(pixel_coords,1));

% Save pixel coordinates
save('forest_pixel_coords.mat', 'pixel_coords', 'maskFinal');

%% --- OPTIONAL: Polygonize & Export (outline) --------------------
% Convert mask to polygons (one polygon per connected component)
BW = maskFinal;
[B,L] = bwboundaries(BW, 'noholes');

% Build polyshapes (pixel coords); convert to map coords later using affine transform
polys = cellfun(@(b) [b(:,2), b(:,1)], B, 'UniformOutput', false); % col,row -> x,y

% Save polygons as MAT
save('forest_polygons_pixels.mat', 'polys');

disp('Pipeline finished. Results saved to forest_pixel_coords.mat and forest_polygons_pixels.mat');
