%% segmentace_lesa_TM25.m
% Robustní extrakce plochy lesa z naskenované mapy TM25_sk3.jpg
% Výstup: forest_pixels - Nx2 pole [row, col]

clear; close all; clc;

% ------- 1) Načtení obrazu -------
I = imread('data/TM25_sk3.jpg');
I = im2double(I);
[H, W, ~] = size(I);

% Zálohovaný zobrazení
figure('Name','Originál'); imshow(I); title('Originál TM25');

% ------- 2) Převod do HSV a základní prahy pro "zelené" -------
hsvI = rgb2hsv(I);
HUE = hsvI(:,:,1);   % 0..1
SAT = hsvI(:,:,2);
VAL = hsvI(:,:,3);

% Prahové hodnoty naladěné na tvůj sken (tuning z histogramu)
hue_min = 0.16;  hue_max = 0.50;   % zelené odstíny
sat_min = 0.06;                  % zelené v mapě mohou mít nižší sytost
val_min = 0.35;                  % příliš tmavé pixely vyloučit

mask_green = (HUE >= hue_min) & (HUE <= hue_max) & (SAT >= sat_min) & (VAL >= val_min);

% malé úpravy - odrušení šumu
mask_green = medfilt2(mask_green, [3 3]);

figure('Name','Hrubá zelená maska'); imshow(mask_green); title('Hrubá zelená maska');

% ------- 3) Detekce tenkých lineárních struktur (vrstevnice, kresba) -------
% Problém: vrstevnice/čáry jsou tenké a bělavě/béžově/brown; detekujeme tenké linky z grayscale
gray = rgb2gray(I);

% Zvýraznění lineárních světlých/kontrastních struktur pomocí top-hat s lineárním SE
line_len = 15;    % délka strukturujícího elementu (naladitelné)
angles = 0:15:165;
tophat_sum = zeros(H,W);

for ang = angles
    se = strel('line', line_len, ang);
    th = imtophat(gray, se);    % top-hat zvýrazní světlé lineární prvky dané orientace
    tophat_sum = tophat_sum + th;
end
% Normalize a prahování top-hat mapy
tophat_sum = mat2gray(tophat_sum);
line_mask = tophat_sum > 0.12;   % práh naladit (0.08..0.18)
line_mask = bwareaopen(line_mask, 20); % odstranit velmi malé artefakty

figure('Name','Detekované linky (top-hat)'); imshow(line_mask); title('Detekované linky (vrstevnice / kresba)');

% ------- 4) Odečtení detekovaných čar z hrubé zelené masky -------
mask_no_lines = mask_green & ~line_mask;

% ------- 5) Odstranění popisků / sítí / malých objektů, zachování větších ploch -------
% Nejprve odstraníme drobné ostrůvky, které nejsou částí lesa
min_area = 200;    % menší objekty odstraníme (naladitelné)
mask_clean = bwareaopen(mask_no_lines, min_area);

% Doplníme díry uvnitř větších oblastí
mask_filled = imfill(mask_clean, 'holes');

% Znovu obnovíme jen díry které chceme zachovat (>50 px)
holes = mask_filled & ~mask_clean;
holes_big = bwareaopen(holes, 50);
forest_mask = (mask_filled & ~holes) | holes_big;

% ------- 6) Doplňkové čištění: odstraň tenké mosty/přeseky způsobené mapovou kresbou -------
% Odebereme objekty s velkým poměrem perimeter/area (velmi "proužkovité")
cc = bwconncomp(forest_mask);
props = regionprops(cc,'Area','Perimeter');
keepIdx = false(size(props));
for k = 1:length(props)
    A = props(k).Area;
    P = props(k).Perimeter;
    if A == 0
        keepIdx(k) = false;
    else
        compactness = (P^2)/(4*pi*A); % =1 pro kruh; vyšší = více pruhů/tenké
        % Podmínka: udržovat objekty se slušnou plochou nebo ne příliš nepravidelné
        keepIdx(k) = (A >= 250) | (compactness < 15);
    end
end
% rekonstrukce masky s vybranými komponentami
labels = labelmatrix(cc);
forest_mask2 = ismember(labels, find(keepIdx));

% ------- 7) Finální morfologické doladění -------
% lehké zaoblení hranic a odstranění malých „ostříže"
forest_mask2 = imopen(forest_mask2, strel('disk',2));
forest_mask2 = imclose(forest_mask2, strel('disk',3));
forest_mask2 = bwareaopen(forest_mask2, 150);

% ------- 8) Výsledek a souřadnice pixelů -------
[row, col] = find(forest_mask2);
forest_pixels = [row, col];    % Nx2

fprintf('Počet pixelů lesa: %d\n', size(forest_pixels,1));
% uložit do souboru .mat
save('forest_pixels_TM25.mat','forest_pixels');

% ------- 9) Vizualizace překrytu na originálu -------
figure('Name','Výsledek překrytu');
imshow(I); hold on;
h = imshow(cat(3, ones(size(forest_mask2)), zeros(size(forest_mask2)), zeros(size(forest_mask2))));
set(h,'AlphaData',0.35*forest_mask2); % červený průhledný překryv
title('Extrahovaný les (červeně) překrytý na originálu');

figure('Name','Maska lesa'); imshow(forest_mask2); title('Finální maska lesa');

% KONEC
