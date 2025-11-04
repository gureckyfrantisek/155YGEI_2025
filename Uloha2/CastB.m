% Plocha lesa očištěná o zbytkovou kresbu

% Segmentovaný obraz mapy při použití Gaborova nebo STD filtru 

% Plocha lesa uložená jako vrstva v GeoPackage v systému UTM *
% * Uvažujte při georeferencování mapy souřadnicový systém S1942 (EPSG:28403)

clc; clear;


I=imread('data\TM25_sk3.jpg');

lab_I=rgb2lab(I);
ab=lab_I(:,:,2);
ab=im2single(ab);
% imshow(ab,[])
Iblur=imgaussfilt(ab,2);
% [L,C]=imsegkmeans(Iblur,3,NumAttempts=10);
% imshow(L,[])

stdflt=stdfilt(ab,ones(3));
ab=fftshift(log(1+abs(fft(ab))));
% [L,C]=imsegkmeans([ab,stdflt],3,NumAttempts=10);
% imshow(L,[])


[numRows,numCols,~] = size(I);

wavelengthMin = 4/sqrt(2);
wavelengthMax = hypot(numRows,numCols);
n = floor(log2(wavelengthMax/wavelengthMin));
wavelength = 2.^(0:(n-2)) * wavelengthMin;

deltaTheta = 45;
orientation = 0:deltaTheta:(180-deltaTheta);

g=gabor(wavelength,orientation);
[g,~]=imgaborfilt(ab,g)
g1=max(g,[],3);

[L,C]=imsegkmeans([ab,g1],5,NumAttempts=10);
imshow(L,[])

