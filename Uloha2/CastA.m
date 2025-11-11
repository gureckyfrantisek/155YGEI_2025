% Implementujte algoritmus pro vyhledávání vzorů v Müllerově mapě Čech na základě obrazové korelace. 
% Vyberte jeden z objektů typu obec s kostelem jako okno pro vyhledávání a na základě vhodné hodnoty 
% korelace vyhledejte všechny pozice obcí s kostelem na mapovém listu. 
% Výsledek odevzdejte jako pixelové souřadnice kostelů.

% zredukovat pocet najitych koralaci na jendom misto pouze na jedno

clc; clear;

mapImage = imread('data\MMC_sk3.jpg');
mapImage=rgb2gray(mapImage);
vzor_x=3125;
vzor_y=3365;
searchWindow = imcrop(mapImage,[vzor_x,vzor_y,35,92]);

correlationThreshold = 0.50; 

% patterns = zeros(35, 92, 5); 

coordinates = [
    3125, 3365;
    795, 1623;
    3620,3220;
    6418,3224;
    2773,5030
];

for i = 1:size(coordinates,1)
    vzor_x = coordinates(i, 1);
    vzor_y = coordinates(i, 2);
    patterns(:, :,i) = imcrop(mapImage, [vzor_x, vzor_y, 35, 92]);
end

averagePattern = mean(patterns, 3);
searchWindow = uint8(averagePattern);
imshow(searchWindow)

correlationMap = normxcorr2(searchWindow(:,:,1), mapImage(:,:,1));

binMap = correlationMap >= correlationThreshold;
localMax = imregionalmax(correlationMap) & binMap;

[yPeak, xPeak] = find(localMax);

sh = size(searchWindow, 1);
sw = size(searchWindow, 2);
churchCoordinatesList = [xPeak - sw + 1, yPeak - sh + 1,];

churchCoordinatesList = unique(churchCoordinatesList, 'rows');
if ~isempty(churchCoordinatesList)
    writematrix(churchCoordinatesList, 'church_coordinates.csv');
end

numChurches = size(churchCoordinatesList, 1);
disp(['Nalezeno kostelů: ', num2str(numChurches)]);

figure;
imshow(mapImage);
hold on;
for i = 1:size(churchCoordinatesList, 1)
    x = churchCoordinatesList(i, 1);
    y = churchCoordinatesList(i, 2);
    rectangle('Position', [x, y, size(searchWindow, 2), size(searchWindow, 1)], ...
        'EdgeColor', 'r', 'LineWidth', 1.5);
end
for j=1:size(coordinates,1)
    rectangle('Position', [coordinates(j,1), coordinates(j,2), size(searchWindow, 2), size(searchWindow, 1)], ...
            'EdgeColor', 'b', 'LineWidth', 1.5);
end

hold off;
title(['Nalezeno kostelů: ', num2str(numChurches)]);
