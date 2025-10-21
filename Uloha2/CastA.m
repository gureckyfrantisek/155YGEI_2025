% Implementujte algoritmus pro vyhledávání vzorů v Müllerově mapě Čech na základě obrazové korelace. 
% Vyberte jeden z objektů typu obec s kostelem jako okno pro vyhledávání a na základě vhodné hodnoty 
% korelace vyhledejte všechny pozice obcí s kostelem na mapovém listu. 
% Výsledek odevzdejte jako pixelové souřadnice kostelů.

% zredukovat pocet najitych koralaci na jendom misto pouze na jedno


clc; clear;

%% copilot !!!!!!!!!NEVERIT!!!!!!!! 
% Načtěte Müllerovu mapu a vyberte okno pro vyhledávání
mapImage = imread('muller_map.png');
searchWindow = imcrop(mapImage); % Uživatel vybere okno
correlationThreshold = 0.8; % Nastavte prahovou hodnotu korelace
% Proveďte obrazovou korelaci a najděte pozice obcí s kostelem
correlationMap = normxcorr2(searchWindow(:,:,1), mapImage(:,:,1)); % Pouze pro jeden kanál
[maxCorr, maxIdx] = max(correlationMap(:)); % Najděte maximální korelaci
if maxCorr > correlationThreshold
    [yPeak, xPeak] = ind2sub(size(correlationMap), maxIdx); % Získejte souřadnice
    churchCoordinates = [xPeak - size(searchWindow, 2) + 1, yPeak - size(searchWindow, 1) + 1]; % Uložení souřadnic
end
% Uložení souřadnic kostelů do pole, pokud je nalezeno více než jedno místo
if maxCorr > correlationThreshold
    churchCoordinatesList = [churchCoordinates]; % Inicializace seznamu souřadnic
    % Procházejte zbytek correlationMap a hledejte další shody
    for i = 1:numel(correlationMap)
        if correlationMap(i) > correlationThreshold
            [y, x] = ind2sub(size(correlationMap), i);
            newCoordinates = [x - size(searchWindow, 2) + 1, y - size(searchWindow, 1) + 1];
            churchCoordinatesList = [churchCoordinatesList; newCoordinates]; % Přidání nových souřadnic
        end
    end
end
% Odstranit duplicity ze seznamu souřadnic kostelů
churchCoordinatesList = unique(churchCoordinatesList, 'rows');
% Uložení souřadnic kostelů do souboru, pokud byly nalezeny
if ~isempty(churchCoordinatesList)
    writematrix(churchCoordinatesList, 'church_coordinates.csv');
end

% Uložení souřadnic kostelů do proměnné pro další zpracování
numChurches = size(churchCoordinatesList, 1);
disp(['Nalezeno kostelů: ', num2str(numChurches)]);