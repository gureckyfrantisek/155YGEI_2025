clc; clear; format short g; close;

%% Input data
image = imread("data/Image1.bmp");
imshow(image)

%% Split the image into separate RGB layers
R = double(image(:, :, 1));
G = double(image(:, :, 2));
B = double(image(:, :, 3));

%%%% JPEG Compression

%% Transform to YCbCr 
Y = 0.2990 * R + 0.5870 * G + 0.1140 * B;
Cb = -0.1687 * R -0.3313 * G + 0.5000 * B + 128;
Cr = 0.5 * R - 0.4187 * G - 0.0813 * B + 128;

% Quantisation matrix: Y
Qy = [  16 11 10 16 24 40 51 61;
        12 12 14 19 26 58 60 55;
        14 13 16 24 40 87 69 56;
        14 17 22 29 51 87 80 62;
        18 22 37 26 68 109 103 77;
        24 35 55 64 81 104 113 92;
        49 64 78 87 103 121 120 101;
        72 92 95 98 112 100 103 99];

% Chrominance matrix: Cb, Cr
Qc = [  17 18 24 47 66 99 99 99
        18 21 26 66 99 99 99 99
        24 26 56 99 99 99 99 99
        47 69 99 99 99 99 99 99
        99 99 99 99 99 99 99 99
        99 99 99 99 99 99 99 99
        99 99 99 99 99 99 99 99
        99 99 99 99 99 99 99 99];

%% Compression factor q
q = 50;

Qyf = 50 * Qy / q;
Qcf = 50 * Qc / q;

%% Store components
Y_old = Y;
Cb_old = Cb;
Cr_old = Cr;

%% Split into 8x8 tiles
% Input size
[m, n] = size(Y);

% Resampling
resSize = [2, 2];       % Try for [2, 1], [3, 3], [4, 4]
for i = 1:resSize(1):m-1
    for j = 1:resSize(2):n-1
        Ysub = Y(i:i+resSize(1)-1, j:j+resSize(2)-1);
        Cbsub = Cb(i:i+resSize(1)-1, j:j+resSize(2)-1);
        Crsub = Cr(i:i+resSize(1)-1, j:j+resSize(2)-1);
        
        Y(i:i+resSize(1)-1, j:j+resSize(2)-1) = mean(Ysub(:));
        Cb(i:i+resSize(1)-1, j:j+resSize(2)-1) = mean(Cbsub(:));
        Cr(i:i+resSize(1)-1, j:j+resSize(2)-1) = mean(Crsub(:));
    end
end

for row = 1:8:m-7
    for col = 1:8:n-7
        Ytile = Y(row : row+7, col : col+7);
        Cbtile = Cb(row : row+7, col : col+7);
        Crtile = Cr(row : row+7, col : col+7);

        % Apply the DCT
        Yt_dct = dct(Ytile);
        Cbt_dct = dct(Cbtile);
        Crt_dct = dct(Crtile);

        % Quantize the transformed rasters
        Yt_q = Yt_dct ./ Qyf;
        Cbt_q = Cbt_dct ./ Qcf;
        Crt_q = Crt_dct ./ Qcf;
        
        % Round the values
        Yt_qr = round(Yt_q);
        Cbt_qr = round(Cbt_q);
        Crt_qr = round(Crt_q);

        % Update the matrix
        Y(row : row+7, col : col+7) = Yt_qr;
        Cb(row : row+7, col : col+7) = Cbt_qr;
        Cr(row : row+7, col : col+7) = Crt_qr;
    end
end

%% Zig zag ordering of pixels
% Works only for square images
if (m == n)
    % Prepare the space
    Yzz = zigzag(Y, m, n);
    Cbzz = zigzag(Cb, m, n);
    Crzz = zigzag(Cr, m, n);
    
    %% Huffman encoding
    [Ycomp, Ydict] = huffman(Yzz);
    [Cbcomp, Cbdict] = huffman(Cbzz);
    [Crcomp, Crdict] = huffman(Crzz);

end

%%%% JPEG Decompression
%% Reverse Zig Zag and Huffman
if (m == n)
    %% Huffman decoding
    Ydecomp = invhuffman(Ycomp, Ydict);
    Cbdecomp = invhuffman(Cbcomp, Cbdict);
    Crdecomp = invhuffman(Crcomp, Crdict);

    %% Inverse zig zag
    % Prepare the space
    Ydezz = invzigzag(Ydecomp, m, n);
    Cbdezz = invzigzag(Cbdecomp, m, n);
    Crdezz = invzigzag(Crdecomp, m, n);
end

for row = 1:8:m-7
    for col = 1:8:n-7
        Ytile = Y(row : row+7, col : col+7);
        Cbtile = Cb(row : row+7, col : col+7);
        Crtile = Cr(row : row+7, col : col+7);

        % De - Quantize the transformed rasters
        Yt_q = Ytile .* Qyf;
        Cbt_q = Cbtile .* Qcf;
        Crt_q = Crtile .* Qcf;

        % Apply the IDCT
        Yt_idct = invdct(Yt_q);
        Cbt_idct = invdct(Cbt_q);
        Crt_idct = invdct(Crt_q);

        % Update the matrix
        Y(row : row+7, col : col+7) = Yt_idct;
        Cb(row : row+7, col : col+7) = Cbt_idct;
        Cr(row : row+7, col : col+7) = Crt_idct;
    end
end

%% Transform to RGB
Rdec = Y + 1.4020 * (Cr-128);
Gdec = Y - 0.3441 * (Cb-128) - 0.7141 * (Cr-128);
Bdec = Y + 1.7720 * (Cb-128) - 0.0001 * (Cr-128);

% Convert back to uint8
Ri = uint8(Rdec);
Gi = uint8(Gdec);
Bi = uint8(Bdec);

% Assemble the image
imageOut(:, :, 1) = Ri;
imageOut(:, :, 2) = Gi;
imageOut(:, :, 3) = Bi;

% Show the image compressed
imshow(imageOut)

%% Standart deviations for RGB
% Calculate errors v
vR = R - Rdec;
vG = G - Gdec;
vB = B - Bdec;

% Sums of squared erros
sumvR = sum(vR.^2);
sumvG = sum(vG.^2);
sumvB = sum(vB.^2);

% Standart deviations
sigR = sqrt(sum(sumvR) / (m*n));
sigG = sqrt(sum(sumvG) / (m*n));
sigB = sqrt(sum(sumvB) / (m*n));

%% Custom functions
function [imgt] = dct(img)
    % Result will have the same size
    imgt = img;

    % Lines
    for u = 0:7
        % Columns
        for v = 0:7
            if u == 0
                Cu = sqrt(2)/2;
            else
                Cu = 1;
            end

            if v == 0
                Cv = sqrt(2)/2;
            else
                Cv = 1;
            end
            
            % Init the sum
            fuv = 0;
            
            for x = 0:7
                for y = 0:7
                    fuv = fuv + 0.25 * Cu * Cv * img(x+1, y+1) * cos((2*x + 1) * u * pi / 16) * cos((2*y + 1) * v * pi / 16);
                end
            end

            % Save the value
            imgt(u+1, v+1) = fuv;
        end
    end

end

function [imgt] = invdct(img)
    % Result will have the same size
    imgt = img;

    % Lines
    for x = 0:7
        % Columns
        for y = 0:7

            % Init the sum
            fxy = 0;
            for u = 0:7
                for v = 0:7
                    % Find Cu and Cv
                    if u == 0
                        Cu = sqrt(2)/2;
                    else
                        Cu = 1;
                    end
        
                    if v == 0
                        Cv = sqrt(2)/2;
                    else
                        Cv = 1;
                    end

                    fxy = fxy + 0.25 * Cu * Cv * img(u+1, v+1) * cos((2*x + 1) * u * pi / 16) * cos((2*y + 1) * v * pi / 16);
                end
            end

            % Save the value
            imgt(x+1, y+1) = fxy;
        end
    end

end

function [row] = zigzag(matrix, m, n)
    row = zeros(m*n, 1);
    added = 0;

    for i = 1:m
        if (mod(i, 2) == 0)
            % Zig
            for j = 0:i-1
                % Add to the front and back at once, 
                % if we're on the last i, add just from one side
                row(added+1) = matrix(1+j, i-j);
                if (i ~= m)
                    row(length(row) - added) = matrix(end-j, (m-i+1)+j);
                end

                added = added + 1;
            end
        else
            % Zag
            for j = 0:i-1
                row(added+1) = matrix(i-j, 1+j);
                if (i ~= m)
                    row(length(row) - added) = matrix((m-i+1)+j, end-j);
                end

                added = added + 1;
            end
        end
    end
end

function [matrix] = invzigzag(row, m, n)
    matrix = zeros(m, n);
    added = 0;

    for i = 1:m
        if (mod(i, 2) == 0)
            % Zig
            for j = 0:i-1
                matrix(1+j, i-j) = row(added+1);
                if (i ~= m)
                    matrix(end-j, (m-i+1)+j) = row(m*n - added);
                end

                added = added + 1;
            end
        else
            % Zag
            for j = 0:i-1
                matrix(i-j, 1+j) = row(added+1);
                if (i ~= m)
                    matrix((m-i+1)+j, end-j) = row(m*n - added);
                end

                added = added + 1;
            end
        end
    end
end

function generateCodes(node, prefix, dict)
    % node = {prob, symbol}  OR  {prob, {left, right}}
    value = node{2};

    % Case 1: leaf node (numeric symbol)
    if isnumeric(value)
        dict(prefix) = value;
        return;
    end

    % Case 2: internal node, recurse
    left = value{1};
    right = value{2};

    % Left gets 1, right gets 0
    generateCodes(left, [prefix '1'], dict);
    generateCodes(right, [prefix '0'], dict);
end

function [compressed, dict] = huffman(data)
    symbols = unique(data);
    counts = zeros(length(symbols), 1);
    for i = 1:length(data)
        index = find(data(i) == symbols);
        counts(index) = counts(index) + 1;
    end

    probs = counts / sum(counts);

    nodes = cell(length(symbols), 1);
    for i = 1:length(symbols)
        nodes{i} = {probs(i), symbols(i)};
    end
    
     while length(nodes) > 1
        % Sort both together
        [probs, order] = sort(probs);
        nodes = nodes(order);
    
        % Take two smallest nodes
        left = nodes{1};
        right = nodes{2};
    
        % Create new merged node
        newNode = {left{1} + right{1}, {left, right}};
    
        % Replace first two nodes with the merged one
        nodes = [nodes(3:end); {newNode}];
    
        % Update probability list
        probs = [probs(3:end); newNode{1}];
    end

    root = nodes{1};
    dict = containers.Map('KeyType', 'char', 'ValueType', 'double');
    generateCodes(root, '', dict);

    encodeDict = containers.Map('KeyType', 'double', 'ValueType', 'char');
    codes = keys(dict);
    for i = 1:length(codes)
        code = codes{i};
        symbol = dict(code);
        encodeDict(symbol) = code;
    end

    % Encode to a string
    compressed = '';
    for i = 1:length(data)
        compressed = [compressed encodeDict(data(i))];
    end
end

function [data] = invhuffman(compressed, dict)
    data = [];
    current = '';

    dictKeys = keys(dict);

    for i = 1:length(compressed)
        current = [current, compressed(i)];

        for k = 1:length(dictKeys)
            key = dictKeys{k};
            if strcmp(current, key)
                data(end+1) = dict(key);
                current = '';
                break; % Found a match, break out of loop
            end
        end
    end
end