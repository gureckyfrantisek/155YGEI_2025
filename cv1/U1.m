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

        % Update the matrix
        Y(row : row+7, col : col+7) = Yt_q;
        Cb(row : row+7, col : col+7) = Cbt_q;
        Cr(row : row+7, col : col+7) = Crt_q;
    end
end

%%%% JPEG Decompression
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