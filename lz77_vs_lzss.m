%%  Source Coding - Final Project
%   - LZ77 vs LZSS -
%   Tommaso Martini (108 15 80)

%   This program compares the performances of LZ77 and LZSS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   BUGS & "TO-FIX"'s
%   - ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;

series_number = 2;

M = 256;  % alphabet cardinality

file_numbers = 1 : 7;

window_size = 10000;

data_matrix = zeros(length(file_numbers), 2);

for file_num = file_numbers
    
    switch file_num
        case 1
            max_win_span = 10000;
        case 2
            max_win_span = 13000;
        case 3
            max_win_span = 10000;
        case 4
            max_win_span = 13000;
        case 5
            max_win_span = 24000;
        case 6
            max_win_span = 11000;
        case 7
            max_win_span = 38000;
        otherwise
            max_win_span = 20000;
    end
    
    base_string_lz77 = './matrici risultati/Risultati LZ77Big/serie';
    mat_name_input = strcat(base_string_lz77, num2str(series_number), '/lz77_res_series', num2str(series_number), '_', num2str(file_num), '.mat');
    load(mat_name_input)
    index = find(windows_span == window_size);
    data_matrix(file_num, 1) = performances(index);
    
    base_string_lzss = './matrici risultati/Risultati LZSSBig/serie';
    mat_name_input = strcat(base_string_lzss, num2str(series_number), '/lzss_res_series', num2str(series_number), '_', num2str(file_num), '.mat');
    load(mat_name_input)
    index = find(windows_span == window_size);
    data_matrix(file_num, 2) = performances(index);
end

figure;
bar(file_numbers, data_matrix);
hold on;
plot(0 : length(file_numbers) + 1, 100 * ones(length(file_numbers) + 2))