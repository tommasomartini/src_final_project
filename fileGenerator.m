close all;
clear all;
clc;

% sequence = 65 * ones(1, 10000);
% sequence = char(sequence);

% alphabet = 'abcdefghijklmnopqrstuvwxyz';
% sequence = repmat(alphabet, [1, 500]);
% length(sequence)

% sequence = randi([0, 255], [1, 10000]);
% sequence = char(sequence);

% sequence = 65 * ones(1, 100);
% sequence = [sequence, 66 * ones(1, 100)];
% sequence = [sequence, 67 * ones(1, 100)];
% sequence = [sequence, 68 * ones(1, 100)];
% sequence = repmat(sequence, [1, 3]);

% sequence = randi([0, 255], [1, 1000]);
% sequence = repmat(sequence, [1, 10]);
% sequence = char(sequence);

cod_file_ID = fopen('22', 'w');
fwrite(cod_file_ID, sequence);
fclose(cod_file_ID);