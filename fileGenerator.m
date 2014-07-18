close all;
clear all;
clc;

% sequence = 65 * ones(1, 10000);
% sequence = char(sequence);

% alphabet = 'abcdefghijklmnopqrstuvwxyz';
% sequence = repmat(alphabet, [1, 500]);
% length(sequence)

sequence = randi([0, 255], [1, 10000]);
sequence = char(sequence);

cod_file_ID = fopen('33', 'w');
fwrite(cod_file_ID, sequence);
fclose(cod_file_ID);