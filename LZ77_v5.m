%%  Source Coding - Final Project
%   - LZ77 Algorithm -
%   Tommaso Martini (108 15 80)

%   v5.0 improvements:
%   - floating bits management. I use a non multiple of 8 number os bits
%   and code in bytes only at last step

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   BUGS & "TO-FIX"'s
%   - ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;

%% Initialization

% Program parameters
verbose_mode = false;
generate_files = false;

% Algorithm parameters
search_window_length = 500;
coding_window_length = 1000;

% Implementation parameters
file_name_input = './cantrbry/cp.html';
file_name_input = './big_files/2';
dictionary_output = 'lz77_dictionary_output_2.txt';
file_name_output = 'lz77_output_2.txt';
M = 256;  % alphabet cardinality

%% Pick a file from the filesystem
stored_file_ID = fopen(file_name_input);
seq = fread(stored_file_ID, Inf, '*uint8');
seq = seq';
msg_length = length(seq);
fclose(stored_file_ID);

%% Encoder

offset_size = ceil(log2(search_window_length));
length_size = ceil(log2(coding_window_length));
symbol_size = ceil(log2(M));

dict_index = 2; % index to span the dictionary

dictionary(1, :) = [0, 0, double(seq(1))];  % first triple

search_index = 1;   % first element of the search window
coding_index = 2;   % first element of the coding window

end_of_file = false;    % when the end of the file is reached this is turned to TRUE

while ~end_of_file
    
    % Pattern matching
    
    pattern = seq(coding_index : min(msg_length - 1, coding_index + coding_window_length - 1));
    m = length(pattern);
    
    if m == 0   % last symbol of the sequence: encode it with a single symbol (so doing I'm not wasting resources)
        % New row in the dictionary
        dictionary(dict_index, :) = [0, 0, double(seq(end))];
        dict_index = dict_index + 1;
    else
        conti = true;
        tmp_pattern = pattern(1);
        
        longest_match = 0;
        match_position = coding_index;
        
        while conti
            search_string = [seq(search_index : coding_index - 1), tmp_pattern(1 : end - 1)];   % the search_string is made by search_window and coding_window
            match_positions = strfind(search_string, tmp_pattern);
            if ~isempty(match_positions)    % some matches found: search again
                longest_match = length(tmp_pattern);
                match_position = search_index + match_positions(1) - 1;
                if length(tmp_pattern) < length(pattern)
                    tmp_pattern = pattern(1 : length(tmp_pattern) + 1);
                else    % I have used the whole pattern
                    conti = false;
                end
            else    % no matches found: stop the cycle
                conti = false;
            end
        end
        
        offset = coding_index - match_position;
        
        % New row in the dictionary
        dictionary(dict_index, :) = [offset, longest_match, double(seq(coding_index + longest_match))];
        dict_index = dict_index + 1;
    end
    
    % Update indeces to scan the file
    coding_index = coding_index + longest_match + 1;
    search_index = max(coding_index - search_window_length, 1); % you cannot start from the char before the first one
    
    if coding_index > length(seq)
        end_of_file = true;
    end
    
    if verbose_mode
        clc;
        fprintf('Dictionary generation progress: %d%%', round(coding_index * 100 / msg_length));
    end
end

%% Dictionary compression

bit_cod_sequence = [];
for dict_row = 1 : size(dictionary, 1)
    
    curr_off = dictionary(dict_row, 1);
    curr_off_bit = de2bi(curr_off);
    curr_off_bit = [curr_off_bit, zeros(1, offset_size - length(curr_off_bit))];
    
    curr_len = dictionary(dict_row, 2);
    curr_len_bit = de2bi(curr_len);
    curr_len_bit = [curr_len_bit, zeros(1, length_size - length(curr_len_bit))];
    
    curr_sym = dictionary(dict_row, 3);
    curr_sym_bit = de2bi(curr_sym);
    curr_sym_bit = [curr_sym_bit, zeros(1, symbol_size - length(curr_sym_bit))];
    
    bit_cod_sequence = [bit_cod_sequence, curr_off_bit, curr_len_bit, curr_sym_bit];
    
    if verbose_mode
        clc;
        disp('Dictionary generation progress: 100%');
        fprintf('Compression data progress: %d%% \n', round(dict_row * 100 / size(dictionary, 1)));
    end
end

num_final_zeros = mod(length(bit_cod_sequence), 8);
if num_final_zeros > 0
    bit_cod_sequence = [bit_cod_sequence, zeros(1, 8 - num_final_zeros)];
end
cod_sequence = [];
while ~isempty(bit_cod_sequence)
    piece = bit_cod_sequence(1 : 8);
    bit_cod_sequence = bit_cod_sequence(9 : end);
    curr_byte = bi2de(piece);
    cod_sequence = [cod_sequence, curr_byte];
end

% Embedding information about the size in bits of offset and length: 3
% bytes
symbol_size_byte = uint8(symbol_size);
offset_size_byte = uint8(offset_size);
length_size_byte = uint8(length_size);
cod_sequence = [symbol_size_byte, offset_size_byte, length_size_byte, cod_sequence];

if verbose_mode
    disp('Coding complete!');
end

if generate_files
    cod_file_ID = fopen(dictionary_output, 'w');
    fwrite(cod_file_ID, cod_sequence);
    fclose(cod_file_ID);
    
    %% Pick up an encoded file
    coded_file_ID = fopen(dictionary_output);
    coded_dictionary = fread(coded_file_ID, Inf, '*uint8');
    coded_dictionary = coded_dictionary';
    fclose(coded_file_ID);
else
    coded_dictionary = cod_sequence;
end

% Extract information about the sizes
symbol_size = double(coded_dictionary(1));
offset_size = double(coded_dictionary(2));
length_size = double(coded_dictionary(3));
coded_dictionary = coded_dictionary(4 : end);

% Bring the dictionary in bit form
bit_coded_dictionary = [];
while ~isempty(coded_dictionary)
    piece = coded_dictionary(1);
    coded_dictionary = coded_dictionary(2 : end);
    curr_bits = de2bi(piece);
    curr_bits = [curr_bits, zeros(1, 8 - length(curr_bits))];
    bit_coded_dictionary = [bit_coded_dictionary, curr_bits];
end
bit_coded_dictionary = double(bit_coded_dictionary);

% decoded_dictionary = zeros(length(coded_dictionary) / (offset_size + length_size + symbol_size), 3);
decoded_dictionary = [];
dictionary_row_index = 1;

while length(bit_coded_dictionary) >= symbol_size + offset_size + length_size
    
    offset_bits = bit_coded_dictionary(1 : offset_size);
    curr_off = bi2de(offset_bits);
    bit_coded_dictionary = bit_coded_dictionary(offset_size + 1 : end);
    
    length_bits = bit_coded_dictionary(1 : length_size);
    curr_len = bi2de(length_bits);
    bit_coded_dictionary = bit_coded_dictionary(length_size + 1 : end);
    
    symbol_bits = bit_coded_dictionary(1 : symbol_size);
    curr_sym = bi2de(symbol_bits);
    bit_coded_dictionary = bit_coded_dictionary(symbol_size + 1 : end);
    
    decoded_dictionary(dictionary_row_index, :) = [curr_off, curr_len, curr_sym];
    dictionary_row_index = dictionary_row_index + 1;
    
    if verbose_mode
        clc;
        disp('Dictionary generation progress: 100%');
        disp('Compression data progress: 100%');
        disp('Coding complete!');
        fprintf('Expanding dictionary progress: %d%%', round(dict_row * 100 / size(decoded_dictionary, 1)));
    end
end

%% Check whether the dictionaries are equal
if ~isequal(dictionary, decoded_dictionary)
    disp('Error! Dictionaries not equal!');
end

%% Decoder
dict_length = size(decoded_dictionary, 1);

decoded_sequence = -1;
decoded_seq_index = 1;  % first EMPTY position

for dictionary_row_index = 1 : dict_length
    
    dictionary_row = decoded_dictionary(dictionary_row_index, :);
    
    offset = dictionary_row(1);
    prefix_length = dictionary_row(2);
    last_symbol = dictionary_row(3);
    
    prefix_from = decoded_seq_index - offset;
    
    for prefix_symbol = 1 : prefix_length
        decoded_sequence(decoded_seq_index) = decoded_sequence(prefix_from + prefix_symbol - 1);
        decoded_seq_index = decoded_seq_index + 1;
    end
    
    decoded_sequence(decoded_seq_index) = last_symbol;
    decoded_seq_index = decoded_seq_index + 1;
    
    if verbose_mode
        clc;
        disp('Dictionary generation progress: 100%');
        disp('Compression data progress: 100%');
        disp('Coding complete!');
        disp('Expanding dictionary progress: 100%');
        fprintf('Decoding data progress: %d%%', round(dictionary_row_index * 100 / dict_length));
    end
end

decoded_sequence = uint8(decoded_sequence);

if verbose_mode
    clc;
    disp('Dictionary generation progress: 100%');
    disp('Compression data progress: 100%');
    disp('Coding complete!');
    disp('Expanding dictionary progress: 100%');
    disp('Decoding data progress: 100%');
    disp('Decoding complete!');
end

if generate_files
    %% Store the decoded file
    dec_file_ID = fopen(file_name_output, 'w');
    fprintf(dec_file_ID, char(decoded_sequence));
    fclose(dec_file_ID);
end

%% Error check

if verbose_mode
    disp('Error check started...');
end

err_counter = 0;
error_locations = [];
for i = 1 : msg_length
    if seq(i) ~= decoded_sequence(i)
        err_counter = err_counter + 1;
        error_locations = [error_locations, i];
    end
end

if err_counter > 0
    fprintf('%d errors found! \n', err_counter);
else
    disp('No errors found!');
end

%% Performances analysis
comp_msg_size = length(cod_sequence)
original_msg_size = msg_length
compression_ratio = comp_msg_size * 100 / original_msg_size
if verbose_mode
    fprintf('Compression: %d %%', compression_ratio);
end
