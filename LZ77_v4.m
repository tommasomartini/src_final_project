%%  Source Coding - Final Project
%   - LZ77 Algorithm -
%   Tommaso Martini (108 15 80)

%   v4.0 improvements:
%   - pattern matching is executed in a pessimistic way, starting from the
%   smallest match. This improves time performances.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   BUGS & "TO-FIX"'s
%   - Fare un pattern matching a bisezione? Provo L, poi 1, poi L - 1, poi
%   2 e cosi' via...
%   + ottengo una riga di dizionario [0, 1, 111] al quarto posto con
%   test_Hodor.txt
%   + conl'esempio di Hodor non metto le r
%   + gestione ultimo simbolo. Qui l'ho corretta, ma non nlle versioni 2 e
%   3!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;

%% Initialization

% Program parameters
verbose_mode = false;

% Algorithm parameters
search_window_length = 10000;
coding_window_length = 2000;

% Implementation parameters
file_name_input = './cantrbry/cp.html';
file_name_input = './big_files/2';
file_name_input = 'sam_test.txt';
% file_name_input = 'test_Hodor.txt';
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

% How many bytes do I need to encode each parameter?
offset_size = ceil(ceil(log2(search_window_length)) / 8);
length_size = ceil(ceil(log2(search_window_length + coding_window_length)) / 8);
symbol_size = ceil(ceil(log2(M)) / 8);

cod_file_ID = fopen(dictionary_output, 'w');

cod_sequence = [];
for dict_row = 1 : size(dictionary, 1)
    
    offset64 = uint64(dictionary(dict_row, 1));
    offset8 = typecast(offset64, 'uint8');
    offset_bytes = offset8(1 : offset_size);
    
    length64 = uint64(dictionary(dict_row, 2));
    length8 = typecast(length64, 'uint8');
    length_bytes = length8(1 : length_size);
    
    symbol64 = uint64(dictionary(dict_row, 3));
    symbol8 = typecast(symbol64, 'uint8');
    symbol_bytes = symbol8(1 : symbol_size);
    
    cod_sequence = [cod_sequence, [offset_bytes, length_bytes, symbol_bytes]];
    
    if verbose_mode
        clc;
        disp('Dictionary generation progress: 100%');
        fprintf('Compression data progress: %d%% \n', round(dict_row * 100 / size(dictionary, 1)));
    end
end

fwrite(cod_file_ID, cod_sequence);
fclose(cod_file_ID);

if verbose_mode
    disp('Coding complete!');
end

%% Pick up an encoded file
coded_file_ID = fopen(dictionary_output);
coded_dictionary = fread(coded_file_ID, Inf, '*uint8');
coded_dictionary = coded_dictionary';
fclose(coded_file_ID);

decoded_dictionary = zeros(length(coded_dictionary) / (offset_size + length_size + symbol_size), 3);
for dict_row = 1 : size(decoded_dictionary, 1)
    
    offset_bytes = coded_dictionary(1 : offset_size);
    coded_dictionary = coded_dictionary(offset_size + 1 : end);
    offset8_dec = [offset_bytes, zeros(1, 8 - offset_size)];
    offset64_dec = typecast(offset8_dec, 'uint64');
    offset_dec = double(offset64_dec);
    
    length_bytes = coded_dictionary(1 : length_size);
    coded_dictionary = coded_dictionary(length_size + 1 : end);
    length8_dec = [length_bytes, zeros(1, 8 - length_size)];
    length64_dec = typecast(length8_dec, 'uint64');
    length_dec = double(length64_dec);
    
    symbol_bytes = coded_dictionary(1 : symbol_size);
    coded_dictionary = coded_dictionary(symbol_size + 1 : end);
    symbol8_dec = [symbol_bytes, zeros(1, 8 - symbol_size)];
    symbol64_dec = typecast(symbol8_dec, 'uint64');
    symbol_dec = double(symbol64_dec);
    
    decoded_dictionary(dict_row, :) = [offset_dec, length_dec, symbol_dec];
    
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

%% Store the decoded file
dec_file_ID = fopen(file_name_output, 'w');
fprintf(dec_file_ID, char(decoded_sequence));
fclose(dec_file_ID);

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
byte_per_triplet = offset_size + length_size + symbol_size;
comp_msg_size = byte_per_triplet * dict_length
original_msg_size = msg_length
compression_ratio = round(comp_msg_size * 100 / original_msg_size);
if verbose_mode
    fprintf('Compression: %d %%', compression_ratio);
end
