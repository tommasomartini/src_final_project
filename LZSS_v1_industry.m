%%  Source Coding - Final Project
%   - LZSS Algorithm -
%   Tommaso Martini (108 15 80)

%   v1.0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   BUGS & "TO-FIX"'s
%   + ho un ultima simbolo nel decoded_dict che interpreto come indice (lorem2)
%   + la decoded sequence viene di lunghezza molto diversa (hodor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;

%% Initialization

% Program parameters
verbose_mode = true;

windows = [10, 100, 1000, 10000, 100000];
perform = zeros(1, length(windows));

for win = 1 : length(windows)
    
    % Algorithm parameters
    search_window_length = windows(win);
    coding_window_length = windows(win);
    
    % Implementation parameters
    file_name_input = './cantrbry/sum';
    % file_name_input = 'test_Hodor.txt';
    dictionary_output = 'dictionary_output.txt';
    file_name_output = 'alice29_output.txt';
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
        
        pattern = seq(coding_index : min(msg_length, coding_index + coding_window_length - 1));
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
            
            if offset == 0  % no matches: encode a symbol
                dictionary(dict_index, :) = [0, 0, double(seq(coding_index))];
            else    % match: encode a pair
                dictionary(dict_index, :) = [1, offset, longest_match];
            end
            
            % New row in the dictionary
            dict_index = dict_index + 1;
        end
        
        % Update indeces to scan the file
        if longest_match == 0   %no matches: go ahead of 1
            coding_index = coding_index + 1;
        else    % match: go ahead of longest_match
            coding_index = coding_index + longest_match;
        end
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
        
        if mod(dict_row - 1, 8) == 0
            indeces_octave = dictionary(dict_row : min(dict_row + 8 - 1, size(dictionary, 1)), 1);
            indeces_octave = indeces_octave';
            if size(dictionary, 1) - dict_row + 1 < 8
                indeces_octave = [indeces_octave, zeros(1, 8 - (size(dictionary, 1) - dict_row + 1))];
            end
            indeces_byte = uint8(bi2de(indeces_octave));
            cod_sequence = [cod_sequence, indeces_byte(1)];
        end
        
        if dictionary(dict_row, 1) == 0     % encode a symbol
            symbol64 = uint64(dictionary(dict_row, 3));
            symbol8 = typecast(symbol64, 'uint8');
            symbol_bytes = symbol8(1 : symbol_size);
            cod_sequence = [cod_sequence, symbol_bytes];
        else    % encode a pair
            offset64 = uint64(dictionary(dict_row, 2));
            offset8 = typecast(offset64, 'uint8');
            offset_bytes = offset8(1 : offset_size);
            
            length64 = uint64(dictionary(dict_row, 3));
            length8 = typecast(length64, 'uint8');
            length_bytes = length8(1 : length_size);
            
            cod_sequence = [cod_sequence, [offset_bytes, length_bytes]];
        end
        
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
    coded_dictionary_length = length(coded_dictionary);
    coded_dictionary = coded_dictionary';
    fclose(coded_file_ID);
    
    % decoded_dictionary = zeros(length(coded_dictionary) / (offset_size + length_size + symbol_size), 3);
    decoded_dictionary = [];
    dictionary_row_index = 1;
    
    current_index = 9;
    indeces_sequence = [];
    
    iinnd = 0;
    
    while ~isempty(coded_dictionary)
        
        if length(coded_dictionary) == 1
            disp('sfsf')
        end
        if current_index > 8
            octave_index = coded_dictionary(1);
            coded_dictionary = coded_dictionary(2 : end);
            indeces_sequence = de2bi(octave_index);
            indeces_sequence = [indeces_sequence, zeros(1, 8 - length(indeces_sequence))];
            current_index = 1;
            
            iinnd = iinnd + 1;
        end
        
        if isempty(coded_dictionary)
            disp('sfsf')
        end
        
        if indeces_sequence(current_index) == 0   % decode a symbol
            symbol_bytes = coded_dictionary(1 : symbol_size);
            coded_dictionary = coded_dictionary(symbol_size + 1 : end);
            symbol8_dec = [symbol_bytes, zeros(1, 8 - symbol_size)];
            symbol64_dec = typecast(symbol8_dec, 'uint64');
            symbol_dec = double(symbol64_dec);
            decoded_dictionary(dictionary_row_index, :) = [0, 0, symbol_dec];
            dictionary_row_index = dictionary_row_index + 1;
        else    % decode a pair
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
            
            decoded_dictionary(dictionary_row_index, :) = [1, offset_dec, length_dec];
            dictionary_row_index = dictionary_row_index + 1;
        end
        current_index = current_index + 1;
        
        if verbose_mode
            clc;
            disp('Dictionary generation progress: 100%');
            disp('Compression data progress: 100%');
            disp('Coding complete!');
            fprintf('Expanding dictionary progress: %d%%', round((coded_dictionary_length -length(coded_dictionary)) * 100 / coded_dictionary_length));
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
        
        info_index = dictionary_row(1);
        
        if info_index == 0
            last_symbol = dictionary_row(3);
            decoded_sequence(decoded_seq_index) = last_symbol;
            decoded_seq_index = decoded_seq_index + 1;
        else
            offset = dictionary_row(2);
            prefix_length = dictionary_row(3);
            prefix_from = decoded_seq_index - offset;
            
            for prefix_symbol = 1 : prefix_length
                decoded_sequence(decoded_seq_index) = decoded_sequence(prefix_from + prefix_symbol - 1);
                decoded_seq_index = decoded_seq_index + 1;
            end
        end
        
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
    
    % error_locations
    
    %% Performances analysis
    comp_msg_size = length(cod_sequence)
    original_msg_size = msg_length
    compression_ratio = round(comp_msg_size * 100 / original_msg_size);
    if verbose_mode
        fprintf('Compression: %d %%', compression_ratio);
    end
    perform(win) = compression_ratio;
end