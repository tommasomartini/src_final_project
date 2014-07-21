%%  Source Coding - Final Project
%   - LZSS Algorithm -
%   Tommaso Martini (108 15 80)

%   Search window changes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   BUGS & "TO-FIX"'s
%   - ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;

prefix_name = 'lzss_search_series1_';

M = 256;  % alphabet cardinality

coding_window_length = 10000;

file_numbers = 1 : 7;
file_numbers = [1, 2, 3, 4, 6];

for file_num = file_numbers
    
    file_num
    
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
    
    %% Pick a file from the filesystem
    
    file_name_input = strcat('./big_files/', num2str(file_num));
    stored_file_ID = fopen(file_name_input);
    theseq = fread(stored_file_ID, Inf, '*uint8');
    theseq = theseq';
    msg_length = length(theseq);
    fclose(stored_file_ID);
    
    % lengths of the windows
    windows_span = 1000 : 1000 : max_win_span;
    performances = zeros(1, length(windows_span));
    
    for win = 1 : length(windows_span)
        
        search_window_length = windows_span(win);
        
        seq = theseq;
        
        % Expressed in bits
        offset_size = ceil(log2(search_window_length - 1));
        length_size = ceil(log2(coding_window_length - 1));
        symbol_size = ceil(log2(M));
        
        pair = offset_size + length_size + 1;   % number of bits to encode a pair
        single = symbol_size + 1;   % number of bits to encode a single symbol
        symbol_thr = ceil(pair / single);  % minimum number of symbol I can encode as a pair. If the match length is shorter than symbol_thr I encode them as single symbols
        
        dict_index = 2; % index to span the dictionary
        dictionary(1, :) = [0, 0, double(seq(1))];  % first triple
        
        search_index = 1;   % first element of the search window
        coding_index = 2;   % first element of the coding window
        
        end_of_file = false;    % when the end of the file is reached this is turned to TRUE
        
        while ~end_of_file
            
            if file_num == 1   % optimistic approach

                pattern = seq(coding_index : min(msg_length, coding_index + coding_window_length - 1));
                m = length(pattern);
                
                if m == 0   % last symbol of the sequence: encode it with a single symbol (so doing I'm not wasting resources)
                    % New row in the dictionary
                    dictionary(dict_index, :) = [0, 0, double(seq(end))];
                    dict_index = dict_index + 1;
                else
                    match_length = length(pattern);
                    search_string = [seq(search_index : coding_index - 1), pattern(1 : end - 1)];   % the search_string is made by search_window and coding_window
                    match_positions = strfind(search_string, pattern);
                    
                    while isempty(match_positions)
                        pattern = pattern(1 : end - 1);
                        search_string = search_string(1 : end - 1);
                        match_length = length(pattern);
                        if isempty(pattern)    % there are no matches
                            match_positions = coding_index; % in this way I know there have been no matches
                        else
                            match_positions = strfind(search_string, pattern);
                        end
                    end
                    
                    match_position = search_index + match_positions(1) - 1;
                    offset = coding_index - match_position;
                    
                    if offset == 0  % no matches: encode a symbol
                        dictionary(dict_index, :) = [0, 0, double(seq(coding_index))];
                    else    % match: encode a pair
                        dictionary(dict_index, :) = [1, offset, match_length];
                    end
                    
                    % New row in the dictionary
                    dict_index = dict_index + 1;
                end
                
                % Update indeces to scan the file
                if match_length == 0   %no matches: go ahead of 1
                    coding_index = coding_index + 1;
                else    % match: go ahead of longest_match
                    coding_index = coding_index + match_length;
                end
                search_index = max(coding_index - search_window_length, 1); % you cannot start from the char before the first one
                
            else    % pessimistic approach
                
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
                        dict_index = dict_index + 1;
                    elseif longest_match < symbol_thr   % match, but it is convenient to encode symbols independently
                        for i = 0 : longest_match - 1
                            dictionary(dict_index, :) = [0, 0, double(seq(coding_index + i))];
                            dict_index = dict_index + 1;
                        end
                    else    % good match: encode a pair
                        dictionary(dict_index, :) = [1, offset, longest_match];
                        dict_index = dict_index + 1;
                    end
                end
                
                % Update indeces to scan the file
                if longest_match == 0   %no matches: go ahead of 1
                    coding_index = coding_index + 1;
                else    % match: go ahead of longest_match
                    coding_index = coding_index + longest_match;
                end
                search_index = max(coding_index - search_window_length, 1); % you cannot start from the char before the first one
            end
            if coding_index > length(seq)
                end_of_file = true;
            end
        end
        
        %% Dictionary compression
        
        bit_cod_sequence = [];
        for dict_row = 1 : size(dictionary, 1)
            
            if mod(dict_row - 1, 8) == 0    % every 8 dict rows I have to insert a byte with flags
                indeces_octave = dictionary(dict_row : min(dict_row + 8 - 1, size(dictionary, 1)), 1);
                indeces_octave = indeces_octave';
                if size(dictionary, 1) - dict_row + 1 < 8   % if I'm left with less than 8 rows...
                    indeces_octave = [indeces_octave, zeros(1, 8 - (size(dictionary, 1) - dict_row + 1))]; % put zeros as last flags
                end
                bit_cod_sequence = [bit_cod_sequence, indeces_octave];
            end
            
            if dictionary(dict_row, 1) == 0     % encode a symbol
                curr_sym = dictionary(dict_row, 3);
                curr_sym_bit = de2bi(curr_sym);
                curr_sym_bit = [curr_sym_bit, zeros(1, symbol_size - length(curr_sym_bit))];
                bit_cod_sequence = [bit_cod_sequence, curr_sym_bit];
            else    % encode a pair
                curr_off = dictionary(dict_row, 2);
                curr_off_bit = de2bi(curr_off);
                curr_off_bit = [curr_off_bit, zeros(1, offset_size - length(curr_off_bit))];
                
                curr_len = dictionary(dict_row, 3);
                curr_len_bit = de2bi(curr_len);
                curr_len_bit = [curr_len_bit, zeros(1, length_size - length(curr_len_bit))];
                bit_cod_sequence = [bit_cod_sequence, curr_off_bit, curr_len_bit];
            end
        end
        
        num_final_zeros = mod(length(bit_cod_sequence), 8);
        bit_cod_sequence = [bit_cod_sequence, zeros(1, 8 - num_final_zeros)];
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
        
        %% Performances analysis
        comp_msg_size = length(cod_sequence);
        original_msg_size = msg_length;
        compression_ratio = comp_msg_size * 100 / original_msg_size;
        performances(win) = compression_ratio;
    end
    
    res_file_name = strcat(prefix_name, num2str(file_num));
    save(res_file_name, 'windows_span', 'performances');
end