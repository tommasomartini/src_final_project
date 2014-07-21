%%  Source Coding - Final Project
%   Tommaso Martini (108 15 80)

%   This program is used to compare the compression varying both search and
%   coding windows

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   BUGS & "TO-FIX"'s
%   - ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;

lz77_prefix_name = 'lz77_comp_s1_';
lzss_prefix_name = 'lzss_comp_s1_';

windows_step = 1000;
windows_start = 1000;

M = 256;  % alphabet cardinality

file_numbers = [1, 2, 3, 4, 6];

optimistic_files = [1];

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
    
    %% Pick a file from the filesystem
    file_name_input = strcat('./big_files/', num2str(file_num));
    stored_file_ID = fopen(file_name_input);
    theseq = fread(stored_file_ID, Inf, '*uint8');
    theseq = theseq';
    msg_length = length(theseq);
    fclose(stored_file_ID);
    
    % lengths of the windows
    search_windows_span = windows_start : windows_step : max_win_span;
    coding_windows_span = windows_start : windows_step : max_win_span;
    lz77_performances = zeros(length(search_windows_span), length(coding_windows_span));
    lzss_performances = zeros(length(search_windows_span), length(coding_windows_span));
    
    %     Double cycle on windows
    for search_win_index = 1 : length(search_windows_span)
        for coding_win_index = 1 : length(coding_windows_span)
            
            search_window_length = search_windows_span(search_win_index);
            coding_window_length = coding_windows_span(coding_win_index);
            
            fprintf('File: %d, srch: %d, cod: %d \n', file_num, search_window_length, coding_window_length);
            
            offset_size = ceil(log2(search_window_length));
            length_size = ceil(log2(coding_window_length));
            symbol_size = ceil(log2(M));
            
            seq = theseq;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% LZ77 Encoder
            
            dict_index = 2; % index to span the dictionary
            
            lz77_dictionary = [];
            lz77_dictionary(1, :) = [0, 0, double(seq(1))];  % first triple
            
            search_index = 1;   % first element of the search window
            coding_index = 2;   % first element of the coding window
            
            if sum(ismember(file_num, optimistic_files)) > 0    % optimistic approach
                end_of_file = false;    % when the end of the file is reached this is turned to TRUE
                while ~end_of_file
                    pattern = seq(coding_index : min(msg_length - 1, coding_index + coding_window_length - 1));
                    m = length(pattern);
                    if m == 0   % last symbol of the sequence: encode it with a single symbol (so doing I'm not wasting resources)
                        lz77_dictionary(dict_index, :) = [0, 0, double(seq(end))];
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
                        
                        % New row in the dictionary
                        lz77_dictionary(dict_index, :) = [offset, match_length, double(seq(coding_index + match_length))];
                        dict_index = dict_index + 1;
                    end
                    
                    % Update indeces to scan the file
                    % Update indeces to scan the file
                    coding_index = coding_index + match_length + 1;
                    search_index = max(coding_index - search_window_length, 1); % you cannot start from the char before the first one
                    
                    if coding_index > length(seq)
                        end_of_file = true;
                    end
                end
                
            else
                end_of_file = false;    % when the end of the file is reached this is turned to TRUE
                while ~end_of_file
                    pattern = seq(coding_index : min(msg_length - 1, coding_index + coding_window_length - 1));
                    m = length(pattern);
                    if m == 0   % last symbol of the sequence: encode it with a single symbol (so doing I'm not wasting resources)
                        % New row in the dictionary
                        lz77_dictionary(dict_index, :) = [0, 0, double(seq(end))];
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
                        lz77_dictionary(dict_index, :) = [offset, longest_match, double(seq(coding_index + longest_match))];
                        dict_index = dict_index + 1;
                    end
                    
                    % Update indeces to scan the file
                    coding_index = coding_index + longest_match + 1;
                    search_index = max(coding_index - search_window_length, 1); % you cannot start from the char before the first one
                    
                    if coding_index > length(seq)
                        end_of_file = true;
                    end
                end
            end
            
            %% LZ77 Dictionary compression
            
            bit_cod_sequence = [];
            for dict_row = 1 : size(lz77_dictionary, 1)
                
                curr_off = lz77_dictionary(dict_row, 1);
                curr_off_bit = de2bi(curr_off);
                curr_off_bit = [curr_off_bit, zeros(1, offset_size - length(curr_off_bit))];
                
                curr_len = lz77_dictionary(dict_row, 2);
                curr_len_bit = de2bi(curr_len);
                curr_len_bit = [curr_len_bit, zeros(1, length_size - length(curr_len_bit))];
                
                curr_sym = lz77_dictionary(dict_row, 3);
                curr_sym_bit = de2bi(curr_sym);
                curr_sym_bit = [curr_sym_bit, zeros(1, symbol_size - length(curr_sym_bit))];
                
                bit_cod_sequence = [bit_cod_sequence, curr_off_bit, curr_len_bit, curr_sym_bit];
            end
            
            num_final_zeros = mod(length(bit_cod_sequence), 8);
            if num_final_zeros > 0
                bit_cod_sequence = [bit_cod_sequence, zeros(1, 8 - num_final_zeros)];
            end
            lz77_cod_sequence = [];
            while ~isempty(bit_cod_sequence)
                piece = bit_cod_sequence(1 : 8);
                bit_cod_sequence = bit_cod_sequence(9 : end);
                curr_byte = bi2de(piece);
                lz77_cod_sequence = [lz77_cod_sequence, curr_byte];
            end
            
            %             symbol_size_byte = uint8(symbol_size);
            %             offset_size_byte = uint8(offset_size);
            %             length_size_byte = uint8(length_size);
            %             lz77_cod_sequence = [symbol_size_byte, offset_size_byte, length_size_byte, lz77_cod_sequence];
            
            %% LZ77 Performances analysis
            comp_msg_size = length(lz77_cod_sequence) + 3;
            lz77_compression_ratio = comp_msg_size * 100 / msg_length;
            lz77_performances(search_win_index, coding_win_index) = lz77_compression_ratio;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% LZSS Encoder
            
            pair = offset_size + length_size + 1;   % number of bits to encode a pair
            single = symbol_size + 1;   % number of bits to encode a single symbol
            symbol_thr = ceil(pair / single);  % minimum number of symbol I can encode as a pair. If the match length is shorter than symbol_thr I encode them as single symbols
            
            dict_index = 2; % index to span the dictionary
            lzss_dictionary(1, :) = [0, 0, double(seq(1))];  % first triple
            
            search_index = 1;   % first element of the search window
            coding_index = 2;   % first element of the coding window
            
            if sum(ismember(file_num, optimistic_files)) > 0    % optimistic approach
                end_of_file = false;    % when the end of the file is reached this is turned to TRUE
                while ~end_of_file
                    
                    pattern = seq(coding_index : min(msg_length, coding_index + coding_window_length - 1));
                    m = length(pattern);
                    
                    if m == 0   % last symbol of the sequence: encode it with a single symbol (so doing I'm not wasting resources)
                        % New row in the dictionary
                        lzss_dictionary(dict_index, :) = [0, 0, double(seq(end))];
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
                            lzss_dictionary(dict_index, :) = [0, 0, double(seq(coding_index))];
                        else    % match: encode a pair
                            lzss_dictionary(dict_index, :) = [1, offset, match_length];
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
                    
                    if coding_index > length(seq)
                        end_of_file = true;
                    end
                end
            else    % pessimistic approach
                end_of_file = false;    % when the end of the file is reached this is turned to TRUE
                while ~end_of_file
                    
                    pattern = seq(coding_index : min(msg_length, coding_index + coding_window_length - 1));
                    m = length(pattern);
                    
                    if m == 0   % last symbol of the sequence: encode it with a single symbol (so doing I'm not wasting resources)
                        % New row in the dictionary
                        lzss_dictionary(dict_index, :) = [0, 0, double(seq(end))];
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
                            lzss_dictionary(dict_index, :) = [0, 0, double(seq(coding_index))];
                            dict_index = dict_index + 1;
                        elseif longest_match < symbol_thr   % match, but it is convenient to encode symbols independently
                            for i = 0 : longest_match - 1
                                lzss_dictionary(dict_index, :) = [0, 0, double(seq(coding_index + i))];
                                dict_index = dict_index + 1;
                            end
                        else    % good match: encode a pair
                            lzss_dictionary(dict_index, :) = [1, offset, longest_match];
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
                    
                    if coding_index > length(seq)
                        end_of_file = true;
                    end
                end
            end
            
            %% LZSS Dictionary compression
            
            bit_cod_sequence = [];
            for dict_row = 1 : size(lzss_dictionary, 1)
                
                if mod(dict_row - 1, 8) == 0    % every 8 dict rows I have to insert a byte with flags
                    indeces_octave = lzss_dictionary(dict_row : min(dict_row + 8 - 1, size(lzss_dictionary, 1)), 1);
                    indeces_octave = indeces_octave';
                    if size(lzss_dictionary, 1) - dict_row + 1 < 8   % if I'm left with less than 8 rows...
                        indeces_octave = [indeces_octave, zeros(1, 8 - (size(lzss_dictionary, 1) - dict_row + 1))]; % put zeros as last flags
                    end
                    bit_cod_sequence = [bit_cod_sequence, indeces_octave];
                end
                
                if lzss_dictionary(dict_row, 1) == 0     % encode a symbol
                    curr_sym = lzss_dictionary(dict_row, 3);
                    curr_sym_bit = de2bi(curr_sym);
                    curr_sym_bit = [curr_sym_bit, zeros(1, symbol_size - length(curr_sym_bit))];
                    bit_cod_sequence = [bit_cod_sequence, curr_sym_bit];
                else    % encode a pair
                    curr_off = lzss_dictionary(dict_row, 2);
                    curr_off_bit = de2bi(curr_off);
                    curr_off_bit = [curr_off_bit, zeros(1, offset_size - length(curr_off_bit))];
                    
                    curr_len = lzss_dictionary(dict_row, 3);
                    curr_len_bit = de2bi(curr_len);
                    curr_len_bit = [curr_len_bit, zeros(1, length_size - length(curr_len_bit))];
                    bit_cod_sequence = [bit_cod_sequence, curr_off_bit, curr_len_bit];
                end
            end
            
            num_final_zeros = mod(length(bit_cod_sequence), 8);
            if num_final_zeros > 0
                bit_cod_sequence = [bit_cod_sequence, zeros(1, 8 - num_final_zeros)];
            end
            lzss_cod_sequence = [];
            while ~isempty(bit_cod_sequence)
                piece = bit_cod_sequence(1 : 8);
                bit_cod_sequence = bit_cod_sequence(9 : end);
                curr_byte = bi2de(piece);
                lzss_cod_sequence = [lzss_cod_sequence, curr_byte];
            end
            
            % Embedding information about the size in bits of offset and length: 3
            % bytes
            symbol_size_byte = uint8(symbol_size);
            offset_size_byte = uint8(offset_size);
            length_size_byte = uint8(length_size);
            lzss_cod_sequence = [symbol_size_byte, offset_size_byte, length_size_byte, lzss_cod_sequence];
            
            %% LZSS Performances analysis
            comp_msg_size = length(lzss_cod_sequence);
            lzss_compression_ratio = comp_msg_size * 100 / msg_length;
            lzss_performances(search_win_index, coding_win_index) = lzss_compression_ratio;
        end
    end
    
    lz77_res_file_name = strcat(lz77_prefix_name, num2str(file_num));
    save(lz77_res_file_name, 'search_windows_span', 'coding_windows_span', 'lz77_performances');
    
    lzss_res_file_name = strcat(lzss_prefix_name, num2str(file_num));
    save(lzss_res_file_name, 'search_windows_span', 'coding_windows_span', 'lzss_performances');
end



