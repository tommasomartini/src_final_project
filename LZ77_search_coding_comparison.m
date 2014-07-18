%%  Source Coding - Final Project
%   - LZ77 Algorithm -
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

M = 256;  % alphabet cardinality

file_numbers = 1 : 7;

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
    search_windows_span = 1000 : 1000 : max_win_span;
    coding_windows_span = 1000 : 1000 : max_win_span;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     search_windows_span = 100 : 500 : 1000;
%     coding_windows_span = 100 : 500 : 1000;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    performances = zeros(length(search_windows_span), length(coding_windows_span));
    
    for search_win = 1 : length(search_windows_span)
        for coding_win = 1 : length(coding_windows_span)
            
            search_window_length = search_windows_span(search_win);
            coding_window_length = coding_win;
            
            seq = theseq;
            
            %% Encoder
            
            dict_index = 2; % index to span the dictionary
            
            dictionary = [];
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
            end
            
            %% Dictionary compression
            
            % How many bytes do I need to encode each parameter?
            offset_size = ceil(ceil(log2(search_window_length)) / 8);
            length_size = ceil(ceil(log2(search_window_length + coding_window_length)) / 8);
            symbol_size = ceil(ceil(log2(M)) / 8);
            
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
                
            end
            
            %% Performances analysis
            byte_per_triplet = offset_size + length_size + symbol_size;
            comp_msg_size = byte_per_triplet * size(dictionary, 1);
            original_msg_size = msg_length;
            compression_ratio = round(comp_msg_size * 100 / original_msg_size);
            
            performances(search_win, coding_win) = compression_ratio;
        end
    end
    
    res_file_name = strcat('performances_results_comp_search_coding_', num2str(file_num));
    save(res_file_name, 'search_windows_span', 'coding_windows_span', 'performances');
    
%     figure;
%     mesh(search_windows_span, coding_windows_span, performances);
end



