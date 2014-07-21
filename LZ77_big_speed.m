%%  Source Coding - Final Project
%   - LZ77 Algorithm -
%   Tommaso Martini (108 15 80)

%   'Big' version 2: speed up lucky cases

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   BUGS & "TO-FIX"'s
%   - ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;

prefix_name = 'lz77_res_series2_';

M = 256;  % alphabet cardinality

file_numbers = 1 : 2;

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
        coding_window_length = 1000;
        
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
                dictionary(dict_index, :) = [offset, match_length, double(seq(coding_index + match_length))];
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
        compression_ratio = comp_msg_size * 100 / original_msg_size;
        
        performances(win) = compression_ratio;
    end
    
    res_file_name = strcat(prefix_name, num2str(file_num));
    save(res_file_name, 'windows_span', 'performances');
end


