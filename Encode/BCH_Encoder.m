function BCH_Encoder(config_path, infile_path, outfile_path)
%     Description: Input k * 1 index messages, encode those messages into n * 1 index codeword.
%     Params:
%         coding_config: coding configuration filepath
%         infile_path: input messages filepath
%         outfile_path: output codewords filepath
%         layer: which layer to encode
%     Data Strcuture:
%            input: text file, k_bch rows bit of length n0, bits are
%            not seperated.
%            output: text file, n_bch rows bit of length n0, bits are
%            seperated by comma in each line.

    % Load configuration from config.json
    coding_config = jsondecode(fileread(config_path));
    n_bch = coding_config.ECC.inner2.n;
    
    % Load messages from infile_path, read by line and append each line to
    % the matrix
    % Read lines from infile_path
    msgs = readlines(infile_path);
    
    % Convert the one-dimensional array of strings into a two-dimensional matrix
    % Each row corresponds to a binary vector extracted from the string
    msgs = cell2mat(cellfun(@(x) int8(x) - '0', msgs, 'UniformOutput', false));
    k_bch = length(msgs(:,1)); 
    n_0 = length(msgs(1,:));
    % Generate BCH generator polynomial
    bchEncoder = comm.BCHEncoder(n_bch, k_bch);
    temp = zeros(n_bch,n_0);
    for i = 1 : n_0
        % Encode the messages
        temp(:,i) = bchEncoder(msgs(:,i));
    end
    % Write codewords to the outfile_path(.txt)
    writematrix(temp,outfile_path);

end