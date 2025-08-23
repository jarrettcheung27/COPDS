function BCH_Codec(config_path, infile_path, outfile_path, mode)
%     Description: Input k * 1 index messages, encode those messages into n * 1 index codeword.
%     Params:
%         coding_config: coding configuration filepath
%         infile_path: input messages filepath
%         outfile_path: output codewords filepath
%         mode: encode or decode
%     Data Strcuture:
%            input: text file, k_bch rows bit of length n0, bits are
%            not seperated.
%            output: text file, n_bch rows bit of length n0, bits are
%            seperated by comma in each line.

    % if args.length < 4
    %    error('Not enough input arguments. Usage: BCH_Codec(config_path, infile_path, outfile_path, mode)');
    % end
    % Load configuration from config.json
    coding_config = jsondecode(fileread(config_path));

    n1_bch = coding_config.ECC.inner1.n;
    k1_bch = coding_config.ECC.inner1.k;
    n2_bch = coding_config.ECC.inner2.n;
    k2_bch = coding_config.ECC.inner2.k;
    
    % Load messages from infile_path, read by line and append each line to
    % the matrix
    % Read lines from infile_path
    msgs = readlines(infile_path);
    
    % Convert the one-dimensional array of strings into a two-dimensional matrix
    % Each row corresponds to a binary vector extracted from the string
    msgs = cell2mat(cellfun(@(x) int8(x) - '0', msgs, 'UniformOutput', false));
    n_0 = length(msgs(1,:));
    
    if mode == "encode" % Encode mode
        % recognize layer
        if length(msgs(:,1)) == k1_bch
            k_bch = k1_bch;
            n_bch = n1_bch;
        elseif length(msgs(:,1)) == k2_bch
            k_bch = k2_bch;
            n_bch = n2_bch;
        else
            error('Input messages length does not match any BCH configuration.');
        end
        % Generate BCH generator polynomial
        bchEncoder = comm.BCHEncoder(n_bch, k_bch);
        temp = zeros(n_bch,n_0);
        for i = 1 : n_0
            % Encode the messages
            temp(:,i) = bchEncoder(msgs(:,i));
        end
        % Write codewords to the outfile_path(.txt)
        writematrix(temp,outfile_path);
  
    else % Decode mode
        % recognize layer
        if length(msgs(:,1)) == n1_bch
            k_bch = k1_bch;
            n_bch = n1_bch;
        elseif length(msgs(:,1)) == n2_bch
            k_bch = k2_bch;
            n_bch = n2_bch;
        else
            error('Input messages length does not match any BCH configuration.');
        end
        % Generate BCH generator polynomial
        bchEncoder = comm.BCDecoder(n_bch, k_bch);
        temp = zeros(k_bch+1,n_0);
        for i = 1 : n_0
            % Encode the messages
            temp(:,i) = bchEncoder(msgs(:,i));
        end
        % Write codewords to the outfile_path(.txt)
        writematrix(temp,outfile_path);
    end
end