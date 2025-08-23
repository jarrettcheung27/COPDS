function BCH_Encoder(config_path, infile_path, outfile_path, layer)
%     Description: Input k * 1 index messages, encode those messages into n * 1 index codeword.
%     Params:
%         coding_config: coding configuration filepath
%         infile_path: input messages filepath
%         outfile_path: output codewords filepath
%         layer: which layer to encode

    % Load configuration from config.json
    coding_config = jsondecode(fileread(config_path));
    n_bch = coding_config['ECC'][layer]['n']
    k_bch = coding_config['ECC'][layer]['k']
    m = coding_config['ECC']['outer']['n']
    % Generate BCH generator polynomial
    bchEncoder = comm.BCHEncoder(double(n_bch), double(k_bch));
    
    % Load messages from infile_path(.txt)
    msgs = readmatrix(infile_path);
    % Encode messages
    encoded_msgs = zeros(m,n_bch);
    for i = 1 : m
        % Encode the message
        encoded_msgs(i,:) = (bchEncoder(msgs(i,:)'))';
    end
    % Write codewords to the outfile_path(.txt)
    writematrix(encoded_msgs, outfile_path);
end