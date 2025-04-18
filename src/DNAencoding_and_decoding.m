function res = DNAencoding_and_decoding (image, encode_rule, decode_rule)
    %DNA编码规则
    rule = cell(8, 1);
    rule{1} = ['A', 'C', 'G', 'T'];
    rule{2} = ['A', 'G', 'C', 'T'];
    rule{3} = ['T', 'G', 'C', 'A'];
    rule{4} = ['T', 'C', 'G', 'A'];
    rule{5} = ['C', 'A', 'T', 'G'];
    rule{6} = ['C', 'T', 'A', 'G'];
    rule{7} = ['G', 'A', 'T', 'C'];
    rule{8} = ['G', 'T', 'A', 'C'];

    % 初始化变量
    [rows, cols] = size(image);
    encoded_image = strings(rows, cols, 4); % 用于存储编码后的DNA序列
    

    % 编码image
    for i = 1:rows
        for j = 1:cols
            % 将像素值转换为8位二进制
            binary_value = dec2bin(image(i, j), 8);
            
            % 使用encode_rule对二进制进行编码
            for k = 1:4
                % 每2位构成一个DNA碱基，使用encode_rule编码
                pair = binary_value(2*k-1:2*k);
                index = bin2dec(pair) + 1; % 1-based index for MATLAB
                encoded_image(i, j, k) = rule{encode_rule(i,j), 1}(index);
            end
        end
    end

    % 初始化解码结果
    decoded_image = zeros(rows, cols, 'uint8'); % 解码后的图像

    % 根据decode_rule进行解码
    for j = 1:cols
        for i = 1:rows
            dna_seq = encoded_image(i, j, :);
            dna_seq = dna_seq(:)'; % 压缩为行向量

            binary_value = '';
            for k = 1:4
                % 找到DNA碱基在decode_rule中的位置，base要转为char类型
                base = char(dna_seq(k));
                index = find(rule{decode_rule(i,j)} == base, 1) - 1; % 找到碱基对应的二进制索引
                binary_value = strcat(binary_value, dec2bin(index, 2));
            end
            
            % 将二进制字符串转换为整数
            decoded_image(i, j) = bin2dec(binary_value);
        end
    end

    % 返回解码后的图像
    res = decoded_image;
end