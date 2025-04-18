function res = DNAcoding_and_compute_rev (encry_image, image_xor, encode_rule, decode_rule)
    %image、image_xor为0~255的整数，encode_rule,decode_rule为1~8的整数
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
    [rows, cols] = size(encry_image);
    encoded_image = strings(rows, cols, 4); % 用于存储编码后的DNA序列
    encoded_xor_image = strings(rows, cols, 4); % 用于存储image_xor编码后的DNA序列
    xor_result = strings(rows, cols, 4); % 存储DNA异或结果

    % 编码image
    for i = 1:rows
        for j = 1:cols
            % 将像素值转换为8位二进制
            binary_value = dec2bin(encry_image(i, j), 8);
            
            % 使用encode_rule对二进制进行编码
            for k = 1:4
                % 每2位构成一个DNA碱基，使用encode_rule编码
                pair = binary_value(2*k-1:2*k);
                index = bin2dec(pair) + 1; % 1-based index for MATLAB
                encoded_image(i, j, k) = rule{decode_rule(i,j), 1}(index);
            end
        end
    end

    % 编码image_xor
    for i = 1:rows
        for j = 1:cols
            % 将像素值转换为8位二进制
            binary_value = dec2bin(image_xor(i, j), 8);
            
            % 使用encode_rule对二进制进行编码
            for k = 1:4
                % 每2位构成一个DNA碱基，使用encode_rule编码
                pair = binary_value(2*k-1:2*k);
                index = bin2dec(pair) + 1; % 1-based index for MATLAB
                encoded_xor_image(i, j, k) = rule{encode_rule(i,j)}(index);
            end
        end
    end

    % 执行DNA级的异或操作
    for i = 1:rows
        for j = 1:cols
            for k = 1:4
                base1 = char(encoded_image(i, j, k));
                base2 = char(encoded_xor_image(i, j, k));
                % DNA碱基异或通过字符串比较，映射同位置不同碱基
                xor_base = dna_xor(base1, base2);
                xor_result(i, j, k) = xor_base;
            end
        end
    end

    % 初始化解码结果
    decoded_image = zeros(rows, cols, 'uint8'); % 解码后的图像

    % 根据decode_rule进行解码
    for j = 1:cols
        for i = 1:rows
            dna_seq = xor_result(i, j, :);
            dna_seq = dna_seq(:)'; % 压缩为行向量

            binary_value = '';
            for k = 1:4
                % 找到DNA碱基在decode_rule中的位置，base要转为char类型
                base = char(dna_seq(k));
                index = find(rule{encode_rule(i,j)} == base, 1) - 1; % 找到碱基对应的二进制索引
                binary_value = strcat(binary_value, dec2bin(index, 2));
            end
            
            % 将二进制字符串转换为整数
            decoded_image(i, j) = bin2dec(binary_value);
        end
    end

    % 返回解码后的图像
    res = decoded_image;
end

% function xor_base = dna_xor(encode_rule, base1, base2)
%     % 定义DNA碱基之间的异或关系
%     rule = cell(8, 1);
%     rule{1} = ['A', 'C', 'G', 'T'];
%     rule{2} = ['A', 'G', 'C', 'T'];
%     rule{3} = ['T', 'G', 'C', 'A'];
%     rule{4} = ['T', 'C', 'G', 'A'];
%     rule{5} = ['C', 'A', 'T', 'G'];
%     rule{6} = ['C', 'T', 'A', 'G'];
%     rule{7} = ['G', 'A', 'T', 'C'];
%     rule{8} = ['G', 'T', 'A', 'C'];
% 
%     %为了符合find()函数的输入，转换为字符
%     base1 = char(base1);
%     base2 = char(base2);
% 
%     idx1 = find(rule{encode_rule} == base1);
%     idx2 = find(rule{encode_rule} == base2);
%     xor_index = bitxor(idx1-1, idx2-1) + 1; % MATLAB中索引是1-based
%     xor_base = rule{encode_rule}(xor_index);
% end