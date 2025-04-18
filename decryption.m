    clc; close all; clear;
% function decryimage = decryption(encry_image)
    %% 导入密文和密钥
    loaddata = load('.\key\encry_image.mat');
    encry_image = loaddata.encry_image;
    loaddata = load('.\key\init_chaos.mat');
    xyzw0 = loaddata.init_chaos;
    loaddata = load('.\key\max.mat');
    max_mat = loaddata.max_mat;
    loaddata = load('.\key\min.mat');
    min_mat = loaddata.min_mat;
    loaddata = load('.\key\subimagesize.mat');
    subimagesize = loaddata.subimagesize;
    loaddata = load('.\key\encry_image_compress.mat');
    encry_image_compress = loaddata.encry_image_compress;


    [M, N] = size(encry_image);         %获取原始图像尺寸
    side = 2 ^ subimagesize;            %子图像边长
    size = side * side;                 %子图像大小
    num = M * N / size;                 %子图像的数量num
    side_num = sqrt(num);
    % % 加噪声
    % encry_image_compress = salt_pepper(encry_image_compress);
    % 剪切攻击
    encry_image_compress = shelter(encry_image_compress, 1/6);

    encry_image((side_num - 1) * side + 1 : side_num * side, ...
            (side_num - 1) * side + 1 : side_num * side) = encry_image_compress;
   
    
    %% 混沌序列生成
    x0 = xyzw0(1); y0 = xyzw0(2); z0 =xyzw0(3); w0 = xyzw0(4);
    [X, Y, Z, W] = fourD_chaotic_system(x0, y0, z0, w0, M*N);

    %处理混沌序列X
    seq = X(1:num);                     % 子图像排序序列，由于生成索引序列
    X= reshape(X, M, N);
    sub_X = cell(side_num, side_num);   % 
    
    for i = 1:side_num
        for j = 1:side_num
            sub_X{i,j} = X((i - 1) * side + 1 : i * side, ...
                (j - 1) * side + 1 : j * side);
        end
    end

    % 处理混沌序列Y
    mask = Y(1 : side * side);
    mask = reshape(mask, side, side);   % 起始随机相位掩模2

    Y = reshape(Y, M, N);
    sub_Y = cell(side_num, side_num);   
    
    for i = 1:side_num
        for j = 1:side_num
            sub_Y{i,j} = Y((i - 1) * side + 1 : i * side, ...
                (j - 1) * side + 1 : j * side);
        end
    end

    %处理混沌序列Z
    Z = reshape(Z, M, N);% Z作为DNA运算的密钥图像
    sub_Z = cell(side_num, side_num);
    
    for i = 1:side_num
        for j = 1:side_num
            sub_Z{i,j} = Z((i - 1) * side + 1 : i * side, ...
                (j - 1) * side + 1 : j * side);
        end
    end

    %处理混沌序列W，作为DNA解码规则
    W = reshape(W, M, N);
    sub_W = cell(side_num, side_num);
    
    for i = 1:side_num
        for j = 1:side_num
            sub_W{i,j} = W((i - 1) * side + 1 : i * side, ...
                (j - 1) * side + 1 : j * side);
        end
    end

    %% DNAcoding解密
    real_part = real(encry_image);
    imaginary_part = imag(encry_image);

    encry_image_DRPE = DNAcoding_and_compute_rev(real_part*255, round(Z*255),round(Y*7)+1, ...
           round(W*7)+1);
    encry_image_DRPE = im2double(encry_image_DRPE);
    
    %归一化还原
    for i = 1:side_num
        for j = 1:side_num
            temp = encry_image_DRPE((i - 1) * side + 1 : i * side, ...
                (j - 1) * side + 1 : j * side);
            orig = temp * (max_mat(i,j) - min_mat(i,j)) + min_mat(i,j);
            encry_image_DRPE((i - 1) * side + 1 : i * side, ...
                (j - 1) * side + 1 : j * side) = orig;
        end
    end

    encry_image_DRPE = encry_image_DRPE+1i*imaginary_part;      %实部加虚部
    %待DRPE解密图像分块
    de_DNA = cell(side_num, side_num);
    for i = 1:side_num
        for j = 1:side_num
            de_DNA{i,j} = encry_image_DRPE((i - 1) * side + 1 : i * side, ...
                (j - 1) * side + 1 : j * side);
        end
    end

    %掩模2分块
    mask2 = cell(side_num *side_num, 1);
    for i = 1:side_num
        for j = 1:side_num
            mask2{(i-1)*side_num+j} = real_part((i - 1) * side + 1 : i * side, ...
                (j - 1) * side + 1 : j * side);
        end
    end

    for i = 1:side_num*side_num-1
        mask2(side_num*side_num-i+1) = mask2(side_num*side_num-i);
    end
    mask2{1,1} = mask;

    mask2 = reshape(mask2, side_num, side_num)';
    

    %% DRPE解密
    plain = cell(side_num, side_num);
    for i = 1:side_num
        for j = 1:side_num
            plain{i,j} = DRPE_inverse(de_DNA{i,j}, sub_X{i,j}, mask2{i,j});
        end
    end
    %% 排序逆过程
    [seq, idx] = sort(seq);
    inv_idx = zeros(numel(idx), 1); % 使用 numel 确保长度匹配
    inv_idx(idx) = 1:numel(idx);
    plain = plain.';
    temp_sub_image = reshape(plain, side_num * side_num, 1);
    
    sorted_sub_image = cell(side_num * side_num, 1);
    for i = 1:num 
        sorted_sub_image{i, 1} = temp_sub_image{inv_idx(i), 1};
    end
    
    sub_image = reshape(sorted_sub_image, side_num, side_num)';

    %zigzag逆过程
    for i = 1:side_num
         for j = 1:side_num
             sub_image{i, j} = inverse_zigzag(sub_image{i, j});
         end
    end
    
    plain_image = zeros(M, N);  %密文图像
    for i = 1:side_num
        for j = 1:side_num
            plain_image((i - 1) * side + 1 : i * side, ...
                (j - 1) * side + 1 : j * side) = sub_image{i,j}; 
        end
    end

    imshow(real(plain_image));
    % 保存解密图像
    imwrite(real(plain_image), '.\output\Cameraman_decry.tif');
    

    

