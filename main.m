clc; close all; clear;
addpath('src', 'utils', 'image');

image = imread("boat_512.tif");
image = imresize(image, [256, 256]);
% image = rgb2gray(image);
% imshow(image);
image = im2double(image);     %转换为双精度浮点数

%生成sha-256
hashdata = hash(image, 'SHA-256');
% disp(hashdata);
% disp(class(hashdata));

%生成4个随机数,作为外部密钥，用于生成混沌初值
randInt = randi([0, 255], 1, 4);

%生成初始值
[x0, y0, z0, w0] = init_val_gene(hashdata, randInt);

%图像基本数据
[M, N] = size(image);       %获取原始图像尺寸
subimagesize = 6;         %随机数控制子图像大小
side = 2 ^ subimagesize;             %子图像边长
size = side * side;         %子图像大小
num = M * N / size;         %子图像的数量num
side_num = sqrt(num);

%% 生成混沌序列
%定义4个混沌序列
X = zeros(M * N,1);
Y = zeros(M * N,1);
Z = zeros(M * N,1);
W = zeros(M * N,1);

% 生成4个混沌序列
[X, Y, Z, W] = fourD_chaotic_system(x0, y0, z0, w0, M*N);
% disp([X, Y, Z, W]);

%%
%将明文图像分割为子图像
sub_image = cell(side_num,side_num);
for i = 1:side_num
    for j = 1:side_num
        sub_image{i,j} = image((i - 1) * side + 1 : i * side, ...
            (j - 1) * side + 1 : j * side);
    end
end

% disp(sub_image{1,1});
%% 处理混沌序列X,Y,Z,W
% 处理混沌序列X，控制子图像重排，作为随机相位掩模1
seq = X(1:num);                     % 用于子图像重排，生成索引序列
X = reshape(X, M, N);
sub_X = cell(side_num, side_num);   % 保存分割的X

for i = 1:side_num
    for j = 1:side_num
        sub_X{i,j} = X((i - 1) * side + 1 : i * side, ...
            (j - 1) * side + 1 : j * side);
    end
end

% 处理混沌序列Y，作为DNA编码规则，初始随机相位掩模2
mask = Y(1 : side * side);
mask = reshape(mask, side, side);   % 初始随机相位掩模2

Y = reshape(Y, M, N);
sub_Y = cell(side_num, side_num);   

for i = 1:side_num
    for j = 1:side_num
        sub_Y{i,j} = Y((i - 1) * side + 1 : i * side, ...
            (j - 1) * side + 1 : j * side);
    end
end

%处理混沌序列Z，作为DNA_xor密钥图像，用同样的编码规则编码
Z = reshape(Z, M, N);
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

%% 加密过程
%对每一幅子图像进行zigzag扫描
for i = 1:side_num
     for j = 1:side_num
         sub_image{i, j} = zigzag(sub_image{i, j});
     end
end
% disp(sub_image{1,1});

%对子图像进行重新排序
[seq, idx] = sort(seq);

sub_image = sub_image';
temp_sub_image = reshape(sub_image, side_num * side_num, 1);

sorted_sub_image = cell(side_num * side_num, 1);
for i = 1:num 
    sorted_sub_image{i, 1} = temp_sub_image{idx(i), 1};
end

sub_image = reshape(sorted_sub_image, side_num, side_num)';
% save('sub_image',"sub_image");
%% DRPE-DNA扩散
encry_images = cell(side_num, side_num);    %保存子密文图像
imaginary_part = cell(side_num, side_num);  %把DRPE后虚部保存下来，最后与实部经过DNA加密后的结果拼接为复值
imaginary_part = cellfun(@(x) zeros(side, side), imaginary_part,'UniformOutput', false);
min_mat = zeros(side_num,side_num);         %最小值矩阵，保存子图像DRPE后实部的最小值，用于解密
man_mat = zeros(side_num,side_num);         %最大值矩阵，保存子图像DRPE后实部的最大值，用于解密
for i = 1 : side_num
    for j = 1 : side_num
       temp = DRPE(sub_image{i, j}, sub_X{i, j}, mask);
       % if i==1&&j==1
       %     save('DRPE_res','temp');
       % end
       real_part = real(temp);
       imaginary_part{i,j} = imag(temp);            %取出虚部
       min_mat(i,j) = min(real_part(:));            %实部最小值
       max_mat(i,j) = max(real_part(:));            %实部最小值
       temp_norm = (real_part-min_mat(i,j))/(max_mat(i,j)-min_mat(i,j));    %实部归一化
       % if i==1&&j==1
       %     save('DRPE_res_norm','temp_norm');
       % end
       encry_images{i,j} = DNAcoding_and_compute(round(temp_norm*255),round(sub_Z{i,j}*255),round(sub_Y{i,j}*7)+1, ...
           round(sub_W{i,j}*7)+1);                  %归一化的实部用于DNA编码运算
       mask = im2double(encry_images{i,j});         %DNA编码结果作为下一次DRPE的掩模2
       encry_images{i,j} = mask;                    %DNA编码结果也作为密文图像的实部
       % mask2 = DNAcoding_image(abs(encry_images{i, j}), sub_Z{i, j}, sub_W{i, j});
       % imshow(mask2);
    end
end

%虚部实部加起来
for i = 1:side_num
    for j = 1:side_num
        encry_images{i,j} = double(encry_images{i,j})+ 1i.*double(imaginary_part{i,j});
    end
end

% 子密文图像合并为密文图像
encry_image = zeros(M, N);  %密文图像
for i = 1:side_num
    for j = 1:side_num
        encry_image((i - 1) * side + 1 : i * side, ...
            (j - 1) * side + 1 : j * side) = encry_images{i,j}; 
    end
end
encry_image_compress = encry_images{side_num,side_num};

init_chaos = [x0,y0,z0,w0];
save('.\key\encry_image.mat', 'encry_image');
save('.\key\init_chaos', 'init_chaos');
save('.\key\min', "min_mat");
save('.\key\max', "max_mat");
save(".\key\subimagesize.mat", 'subimagesize');
save(".\key\encry_image_compress.mat", 'encry_image_compress');

%% 密文分析,密文要用uint8的0~255整数表示
encrypimg = real(encry_image);
encrypimg = uint8(encrypimg * 255);
%保存图像的矢量图
imwrite(image, '.\output\Cameraman_plain.tif');
imwrite(encrypimg, '.\output\Cameraman_encry.tif');
%显示图像
figure;
imshow(encrypimg);
%直方图
im_histogram(image);
im_histogram(encrypimg);
im_histogram(uint8(real(encry_image_compress)*255))
%熵
entro = entropy(encrypimg);
disp(entro);
%相关系数
correlation = im_correlation(image*255,4000);
disp(correlation);
correlation = im_correlation(encrypimg,4000);
disp(correlation);
%像素数变化率NPCR和统一平均变化强度UACI
[NPCR, UACI] = NPCR_UACI(image*255, encrypimg);
disp([NPCR, UACI]);



        