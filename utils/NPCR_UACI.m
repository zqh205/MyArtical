function [npcr, uaci] = NPCR_UACI(original_image, encrypted_image)
    % 计算 NPCR 和 UACI 的函数
    % 输入:
    %   original_image - 原始图像 (二维矩阵)
    %   encrypted_image - 加密后的图像 (二维矩阵)
    % 输出:
    %   npcr - 像素数变化率
    %   uaci - 统一平均变化强度

    %转double
    original_image = double(original_image);
    encrypted_image = double(encrypted_image);

    % 获取图像的尺寸
    [H, W] = size(original_image);

    % 初始化 NPCR 和 UACI 的计算变量
    npcr_sum = 0;
    uaci_sum = 0;

    % 遍历图像的每个像素
    for i = 1:H
        for j = 1:W
            % 计算像素值是否发生变化
            if original_image(i, j) ~= encrypted_image(i, j)
                npcr_sum = npcr_sum + 1;
            end

            % 计算像素值的差值
            uaci_sum = uaci_sum + abs(original_image(i, j) - encrypted_image(i, j));
        end
    end

    % 计算 NPCR
    npcr = ((npcr_sum) / (H * W));

    % 计算 UACI
    uaci = uaci_sum / (H * W) /255;

    % 返回 NPCR 和 UACI
end