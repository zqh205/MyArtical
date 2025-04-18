function restored_image = inverse_zigzag(reordered_image)
    [rows, cols] = size(reordered_image);

    % 获取Zigzag扫描的逆索引
    zigzag_order = zigzag_index(rows, cols);

    % 还原图像
    restored_image = zeros(rows, cols);
    restored_image(zigzag_order) = reshape(reordered_image, [], 1);

    % % 显示结果
    % figure;
    % subplot(1, 3, 3), imshow(restored_image, []), title('Restored Image');
end

function index = zigzag_index(rows, cols)
    % 创建zigzag扫描的索引矩阵
    index = zeros(rows, cols);
    index_order = 1;
    for s = 1:(rows + cols - 1)
        if mod(s, 2) == 0
            for i = max(1, s-cols+1):min(s, rows)
                j = s - i + 1;
                index(i, j) = index_order;
                index_order = index_order + 1;
            end
        else
            for j = max(1, s-rows+1):min(s, cols)
                i = s - j + 1;
                index(i, j) = index_order;
                index_order = index_order + 1;
            end
        end
    end
    % 将2D索引矩阵转换为1D索引向量
    [~, index] = sort(index(:));
end