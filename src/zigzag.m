function reordered_image = zigzag(image)
    [rows, cols] = size(image);

    % 获取Zigzag扫描的索引
    zigzag_order = zigzag_index(rows, cols);

    % 按Zigzag扫描顺序获取像素，并按行排列
    flattened_image = image(zigzag_order);
    reordered_image = reshape(flattened_image, rows, cols);
    
    % % 显示结果
    % figure;
    % subplot(1, 3, 1), imshow(image, []), title('Original Image');
    % subplot(1, 3, 2), imshow(reordered_image, []), title('Zigzag Reordered Image');
end

function index = zigzag_index(rows, cols)
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
    [~, index] = sort(index(:));
end