function grayChangeRate = calculateGrayChangeRate(image)
    % 确保输入图像是灰度图
    if size(image, 3) == 3
        image = rgb2gray(image);
    end
    
    % 获取图像的大小
    [rows, cols] = size(image);
    
    % 初始化灰度变化合计
    totalGrayChange = 0;
    validCount = 0;
    
    % 遍历图像的每一个像素
    for i = 1:rows
        for j = 1:cols
            % 当前像素的灰度值
            currentPixel = double(image(i, j));
            
            % 检查右邻居像素
            if j < cols
                rightPixel = double(image(i, j + 1));
                totalGrayChange = totalGrayChange + abs(currentPixel - rightPixel);
                validCount = validCount + 1;
            end
            
            % 检查下邻居像素
            if i < rows
                bottomPixel = double(image(i + 1, j));
                totalGrayChange = totalGrayChange + abs(currentPixel - bottomPixel);
                validCount = validCount + 1;
            end
        end
    end
    
    % 计算平均灰度变化
    grayChangeRate = totalGrayChange / validCount;
end