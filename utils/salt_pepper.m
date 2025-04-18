function noisyImage = salt_pepper(img)
% 添加椒盐噪声  
    noisyImage = imnoise(img, 'salt & pepper', 0.1); % 2% 的噪声  
end