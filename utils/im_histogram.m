function im_histogram(encrypted_image) 
    gray_hist = imhist(encrypted_image);   
    
    % 绘制直方图  
    figure;  
 
    bar(gray_hist, 'k');   
    ylim([0,600]);

end