function [C] = DRPE(I, A1, B1)
    %加密过程
    % A = im2double(I);
    A = I;

    % figure;
    % imshow(A);
    % title("DRPE内的原图");

    [m, n] = size(I);
    A11 = exp(1i*2*pi*A1);   %掩模1
    A111 = A.*A11;
    
    B = fft2(A111);     %傅里叶变换
    B11 = exp(1i*2*pi*B1);   % 掩模2
    B111 = B.*B11;
    C = ifft2(B111);    %傅里叶逆变换
    
    C1 = abs(C);      % 加密图像
end
