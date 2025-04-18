function [F] = DRPE_inverse(I, A1, B1)
    %解密过程
    A = im2double(I);
    [m, n] = size(I);
  
    %%
    C = A;
    D = fft2(C);    %傅里叶变换
    D1 = D.*exp(-1i*2*pi*B1);   %乘以掩模2的共轭
    D11 = ifft2(D1);        %傅里叶逆变换
    D111 = D11.*exp(-1i*2*pi*A1);   %乘以掩模1的共轭
    F = D111;       %解密图像

    % figure;
    % imshow(F);
    % title("解密函数内的解密图像")

end
