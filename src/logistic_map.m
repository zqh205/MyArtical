function f = logistic_map(num)
    num = 65536;
%测试代码
    % clc;close all;clear;
    % num = 256;
    % smallImage_side_len = 16;

    mu = 3.6;       %参数
    x = 0.4;        %初值
    n1 = 1000;      %门限值
    f = zeros(1,num);

    for i = 1:n1 + num
        x = mu.*x.*(1-x);         %Logistic映射
    end
    
    f = f(n1+1:n1+num);
    %对混沌序列标准归一化
    normalized_f = (f - min(f)) / (max(f) - min(f));
    
    % % 测试代码
    % figure;
    % plot(1:num, normalized_f(1,:), 'k');
    % xlabel('N')
    % % ylabel('x')
    % histogram(normalized_f, 50);
end