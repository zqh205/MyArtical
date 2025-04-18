function res = shelter(img, x)
    % x为剪切掉的大小，例如1/4
    h = length(img(1,:));  
    w = length(img(:, 1));

    temp1 = zeros(fix(w * x), h);  
    temp2 = ones(w - fix(w * x) , h);

    res = [temp1; temp2];  
    res = res .* img;
end