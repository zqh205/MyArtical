function [x0, y0, z0, w0] = init_val_gene(hashdata, randInt)
    h1 = 0.0;
    h2 = 0.0;
    h3 = 0.0;
    h4 = 0.0;
    
    for i = 1:32
        if i <= 8
            h1 = bitxor(h1, hashdata(i, 1));
        end
        if (i > 9) && (i <= 16)
            h2 = bitxor(h2, hashdata(i, 1));
        end
        if i > 17 && i <= 24
            h3 = bitxor(h3, hashdata(i, 1));
        end
        if i > 25 && i <= 32
            h4 = bitxor(h4, hashdata(i, 1));
        end
    end
    
    
    %引入随机数影响
    h1 = bitxor(h1, randInt(1,1));
    h2 = bitxor(h2, randInt(1,2));
    h3 = bitxor(h3, randInt(1,3));
    h4 = bitxor(h4, randInt(1,4));

    h1 = double(h1) / 255.0;
    h2 = double(h2) / 255.0;
    h3 = double(h3) / 255.0;
    h4 = double(h4) / 255.0;
    
    x0 = mod((h1 + h2) * 10 ^ 14, 256) / 255.0;
    y0 = mod((h2 + h3) * 10 ^ 14, 256) / 255.0;
    z0 = mod((h3 + h4) * 10 ^ 14, 256) / 255.0;
    w0 = mod((h1 + h2 + h3 + h4) * 10 ^ 14, 256) /255.0;
    
    % disp([x0, y0, z0, w0]);
end