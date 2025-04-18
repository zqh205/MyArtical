function [X,Y,Z,W] = fourD_chaotic_system(x0, y0, z0, w0, size_of_image)
    % ẋ = ax + byz
    % ẏ = cy + dxz
    % ż = exy + kz + mxw
    % ẇ = ny

    %参数包括a,b,c,d,e,k(-26,-7),m,n,状态参数包括x,y,z,w       
    a = 8;
    b = -1;
    c = -40;
    d = 1;
    e = 2;
    % k = randi([-26, -7]);
    k = -14;
    m = 1;
    n = -2;
    
    % x0 = 0.5118;    
    % y0 = 0.6103;    
    % z0 = 0.8662;    
    % w0 = 0.3741;
    
    % x0 = 0.4;
    % y0 = 0.13;
    % z0 = 0.2;
    % w0 = 0.8;
    % size_of_image = 3000;

    total_iterations = 4000 + size_of_image;  % 总的迭代步数，前 4000 要舍弃
    % T = total_iterations * 1;  % 根据迭代和步长设定总时间
    dt = 1;   % 时间步长

    % 延迟时间步数
    tau = 2;  % 延迟时间，以秒计
    delay_steps = floor(tau / dt);  % 将延迟时间转换为时间步数

    % 初始化，创建存储 x, y, z 的数组
    x = zeros(1, total_iterations);  % 用于存储 x 值的数组
    y = zeros(1, total_iterations);  % 用于存储 y 值的数组
    z = zeros(1, total_iterations);  % 用于存储 z 值的数组
    w = zeros(1, total_iterations);  % 用于存储 w 值的数组

    % 设置初始值
    x(1) = x0;  % x 的初始条件
    y(1) = y0;  % y 的初始条件
    z(1) = z0;  % z 的初始条件
    w(1) = w0;  % w 的初始条件

    % 数值积分主循环
    for t = 1:total_iterations-1
        % ẋ = ax + byz
        % ẏ = cy + dxz
        % ż = exy + kz + mxw
        % ẇ = ny
        % 洛伦兹系统的状态方程计算
        dx = a * x(t) + b * y(t) *z(t);  % x 的微分方程
        dy = c * y(t) + d * x(t) * z(t);  % y 的微分方程
        dz = e * x(t) * y(t) + k * z(t) + m * x(t) * w(t);  % z 的微分方程
        dw = n * y(t);

        % x(t + 1) = x(t) + dx * dt;  
        % y(t + 1) = y(t) + dy * dt;  
        % z(t + 1) = z(t) + dz * dt;  
        % w(t + 1) = w(t) + dw * dt;

        % 使用欧拉法更新状态
        x_temp = x(t) + dx * dt;  % 暂时更新后的 x 值
        y_temp = y(t) + dy * dt;  % 暂时更新后的 y 值
        z_temp = z(t) + dz * dt;  % 暂时更新后的 z 值
        w_temp = w(t) + dw * dt;

        % 处理延迟
        if t <= delay_steps
            % 在延迟步数内，用初始值补全
            x_tau = x0;
            y_tau = y0;
            z_tau = z0;
            w_tau = w0;
        else
            % 获取延迟 τ 的值
            x_tau = x(t - delay_steps);
            y_tau = y(t - delay_steps);
            z_tau = z(t - delay_steps);
            w_tau = w(t - delay_steps);
        end

        % 计算修改函数值
        fx = (x_temp + x_tau) / 2 + cos(t * dt);  % x 的修正函数
        fy = (y_temp + y_tau) / 2 + cos(t * dt);  % y 的修正函数
        fz = (z_temp + z_tau) / 2 + cos(t * dt);  % z 的修正函数
        fw = (w_temp + w_tau) / 2 + cos(t * dt);

        % 模运算后更新状态
        x(t + 1) = mod(fx, 1);  % 用 fx 的模值更新 x(t+1)
        y(t + 1) = mod(fy, 1);  % 用 fy 的模值更新 y(t+1)
        z(t + 1) = mod(fz, 1);  % 用 fz 的模值更新 z(t+1)
        w(t + 1) = mod(fw, 1);  % 用 fw 的模值更新 w(t+1)
    end

    % 舍弃前 1000 个点，只保留之后的 10000 个
    valid_x = x(4001:end);
    valid_y = y(4001:end);
    valid_z = z(4001:end);
    valid_w = w(4001:end);

    X = valid_x(1:size_of_image)';
    Y = valid_y(1:size_of_image)';
    Z = valid_z(1:size_of_image)';
    W = valid_w(1:size_of_image)';

    % 评估指标
    % % 自相关系数
    % autocorr(valid_w);
    % 
    % 
    % % 绘制结果
    % %绘制X-Y坐标图
    % figure('NumberTitle', 'off', 'Name', '水平方向');%figure1是图像水平方向相关系数，纵坐标为：x+1,y
    % plot(round(valid_x*255),round(valid_y*255),'b.','linewidth',3,'markersize',3); 
    % axis([0 300 0 300]);
    % xlabel('x');  % x 轴标签
    % ylabel('y');  % y 轴标签

    % %绘制时序图
    % figure;  % 创建新图
    % subplot(4,1,1);
    % plot((1:length(valid_x)) * dt, valid_x);  % 绘制 x 的时间序列
    % title('x(t)');  % 图标题
    % xlabel('Time');  % x 轴标签
    % ylabel('x(t)');  % y 轴标签
    % 
    % subplot(4,1,2);
    % plot((1:length(valid_y)) * dt, valid_y);  % 绘制 y 的时间序列
    % title('y(t)');
    % xlabel('Time');
    % ylabel('y(t)');
    % 
    % subplot(4,1,3);
    % plot((1:length(valid_z)) * dt, valid_z);  % 绘制 z 的时间序列
    % title('z(t)');
    % xlabel('Time');
    % ylabel('z(t)');
    % 
    % subplot(4,1,4);
    % plot((1:length(valid_w)) * dt, valid_w);  % 绘制 z 的时间序列
    % title('w(t)');
    % xlabel('Time');
    % ylabel('w(t)');

    % % 绘制直方图来测试混沌序列的均匀性
    % figure;
    % subplot(4,1,1);
    % histogram(X, 50);
    % title('Histogram of x(t)');
    % xlabel('Value');
    % ylabel('Frequency');
    % 
    % subplot(4,1,2);
    % histogram(Y, 50);
    % title('Histogram of y(t)');
    % xlabel('Value');
    % ylabel('Frequency');
    % 
    % subplot(4,1,3);
    % histogram(Z, 50);
    % title('Histogram of z(t)');
    % xlabel('Value');
    % ylabel('Frequency');
    % 
    % subplot(4,1,4);
    % histogram(W, 50);
    % title('Histogram of w(t)');
    % xlabel('Value');
    % ylabel('Frequency');


    % % 绘图  
    % figure;  
    % 
    % subplot(2, 2, 1);  
    % plot(T, X);  
    % title('x vs time');  
    % xlabel('time');  
    % ylabel('x');  
    % 
    % subplot(2, 2, 2);  
    % plot(T, Y);  
    % title('y vs time');  
    % xlabel('time');  
    % ylabel('y');  
    % 
    % subplot(2, 2, 3);  
    % plot(T, Z);  
    % title('z vs time');  
    % xlabel('time');  
    % ylabel('z');  
    % 
    % subplot(2, 2, 4);  
    % plot(T, W);  
    % title('w vs time');  
    % xlabel('time');  
    % ylabel('w');  
    % 
    % % 三维相图  
    % figure;  
    % plot3(X, Y, Z);  
    % title('3D Phase Space');  
    % xlabel('x');  
    % ylabel('y');  
    % zlabel('z');  
    % grid on; 
end