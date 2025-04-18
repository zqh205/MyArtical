function [x, y, z] = modified_lorenz_chaos(x0, y0, z0, size_of_image)
    % 参数设置
    a = 20;  % 洛伦兹系统中 σ 参数,大于0即可，影响较小
    b = 50;  % 洛伦兹系统中 ρ 参数,一般大于24.7，表现混沌
    c = 8;   % 洛伦兹系统中 β 参数，经典值8/3，其他值比如8也可以

    % % 初始条件
    % x0 = 0.5;  % x 的初始值
    % y0 = 0.3;  % y 的初始值
    % z0 = 0.3;  % z 的初始值
    % size_of_image = 256 * 256;

    % 时间设置
    total_iterations = 104000;  % 总的迭代步数，前 1000 要舍弃
    T = total_iterations * 1;  % 根据迭代和步长设定总时间
    dt = 1;   % 时间步长

    % 延迟时间步数
    tau = 1;  % 延迟时间，以秒计
    delay_steps = floor(tau / dt);  % 将延迟时间转换为时间步数

    % 初始化，创建存储 x, y, z 的数组
    x = zeros(1, total_iterations);  % 用于存储 x 值的数组
    y = zeros(1, total_iterations);  % 用于存储 y 值的数组
    z = zeros(1, total_iterations);  % 用于存储 z 值的数组

    % 设置初始值
    x(1) = x0;  % x 的初始条件
    y(1) = y0;  % y 的初始条件
    z(1) = z0;  % z 的初始条件

    % 数值积分主循环
    for t = 1:total_iterations-1
        % 洛伦兹系统的状态方程计算
        dx = a * (y(t) - x(t));  % x 的微分方程
        dy = x(t) * (b - z(t)) - y(t);  % y 的微分方程
        dz = x(t) * y(t) - c * z(t);  % z 的微分方程

        % 使用欧拉法更新状态
        x_temp = x(t) + dx * dt;  % 暂时更新后的 x 值
        y_temp = y(t) + dy * dt;  % 暂时更新后的 y 值
        z_temp = z(t) + dz * dt;  % 暂时更新后的 z 值

        % 处理延迟
        if t <= delay_steps
            % 在延迟步数内，用初始值补全
            x_tau = x0;
            y_tau = y0;
            z_tau = z0;
        else
            % 获取延迟 τ 的值
            x_tau = x(t - delay_steps);
            y_tau = y(t - delay_steps);
            z_tau = z(t - delay_steps);
        end

        % 计算修改函数值
        fx = (x_temp + x_tau) / 2 + sin(t * dt);  % x 的修正函数
        fy = (y_temp + y_tau) / 2 + sin(t * dt);  % y 的修正函数
        fz = (z_temp + z_tau) / 2 + sin(t * dt);  % z 的修正函数

        % 模运算后更新状态
        x(t + 1) = mod(fx, 1);  % 用 fx 的模值更新 x(t+1)
        y(t + 1) = mod(fy, 1);  % 用 fy 的模值更新 y(t+1)
        z(t + 1) = mod(fz, 1);  % 用 fz 的模值更新 z(t+1)
    end

    % 舍弃前 1000 个点，只保留之后的 10000 个
    valid_x = x(4001:end);
    valid_y = y(4001:end);
    valid_z = z(4001:end);

    x = valid_x(1:size_of_image)';
    y = valid_y(1:size_of_image)';
    z = valid_z(1:size_of_image)';
    % 自相关系数
    % autocorr(valid_x);

    % 作为图像输出
    % z = reshape(z, 256,256);
    % figure;
    % imshow(z);
    % title('z');
    
    % 均值、方差
    % meanValue = mean(x);
    % varianceValue = var(x);
    % 
    % disp(['Mean: ', num2str(meanValue)]);
    % disp(['Variance: ', num2str(varianceValue)]);



    %% 绘制结果
    % %绘制X-Y坐标图
    % figure('NumberTitle', 'off', 'Name', '水平方向');%figure1是图像水平方向相关系数，纵坐标为：x+1,y
    % plot(round(valid_x*255),round(valid_y*255),'b.','linewidth',3,'markersize',3); 
    % axis([0 300 0 300]);

    % figure;  % 创建新图
    % subplot(3,1,1);
    % plot((1:length(valid_x)) * dt, valid_x);  % 绘制 x 的时间序列
    % title('Modified x(t)');  % 图标题
    % xlabel('Time');  % x 轴标签
    % ylabel('x(t)');  % y 轴标签
    % 
    % subplot(3,1,2);
    % plot((1:length(valid_y)) * dt, valid_y);  % 绘制 y 的时间序列
    % title('Modified y(t)');
    % xlabel('Time');
    % ylabel('y(t)');
    % 
    % subplot(3,1,3);
    % plot((1:length(valid_z)) * dt, valid_z);  % 绘制 z 的时间序列
    % title('Modified z(t)');
    % xlabel('Time');
    % ylabel('z(t)');
    % 
    % % 绘制直方图来测试混沌序列的均匀性
    % figure;
    % subplot(3,1,1);
    % histogram(x, 50);
    % title('Histogram of x(t)');
    % xlabel('Value');
    % ylabel('Frequency');
    % 
    % subplot(3,1,2);
    % histogram(y, 50);
    % title('Histogram of y(t)');
    % xlabel('Value');
    % ylabel('Frequency');
    % 
    % subplot(3,1,3);
    % histogram(z, 50);
    % title('Histogram of z(t)');
    % xlabel('Value');
    % ylabel('Frequency');
end