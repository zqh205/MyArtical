%图像相邻像素相关性
 
%随机选取图像A中相邻的N个像素点，分别计算在水平，垂直，正对角和反对角方向上的相关系数
 
%注意是算密文图像 R G B分量分开算
 
function r = im_correlation(A,N)
    %像素值要变为0~255
 
    A=double(A);[m,n]=size(A);r=zeros(1,4);
     
    x1=mod(floor(rand(1,N)*10^10),m-1)+1;
     
    x2=mod(floor(rand(1,N)*10^10),m)+1;
     
    x3=mod(floor(rand(1,N)*10^10),m-1)+2;
     
    y1=mod(floor(rand(1,N)*10^10),n-1)+1;
     
    y2=mod(floor(rand(1,N)*10^10),n)+1;
     
    u1=zeros(1,N);u2=zeros(1,N);u3=zeros(1,N);u4=zeros(1,N);
     
    v1=zeros(1,N);v2=zeros(1,N);v3=zeros(1,N);v4=zeros(1,N);
     
    for i=1:N
     
        u1(i)=A(x1(i),y2(i));v1(i)=A(x1(i)+1,y2(i));
     
        u2(i)=A(x2(i),y1(i));v2(i)=A(x2(i),y1(i)+1);
     
        u3(i)=A(x1(i),y1(i));v3(i)=A(x1(i)+1,y1(i)+1);
     
        u4(i)=A(x3(i),y1(i));v4(i)=A(x3(i)-1,y1(i)+1);
     
    end
     
    r(1)=mean((u1-mean(u1)).*(v1-mean(v1)))/(std(u1,1)*std(v1,1));
     
    r(2)=mean((u2-mean(u2)).*(v2-mean(v2)))/(std(u2,1)*std(v2,1));
     
    r(3)=mean((u3-mean(u3)).*(v3-mean(v3)))/(std(u3,1)*std(v3,1));
     
    r(4)=mean((u4-mean(u4)).*(v4-mean(v4)))/(std(u4,1)*std(v4,1));
     
     
     
    figure('NumberTitle', 'off', 'Name', '水平方向');%figure1是图像水平方向相关系数，纵坐标为：x+1,y
     
    plot(u1,v1,'b.','linewidth',3,'markersize',3);
     
    axis([0 300 0 300]);
     
    set(gca,'XTick',0:50:300,'YTick',0:50:300,'fontsize',22,'fontname','times new roman');
     
    set(gca,'XTickLabel',{'0','50','100','150','200','250','300'});
     
    set(gca,'YTickLabel',{'0','50','100','150','200','250','300'});
     
    % xlabel('Pixel gray value on location(\itx\rm,\ity\rm)');
    % 
    % ylabel('Pixel gray value on location(\itx\rm+1,\ity\rm)');
     
    figure('NumberTitle', 'off', 'Name', '垂直方向');%figure2是图像垂直方向相关系数，纵坐标为：x,y+1
     
    plot(u2,v2,'b.','linewidth',3,'markersize',3);
     
    axis([0 300 0 300]);
     
    set(gca,'XTick',0:50:300,'YTick',0:50:300,'fontsize',22,'fontname','times new roman');
     
    set(gca,'XTickLabel',{'0','50','100','150','200','250','300'});
     
    set(gca,'YTickLabel',{'0','50','100','150','200','250','300'});
     
    % xlabel('Pixel gray value on location(\itx\rm,\ity\rm)');
    % 
    % ylabel('Pixel gray value on location(\itx\rm,\ity\rm+1)');
     
    figure('NumberTitle', 'off', 'Name', '正对角方向');%figure3是图像正对角线方向相关系数，纵坐标为：x+1,y+1
     
    plot(u3,v3,'b.','linewidth',3,'markersize',3);
     
    axis([0 300 0 300]);
     
    set(gca,'XTick',0:50:300,'YTick',0:50:300,'fontsize',22,'fontname','times new roman');
     
    set(gca,'XTickLabel',{'0','50','100','150','200','250','300'});
     
    set(gca,'YTickLabel',{'0','50','100','150','200','250','300'});
     
    % xlabel('Pixel gray value on location(\itx\rm,\ity\rm)');
    % 
    % ylabel('Pixel gray value on location(\itx\rm+1,\ity\rm+1)');
     
    figure('NumberTitle', 'off', 'Name', '反对角方向');%figure4是图像反对角方向相关系数，纵坐标为：x-1,y+1

    plot(u4,v4,'b.','linewidth',3,'markersize',3);

    axis([0 300 0 300]);

    set(gca,'XTick',0:50:300,'YTick',0:50:300,'fontsize',22,'fontname','times new roman');

    set(gca,'XTickLabel',{'0','50','100','150','200','250'});

    set(gca,'YTickLabel',{'0','50','100','150','200','250'});
     
    % xlabel('Pixel gray value on location(\itx\rm,\ity\rm)');
    % 
    % ylabel('Pixel gray value on location(\itx\rm-1,\ity\rm+1)');
     
end