%% 读入数据
clc
clear
% 数据特征提取主程序对小时数据进行提取

date_1 = csvread('D:\date_code\date\B_Hour.csv',1,1);
date_2 = date_1(:, 2:3)'; % 提取最大值和最小值两列
date_3 = date_2';
date_0 = (mean(date_2))';%取小时内最大值最小值平均
 r = floor(date_0*100)/100;

%% 平滑
[THR,SORH,KEEPAPP,CRIT]=ddencmp('den','wp',date_bhm);
xx=wdencmp('gbl',date_bhm,'db6',3,THR,SORH,KEEPAPP);
date_s=smoothdata(xx,'movmean',4); %使用两次平滑减小临近位置的极值

%% 极值查找与提取
datemin_index(:, 1) = find(diff(sign(diff(date_s)))>0)+1;%min
datemin_index(:, 2) = 1;
datemax_index(:, 1) = find(diff(sign(diff(date_s)))<0)+1;
datemax_index(:, 2) = 2;
extremum_index = sortrows([datemin_index; datemax_index]);
 c = diff(extremum_index(:, 1)); % 在这里查看极大值间隔的分布特性

%%
dd=[];
dd((1:length(b)))=b';
dd((length(b)+1):(length(b)+length(c)))=c';
d=dd';
[e,f]=sort(d);%e与f合并在一起，再增加一列判断e的值为极大值还是极小值
test_p(1,1)=1;test_p(1,2)=2;test_p(1,3)=x(1);
test_p(2:length(e)+1,1)=e;
for i=1:length(e)
    if f(i)<=(length(b))
        test_p(i+1,2)=1;%min
    else
        test_p(i+1,2)=2;
    end
end
test_p(2:length(e)+1,3)=x(test_p(2:length(e)+1,1));
y=diff(test_p(:,3));%这里用的是一个差分来去到想要的极值点，存在一定的问题，详细问题见笔记论文paper1
ya=abs(y);
s=2;
yb=[];yb(1,:)=test_p(1,:);
for i=1:length(ya)
    if ya(i)>0.008 % 连续极值间隔应该大于这个值  
        yb(s,:)=test_p(i+1,:);
        s=s+1;
    end
end
yc=[];
yd=1;
for i=1:length(yb(:,1))-1
    if yb(i,2)==2
        if yb(i+1,2)==1
            yc(yd,:)=yb(i,:);
             yd=yd+1;
        end
    elseif yb(i,2)==1
       if yb(i+1,2)==2
            yc(yd,:)=yb(i,:);
             yd=yd+1;
        end
    end
end %yc 为最终提取的结果
%% 极值点数据间数据提取
fss=1;
fsr=1;
sd=[];
t=1;
q=1;
ws=[];
test_ss=[];
fsd=length(test_p(:,1));
for k=1:length(test_p(:,1))-1
     if test_p(k,2)==1
         ws(q,1)=test_p(k,1);
         ws(q,2)=1;
         q=q+1;
         ye=sign(x(test_p(k,1):end)-x(test_p(k,1)));
        qq=find(diff(ye)<0)+test_p(k,1);
        for g=1:length(qq)
            for i=1:(fsd-1)
                j=test_p(i,1);
                u=test_p(i+1,1);
                if  (j<qq(g))&&(qq(g)<u)
                       if ((qq(g)-j)>12)&&(qq(g)+6<93338)&&((u-qq(g))>5)
                           %此处控制所取中间点需要满足取点要求
                           if test_p(i,2)==2%max到min中间
                                test_ss(t,1)=qq(g);
                                ws(q,1)=qq(g);
                                ws(q,2)=3;
                                test_ss(t,2)=3;
                                test_ss(t,3:15)=a(qq(g)-6:qq(g)+6);
                                t=t+1;
                                q=q+1;
                           end
                       end
                end
            end
        end
     elseif test_p(k,2)==2
         ws(q,1)=test_p(k,1);
         ws(q,2)=2;
         q=q+1;
         ye=sign(x(test_p(k,1):end)-x(test_p(k,1)));
        qq=find(diff(ye)>0)+test_p(k,1);
        for g=1:length(qq)
            for i=1:(fsd-1)
                j=test_p(i,1);
                u=test_p(i+1,1);
                if  (j<qq(g))&&(qq(g)<u)
                    fss=fss+1;
                       if ((qq(g)-j)>12)&&(qq(g)+6<93338)&&((u-qq(g))>5)
                           if test_p(i,2)==1%min到max中间
                                test_ss(t,1)=qq(g);
                                ws(q,1)=qq(g);
                                ws(q,2)=4;
                                test_ss(t,2)=4;
                                test_ss(t,3:15)=a(qq(g)-6:qq(g)+6);
                                t=t+1;
                                q=q+1;
                           end
                      end
                end
            end
        end
     end
end
test_w=[];
test_w(1:length(yc),:)=yc(:,1:2);
test_w((length(yc)+1):(length(yc)+length(test_ss(:,1))),:)=test_ss(:,1:2);
test_ws=sortrows(test_w,1);%按指定的一列排序
test_wss=unique(test_ws,'rows') ;%去掉重复项
xa=[];xb=[];xc=1;xd=1;
for i=1:length(test_wss(:,1))
    if test_wss(i,2)==3
        xa(xc)=test_wss(i,1);
        xc=xc+1;
    elseif test_wss(i,2)==4
        xb(xd)=test_wss(i,1);
        xd=xd+1;
    end
end %test_wss为输出数据，数据存取的时候重命名为date_index


    
