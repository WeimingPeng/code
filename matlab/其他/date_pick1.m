%% 数据读入
clc
clear
% 读入数据统一为2004年4月5日0分0时0秒
date_ah = csvread("D:\date_code\date\A_Hour.csv",1,1);
date_ahm =mean(date_ah(1:74000, 2:3), 2); % 数据中极大值与极小值之间的平均
date_ahtest = mean(date_ah(74001:end,2:3),2);

%% 归一化
[date_a1, re] = mapminmax(date_ahm',0,1); % 线性缩放到0，1之间
date_a2 = mapminmax("apply",date_ahtest',re); %对测试数据按训练数据归一化 

%% 下面将对训练集进行提取
%% 平滑
[THR,SORH,KEEPAPP,CRIT]=ddencmp('den','wp',date_a1');
xx=wdencmp('gbl',date_a1','db6',3,THR,SORH,KEEPAPP);
date_s=smoothdata(xx,'movmean',4); % 使用了两次平滑减小临近位置的极值

%% 极值查找与提取
datemin_index(:, 1) = find(diff(sign(diff(date_s)))>0)+1;%min
datemin_index(:, 2) = 1;
datemax_index(:, 1) = find(diff(sign(diff(date_s)))<0)+1;
datemax_index(:, 2) = 2;
extremum_index = sortrows([datemin_index; datemax_index]);
c = diff(extremum_index(:, 1));  % 在这里查看极大值间隔的分布特性

%% 按某一数据长度提取数据。目前定位一天
% 下面将对临界收敛数据进行处理 ;即所提取数据为临街点后收敛数据
date1 = extremum_index(3:end,:);
apart = cell(length(date1(:,1)),2);
for j = 1:length(date1(:,1))
    apart{j,1} = date1(j,:);
    apart{j,2} = date_ahm(date1(j,1)-24:date1(j,1));
end
% 下面将提取发散数据
% step1 筛选出连续发散区间，滚动取值
% 提取连续变化说我数据
t1 = 1;
for i = 2:length(date1(:,1))
    if date1(i,1) - date1(i-1,1)>30
        date_1(t1, :) =[date1(i-1), date1(i,1), date1(i-1,2)];
        t1 = t1+1;
    end
end
% 提取目标数据
t2 =1;
for i = 1: length(date_1(:,1))
    for j = 1:(date_1(i,2)-date_1(i,1)-29)
        bpart{t2,1} = date_1(i,3); % 标记代号1是上升区间，标记代号2为下降区间
        bpart{t2,2} = date_1(i,1)+24+j;
        bpart{t2,3} = date_ahm(date_1(i,1)+j : date_1(i,1)+j+24);
        t2 =t2+1;
    end
end
