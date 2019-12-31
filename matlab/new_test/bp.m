%% #1 数据读入
clc
clear
data_bhh = csvread("D:\date_code\data\C_Hour.csv",1,1); 

%归一化
[data_bh1, re] = mapminmax(data_bhh',0,1);
data_bh = data_bh1';

data1 = data_bh(1:80016, 1); % 训练数据
data2 = data_bh(80017:end, 1); % 测试数据
% 采样间隔为一小时，这里没有选用极大值与极小值是因为其采样间隔不一

%% #2 平滑
% 小波平滑
[THR,SORH,KEEPAPP,CRIT]=ddencmp('den','wp',data1);
data_1 = wdencmp('gbl',data1,'db6',3,THR,SORH,KEEPAPP);

% 移动平均3-5位（使用三次）
data_2 = smoothdata(data1,'movmean',3);
data_3 = smoothdata(data_2,'movmean',3);
data_4 = smoothdata(data_3,'movmean',3);
data_mv = data_4;
data_wv = data_1;

%% # 数据分段
% 训练集
for i = 1:length(data_mv)-24
    p_train(1:24,i) = data_mv(i:i+23);
    t_train(1,i) = data_mv(i+24);
end
% 测试集
for i = 1:length(data2)-35
    p_test(1:24,i) = smoothdata(data2(i:23+i),'movmean',3);
    t_test(1,i) = data2(24+i);
    test_dat(1:12,i) = data2(i+24:i+35);
end

%% 建立bp网络
net = newff(p_train,t_train,[30,25],{'logsig','logsig'},'traingdm');
net.trainParam.epochs = 10000;
net.trainParam.goal = 1e-5;
net.trainParam.lr = 0.01;
net.trainParam.mc = 0.9;

net= train(net, p_train,t_train);

%% 仿真
for i = 1:length(t_test)
    p2 = p_test(:,i);
    t = 1;
    for j =1:12 %这里计算趋势的长度，与分段对应
        p1 = sim(net,p2);
        p2 = [p2(2:end); p1];
        p_out(t, 1) = p1;
        t = t+1;
    end
    pr_out(1:12,i) = p_out;
end

%% tips:
% Dear programmer:
% When I wrote this code, only god and I kenw how it work. Now only god knows it!
% Good lucky!