%% 数据读入
clc
clear
% 读入数据
date_ah = csvread("D:\date_code\date\B_Hour.csv",1,1);
date_ahm =mean(date_ah(1:74000, 2:3), 2); % 数据中极大值与极小值之间的平均
date_ahtest = mean(date_ah(74001:end,2:3),2);

%% 数据提取
% 训练数据提取
p_train1 = pick(date_ahm,date_ahm); % 测试数据提取
p_test1 = pick(date_ahtest,date_ahtest); % 这里提取改为只返回数据，归一化和反归一化化部分，不放在此处

%% 构建训练集
% 归一化
% step1 训练集归一化
[THR,SORH,KEEPAPP,CRIT]=ddencmp('den','wp',date_ahm');
xx=wdencmp('gbl',date_ahm','db6',3,THR,SORH,KEEPAPP);
date_s=smoothdata(xx,'movmean',4); % 使用了两次平滑减小临近位置的极值
date_f = diff(date_s);
[date_diff, re] = mapminmax(date_f,0,1);
t1 =randperm(size(p_train1,1));
train_date = p_train1(t1,:); % 训练集将会被打乱
p_train11 = train_date(:,5)';
for i = 1:size(p_train11,2)
    t2 = mapminmax('apply',p_train11{1,i}',re);
    p_train{1,i} =t2';
end
% step2 测试集归一化
for i = 1:size(p_test1, 1)
    t2 = mapminmax('apply',p_test1{i,5}',re);
    p_test{1,i} =t2';
end
% 重写训练输出
for i = 1:size(train_date,1)
    switch train_date{i, 6}
        case 1
            t_train{1,i} = [1;0];
        case 2
            t_train{1,i} = [0;1];
    end
end


%% 建立神经网络
net = newrb(p_train,t_train, 0.2, 0.75);

%% 预测
t_pr =sim(net,p_test);

%% 效果评估
%将结果转换为标记
for i =1:size(t_pr,2)
    p_out(i,1) = find(t_pr{1,i} ==max(t_pr{1,i}));
end
 % 正确率
 t_right = 0;
 for i =1:size(t_pr,2)
    if p_out(i,1)-p_test1{i,6}==0
        t_right = t_right+1;
    end
end
right_ridio = t_right/size(t_pr,2);

