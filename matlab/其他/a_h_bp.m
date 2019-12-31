%% 数据读入与数据集划分
load('date_1.mat')
date_c = date_1(:,3:14)';
[ p_train ,ps] = mapminmax(date_c(1:11,1:5000));
[t_train,ts] = mapminmax(date_c(12,1:5000));
p_test = mapminmax(date_c(1:11,5001:end));
%% 建立bp神经网络
net = newff(p_train,t_train,[25,15],{'tansig','tansig'},'trainlm');
net.trainParam.epochs = 1000;
net.trainParam.lr = 1e-8;
net.trainParam.goal = 1e-4;
net = train(net,p_train,t_train);