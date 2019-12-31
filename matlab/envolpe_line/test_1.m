%本代码用于测试差分包络线
%% P1 数据读入
clc
clear
close
date_bh = csvread("D:\date_code\data\A_Hour.csv",1,1);
data_in =mean(date_bh(:, 2:3), 2); % 数据中极大值与极小值之间的平均

%% 数据提取
data_diff = diff(data_in);
diff_up= data_diff(data_diff >= 0);
diff_down = data_diff(data_diff < 0);

%% 测试1000个数据
data1 = data_in(1:1000);
data1_diff = diff(data1);
up_index = find(data1_diff >= 0);
down_index = find(data1_diff < 0);
% 采样包络线
p1 = 1:2:length(up_index);
p2 = 1:2:length(down_index);
p3 = up_index(p1);
p4 = down_index(p2);
%% 简单测试
p5 = p3(p3 < 105);
p6 = p4(p4 < 105);
p7 = [16,36,40,43,52,58,66,82];
subplot(2,1,1);
plot(p5,data1_diff(p5))
hold on
plot(p6,data1_diff(p6))
hold off
subplot(2,1,2)
plot(data_in(1:105))
hold on
plot(p7,data_in(p7))
%% 测试绘图
figure(1)
plot(data1)
figure(2)
subplot(2,1,1)
plot(p3,data1_diff(p3))
hold on
plot(p4,data1_diff(p4))
%plot(data1_diff,'y')
hold off
subplot(2,1,2)
plot(data1)
