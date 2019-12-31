%% 编程测试
clc
clear
close
data_bh = csvread("D:\date_code\data\A_Hour.csv",1,1);
data_in =mean(data_bh(:, 2:3), 2); % 数据中极大值与极小值之间的平均
%% 测试

data1 = data_bh(1:1000,1);
dat_diff = diff(data1);
diff_d = diff(dat_diff);
[a,b,c,d] = jbtest(diff_d);
disp(num2str([a,b,c,d]))