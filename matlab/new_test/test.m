%% 测试datapick当数据为1：1000时

plot(data)
p1 = smoothdata(data,'movmean',3);
hold on
plot(p1(1:1000))
hold off
%% 测试extra_index1间隔
a = diff(extra_index1(:,1));
b = find(diff(extra_index1(1:end-1,2)) == 0);
c = diff(extra_index1(:,2));
disp(num2str(b));
% 对异常数据适度画图
plot(extra_index1(8472):extra_index1(8475),data1(extra_index1(8472):extra_index1(8475)))
e = find(a ==1);

%% 对得到极值数处理
p1 = data2(extra_index1(:,1));
figure(1)
plot(p1)
