%% P1 数据读入
clc
clear
data1 = csvread("D:\date_code\data\C_Hour.csv",1,1); 
data2 = data1(1:end,1) ;
% 采样间隔为一小时，这里没有选用极大值与极小值是因为其采样间隔不一

%% P2 归一化&平滑
% 使用多阶平滑，根据验证三阶平滑比较合适
[data3,re_data2] = mapminmax(data2', 0,1);
data_1 = smoothdata(data3','movmean',3);
data_2 = smoothdata(data_1,'movmean',3);
data_3 = smoothdata(data_2,'movmean',3);
data_4 = smoothdata(data_3,'movmean',3);
data_5 = smoothdata(data_4,'movmean',3);

%% P3 数据线性粗分段
% 找到最大值和最小值点（局部）
p1 = diff(data_3);
extra_index1(1,:) = [1, 1];
t = 2;
for i = 2:length(data_3)-1
	if p1(i-1) >= 0 && p1(i) < 0
		extra_index1(t, :) = [i,1];%极大值
		t = t+1;
	elseif p1(i-1) < 0 && p1(i) >= 0
		extra_index1(t, :) = [i,2]; %极小值
		t = t+1;
	end
end

%% 计算转折点
% 按极值分段(在平滑后的数据上)
t = 1;
for i = 2:size(extra_index1, 1)
	data_part1{t, 1} = data_3(extra_index1(i-1,1):extra_index1(i,1));
	data_part1{t, 2} = [extra_index1(i-1,1),extra_index1(i,1)];
	t = t+1;
end
% 计算重要点（转折）,转折点为距离直线的最大偏离距离点
t = 1;
for i = 1:size(data_part1, 1)
	if data_part1{i, 2}(2) - data_part1{i, 2}(1) > 5
		p1 = data_part1{i, 1};
		for j = 2:length(p1)-1
			p2(j-1) = abs(p1(1)+(length(p1)-1)*(p1(j)-p1(1))/(p1(end)-p1(j))-j);
		end
		[p3, p4] = max(p2);
		if p3 > 0.01 % 需要设定一个最大偏离距离，设定一个阀值
			turnpoint_index(t, :) = [data_part1{i,2}(1)+p4, 3];
			t = t+1;
		end
	end
end
% 与extra_index1合并，形成初始分段索引
part_index1 = sortrows([extra_index1; turnpoint_index],1);

%% 分段索引处理
%这里需要对斜率变化率进行计算，需要设置变化率阀值
t = 1; i = 1; r =1;
p1 = part_index1(:,1);
% 这里可能会有标记位置奇数与偶数的问题，可以考虑将代码执行多次
while i <= length(p1)-2
% 	p2 = data3(p1(i):p1(i+1));
% 	p3 = data3(p1(i):p1(i+2));
% 	p4 = 1:length(p2);
% 	p5 = 1:length(p3);
% 	k1 = regress(p2', [p4'*0.1, ones(length(p2),1)]);
% 	k2 = regress(p3', [p5'*0.1, ones(length(p3),1)]);
	k1 = (data3(p1(i+1)) - data3(p1(i)))/(p1(i+1)-p1(i));
	k1 = (data3(p1(i+2)) - data3(p1(i)))/(p1(i+2)-p1(i));
	% 计算斜率变化率
	change_ratio = abs((k2(1) - k1(1))/k1(1));
    f_change(r) = change_ratio;
    r =r+1;
	if change_ratio < 50 % 本处需要设定阀值，阀值为斜率变化率
		part_index(t, :) = part_index1(i,:);
		t = t+1; 
        i = i+1;
	else
		part_index(t,:) = part_index1(i,:);
		t = t+1; 
        i = i+2;
	end
end

%% tips:
% Dear programmer:
% When I wrote this code, only god and I kenw how it work. Now only god knows it!
% Good lucky!
	
