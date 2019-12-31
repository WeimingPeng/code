%% P1 数据读入
clc
clear
data_bh = csvread("D:\date_code\data\C_Hour.csv",1,1);
data_in = data_bh(:, 1); % 数据中极大值与极小值之间的平均

%% P2 平滑
[THR,SORH,KEEPAPP,CRIT]=ddencmp('den','wp',data_in);
data = wdencmp('gbl',data_in,'db6',3,THR,SORH,KEEPAPP);
% data=smoothdata(p1,'movmean',4); %使用两次平滑减小临近位置的极值
% data_s = smoothdata(data,'movmean',4);

%% P3 建立标记集data_mark
p1 = diff(data);
data_mark(1) = 2;
for i = 2:length(data)-1
	if p1(i-1)>0 && p1(i)<0
		data_mark(i) = 1; %极大值
	elseif p1(i-1)<0 && p1(i)>0
		data_mark(i) = 2; %极小值
	elseif p1(i-1)<0 && p1(i)<0
		data_mark(i) = 3; %下降区间
	elseif p1(i-1)>0 && p1(i)>0
		data_mark(i) = 4; %上升区间
	end
end
data_mark = [data_mark,2]; 
% 极值点的索引
extra_index = find(data_mark < 3);

%% P4 data_mark处理，对区间小于5的部分处理。
p2 = diff(extra_index);
t = 2;
r=0;g=0;f=0;
data_mark1 = data_mark;
while t<length(extra_index)-1
    if p2(t)<5
        if data_mark(extra_index(t)-1) == data_mark(extra_index(t+1)+1)
            p3 = data_mark(extra_index(t)-1);
        r=r+1;
        elseif data_mark(extra_index(t)-1) == 3 || data_mark(extra_index(t)-1) == 4
            p3 = data_mark(extra_index(t)-1);
            g=g+1;
        elseif data_mark(extra_index(t+1)+1) == 3 || data_mark(extra_index(t+1)+1) == 4
            p3 = data_mark(extra_index(t+1)+1);
            f=f+1;
        end
        data_mark1(extra_index(t):extra_index(t+1)) = p3*ones(1,p2 (1,t)+1);
        t = t+2;
    else
        t=t+1;
    end
end
extra_index1 = find(data_mark1 < 3);

%% P5 查找极值点等值数据
% 对数据进行适度取位后，按指定阀值确定等值数据
data1 = floor(data*10000)/10000; %取小数后4位
p=1;
for j = 2:length(extra_index1)
    p1 = extra_index1(j);        
    p2= find(abs(data1(p1:end)-data1(p1))<0.0001);
    datamark1{p, 1} = p2+p1-1;
    p=p+1;
end

%% P6 正例提取，计算出每个极值区间的特征值
for i = 1:length(extra_index1)-1
    p1 = data_bh(extra_index1(i):extra_index1(i+1));
    postive_feat(i,1:4) = datafeature(p1');
    postive_feat(i,5) = extra_index1(i+1);
end

%% P7负例提取
% 使用的是三段的提取方法
p7_out = [];
extra_mark(1,:) = diff(extra_index);
for i =2:length(extra_index)-1
	if extra_mark(i)<10
		p1 = extra_index(i)+floor(rand(1)*(extra_mark(i)-2)+2);
		p2 = data_bh(extra_index(i):p1);
		p3(1,1:4) = datafeature(p2');
		p3(5) = p1;
		p7_out = [p7_out;p3];
	elseif extra_mark(i)>=10&&extra_mark(i)<20
		p1 = extra_index(i)+floor(rand(1)*(extra_mark(i)/2-2)+2);
		p2 = data_bh(extra_index(i):p1);
		p3(1,1:4) = datafeature(p2');
		p3(5) = p1;
		p7_out = [p7_out;p3];

		p4 = extra_index(i)+floor((rand(1)+1)*(extra_mark(i)/2));
		p5 = data_bh(extra_index(i):p4);
		p6(1,1:4) = datafeature(p5');
		p6(5) = p4;
		p7_out = [p7_out;p6];
	elseif extra_mark(i) >= 20
		p1 = extra_index(i)+floor(rand(1)*(extra_mark(i)/3-2)+2);
		p2 = data_bh(extra_index(i):p1);
		p3(1,1:4) = datafeature(p2');
		p3(5) = p1;
		p7_out = [p7_out;p3];

		p4 = extra_index(i)+floor((rand(1)+1)*(extra_mark(i)/3));
		p5 = data_bh(extra_index(i):p4);
		p6(1,1:4) = datafeature(p5');
		p6(5) = p4;
		p7_out = [p7_out;p6];

		p7 = extra_index(i)+floor((rand(1)+2)*(extra_mark(i)/3));
		p8 = data_bh(extra_index(i):p7);
		p9(1,1:4) = datafeature(p8');
		p9(5) = p7;
		p7_out = [p7_out;p9];
	end
end

%% P8宋的验证数据提取
p_out = [];
for i = 1:length(extra_index1)-1
	p1 = datamark1{i,1}(2:end);
	p_out = [p_out;p1];
end
t_index = unique(p_out,'sorted');


%% 特征提取函数
% input必须为列向量
function a = datafeature(b)
    % 斜率k
    p1 = 1:length(b);
    p2 = regress(b, [ p1'*01,ones(length(b),1)]);
    a(1,1) = p2(1);
    %差分方差
    a(1,2) = var(diff(b));
    %差分平均
    a(1,3) = mean(diff(b));
    %长度
    a(1,4) = length(b);
end

