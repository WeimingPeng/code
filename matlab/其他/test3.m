%% P1 ���ݶ���
clc
clear
date_bh = csvread("D:\date_code\data\A_Hour.csv",1,1);
date_p1 =mean(date_bh(:, 2:3), 2); % �����м���ֵ�뼫Сֵ֮���ƽ��
date_bhm =date_p1(1: 75000);
 % r = floor(date_p1*100)/100; % ��������һЩС��λ

%% P2 ƽ��
[THR,SORH,KEEPAPP,CRIT]=ddencmp('den','wp',date_bhm);
p1=wdencmp('gbl',date_bhm,'db6',3,THR,SORH,KEEPAPP);
date_s=smoothdata(p1,'movmean',4); %ʹ������ƽ����С�ٽ�λ�õļ�ֵ

%% ��ֵ��������ȡ
datemin_index(:, 1) = find(diff(sign(diff(date_s)))>0)+1;%min
datemin_index(:, 2) = 1;
datemax_index(:, 1) = find(diff(sign(diff(date_s)))<0)+1;
datemax_index(:, 2) = 2;
a_index = sortrows([datemin_index; datemax_index]);
b_index =a_index;

%% �޳�һЩ������
% step1 ���2����3����4������
for i = 1:size(a_index, 1)-4
    if a_index(i+2,1) - a_index(i+1,1)<5 && a_index(i+1,1) - a_index(i,1)<5 %�������������䳤��С��5
        b_index(i+1, :) = [0, 0];
        b_index(i+2, 2) = 3;
    elseif a_index(i+3,1) - a_index(i+2,1)<6 && ...
            a_index(i+2,1) - a_index(i+1,1)<6&& a_index(i+1,1) - a_index(i,1)<6 %�������������䳤��С��6
        b_index(i+1,:) = [0, 0];
        b_index(i+2,:) = [0, 0];
        b_index(i+3,2) = 4;
    elseif a_index(i+4,1) - a_index(i+3,1)<6 && ...
            a_index(i+3,1) - a_index(i+2,1)<6&& a_index(i+2,1) - a_index(i+1,1)<6&&...
            a_index(i+1,1) - a_index(i,1)<6 % ò��û������������4�������䳤��С��6
    end  
end
b_index1 = sortrows(b_index);
t = 1;
for i = 1:size(b_index1,1)
    if b_index1(i,2) > 0
		b_index2(t,:) = b_index1(i,:);
		t = t+1;
    end
end
b_index =b_index2;
c_index =b_index;
%% ������γ����쳣
for i =2:size(b_index,1)
    if b_index(i,1)-b_index(i-1,1)<5
        c_index(i,:) = [0,0];
        c_index(i-1,:) = [0,0];
    end
end
c_index1 = sortrows(c_index);
t=1;
for i = 1:size(c_index1,1)
    if c_index1(i,2) > 0
        c_index2(t,:) = c_index1(i,:);
        t = t+1;
    end
end
extre_index = c_index2;

%% ԭʼ���ݷֶ�
% ��ԭʼ���ݰ���ֵ����ΪС�Σ���ʼ��������һ����ֵ��
date_part{1,1} = date_bhm(1:extre_index(1,1)); %ԭʼ����
a = smoothdata(date_part{1,1},'movmean',3); %ԭʼ�����ƶ�ƽ��
date_part{1,2} =diff(a); 
date_part{1,3} = length(a);
date_part{1,4} = extre_index(1,:);
for i = 2:size(extre_index,1)
    date_part{i,1} = date_bhm(extre_index(i-1,1)+1:extre_index(i,1));
    a = smoothdata(date_part{i,1},'movmean',3);
    date_part{i,2} =diff(a);
    date_part{i,3} = extre_index(i,1) - extre_index(i-1,1);
    date_part{i,4} = extre_index(i,:);
end

%% �ֶ�����������ȡ������������ 
for i = 2:size(date_part,1)
    date_feature{i-1,1} = length(date_part{i-1,1});
    date_feature{i-1,2} = length(date_part{i,1});
    e1 = date_part{i-1,1};
    ee = 1:length(e1);
    e2 = regress(e1,[ee'*0.1, ones(length(date_part{i-1,1}),1)]);
    date_feature{i-1,3} = e2(1);
    e1 = date_part{i,1};
    ee1 = 1:length(e1);
    e2 = regress(e1,[ee1'*0.1,ones(length(date_part{i,1}),1)]);
    date_feature{i-1,4} =e2(1);
    date_feature{i-1,5} = mean(date_part{i-1,2});
    date_feature{i-1,6} = mean(date_part{i,2});
    date_feature{i-1,7} = date_part{i,4}(1);
    date_feature{i-1,8} = date_part{i,4}(2);
end

%% ����������ȡ f
f=1;
for i = 2:size(date_part,1)
	if date_part{i,3} < 10 %���ڵ�������С��10����ȡһ������
		f1 = 2+floor(rand(1)*(date_part{i,3}-2));
		f2 = date_part{i,1}(1:f1);
		date_feature1{f,1} = length(date_part{i,1});
		date_feature1{f,2} = f1;
    	f3 = date_part{i-1,1};
    	f4 = 1:length(f3);
    	f5 = regress(f3,[f4'*0.1, ones(length(date_part{i-1,1}),1)]);
    	date_feature1{f,3} = f5(1);
    	f6 = 1:length(f2);
    	f7 = regress(f2,[f6'*0.1,ones(length(f2),1)]);
    	date_feature1{f,4} =f7(1);
    	date_feature1{f,5} = mean(date_part{i-1,2});
		date_feature1{f,6} = mean(date_part{i,2}(1:f1-1));
        date_feature1{f,7} = date_part{i,4}(1)+f1;
		date_feature1{f,8} = 5;
		f = f+1;

	elseif date_part{i,3} < 20&& date_part{i,3} >= 10
		f1 = 2+floor(rand(1)*(date_part{i,3}/2-2));
		f2 = date_part{i,1}(1:f1);
		date_feature1{f,1} = length(date_part{i,1});
		date_feature1{f,2} = f1;
    	f3 = date_part{i-1,1};
    	f4 = 1:length(f3);
    	f5 = regress(f3,[f4'*0.1, ones(length(date_part{i-1,1}),1)]);
    	date_feature1{f,3} = f5(1);
    	f6 = 1:length(f2);
    	f7 = regress(f2,[f6'*0.1,ones(length(f2),1)]);
    	date_feature1{f,4} =f7(1);
    	date_feature1{f,5} = mean(date_part{i-1,2});
		date_feature1{f,6} = mean(date_part{i,2}(1:f1-1));
        date_feature1{f,7} = date_part{i,4}(1)+f1;
		date_feature1{f,8} = 5;
		f = f+1;

		f1 = floor((1+rand(1))*date_part{i,3}/2);
		f2 = date_part{i,1}(1:f1);
		date_feature1{f,1} = length(date_part{i,1});
		date_feature1{f,2} = f1;
    	f3 = date_part{i-1,1};
    	f4 = 1:length(f3);
    	f5 = regress(f3,[f4'*0.1, ones(length(date_part{i-1,1}),1)]);
    	date_feature1{f,3} = f5(1);
    	f6 = 1:length(f2);
    	f7 = regress(f2,[f6'*0.1,ones(length(f2),1)]);
    	date_feature1{f,4} =f7(1);
    	date_feature1{f,5} = mean(date_part{i-1,2});
		date_feature1{f,6} = mean(date_part{i,2}(1:f1-1));
        date_feature1{f,7} = date_part{i,4}(1)+f1;
		date_feature1{f,8} = 5;
		f = f+1;
	elseif date_part{i,3} >= 20
		f1 = 2+floor(rand(1)*(date_part{i,3}/3-2));
		f2 = date_part{i,1}(1:f1);
		date_feature1{f,1} = length(date_part{i,1});
		date_feature1{f,2} = f1;
    	f3 = date_part{i-1,1};
    	f4 = 1:length(f3);
    	f5 = regress(f3,[f4'*0.1, ones(length(date_part{i-1,1}),1)]);
    	date_feature1{f,3} = f5(1);
    	f6 = 1:length(f2);
    	f7 = regress(f2,[f6'*0.1,ones(length(f2),1)]);
    	date_feature1{f,4} =f7(1);
    	date_feature1{f,5} = mean(date_part{i-1,2});
		date_feature1{f,6} = mean(date_part{i,2}(1:f1-1));
        date_feature1{f,7} = date_part{i,4}(1)+f1;
		date_feature1{f,8} = 5;
		f = f+1;

		f1 = floor((rand(1)+1)*date_part{i,3}/3);
		f2 = date_part{i,1}(1:f1);
		date_feature1{f,1} = length(date_part{i,1});
		date_feature1{f,2} = f1;
    	f3 = date_part{i-1,1};
    	f4 = 1:length(f3);
    	f5 = regress(f3,[f4'*0.1, ones(length(date_part{i-1,1}),1)]);
    	date_feature1{f,3} = f5(1);
    	f6 = 1:length(f2);
    	f7 = regress(f2,[f6'*0.1,ones(length(f2),1)]);
    	date_feature1{f,4} =f7(1);
    	date_feature1{f,5} = mean(date_part{i-1,2});
		date_feature1{f,6} = mean(date_part{i,2}(1:f1-1));
        date_feature1{f,7} = date_part{i,4}(1)+f1;
		date_feature1{f,8} = 5;
		f = f+1;

		f1 = floor((rand(1)+2)*date_part{i,3}/3);
		f2 = date_part{i,1}(1:f1);
		date_feature1{f,1} = length(date_part{i,1});
		date_feature1{f,2} = f1;
    	f3 = date_part{i-1,1};
    	f4 = 1:length(f3);
    	f5 = regress(f3,[f4'*0.1, ones(length(date_part{i-1,1}),1)]);
    	date_feature1{f,3} = f5(1);
    	f6 = 1:length(f2);
    	f7 = regress(f2,[f6'*0.1,ones(length(f2),1)]);
    	date_feature1{f,4} =f7(1);
    	date_feature1{f,5} = mean(date_part{i-1,2});
		date_feature1{f,6} = mean(date_part{i,2}(1:f1-1));
        date_feature1{f,7} = date_part{i,4}(1)+f1;
		date_feature1{f,8} = 5;
		f = f+1;
	end
end
