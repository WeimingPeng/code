%% ���ݶ���
clc
clear
% ��������ͳһΪ2004��4��5��0��0ʱ0��
date_ah = csvread("D:\date_code\date\A_Hour.csv",1,1);
date_ahm =mean(date_ah(1:74000, 2:3), 2); % �����м���ֵ�뼫Сֵ֮���ƽ��
date_ahtest = mean(date_ah(74001:end,2:3),2);

%% ��һ��
[date_a1, re] = mapminmax(date_ahm',0,1); % �������ŵ�0��1֮��
date_a2 = mapminmax("apply",date_ahtest',re); %�Բ������ݰ�ѵ�����ݹ�һ�� 

%% ���潫��ѵ����������ȡ
%% ƽ��
[THR,SORH,KEEPAPP,CRIT]=ddencmp('den','wp',date_a1');
xx=wdencmp('gbl',date_a1','db6',3,THR,SORH,KEEPAPP);
date_s=smoothdata(xx,'movmean',4); % ʹ��������ƽ����С�ٽ�λ�õļ�ֵ

%% ��ֵ��������ȡ
datemin_index(:, 1) = find(diff(sign(diff(date_s)))>0)+1;%min
datemin_index(:, 2) = 1;
datemax_index(:, 1) = find(diff(sign(diff(date_s)))<0)+1;
datemax_index(:, 2) = 2;
extremum_index = sortrows([datemin_index; datemax_index]);
c = diff(extremum_index(:, 1));  % ������鿴����ֵ����ķֲ�����

%% ��ĳһ���ݳ�����ȡ���ݡ�Ŀǰ��λһ��
% ���潫���ٽ��������ݽ��д��� ;������ȡ����Ϊ�ٽֵ����������
date1 = extremum_index(3:end,:);
apart = cell(length(date1(:,1)),2);
for j = 1:length(date1(:,1))
    apart{j,1} = date1(j,:);
    apart{j,2} = date_ahm(date1(j,1)-24:date1(j,1));
end
% ���潫��ȡ��ɢ����
% step1 ɸѡ��������ɢ���䣬����ȡֵ
% ��ȡ�����仯˵������
t1 = 1;
for i = 2:length(date1(:,1))
    if date1(i,1) - date1(i-1,1)>30
        date_1(t1, :) =[date1(i-1), date1(i,1), date1(i-1,2)];
        t1 = t1+1;
    end
end
% ��ȡĿ������
t2 =1;
for i = 1: length(date_1(:,1))
    for j = 1:(date_1(i,2)-date_1(i,1)-29)
        bpart{t2,1} = date_1(i,3); % ��Ǵ���1���������䣬��Ǵ���2Ϊ�½�����
        bpart{t2,2} = date_1(i,1)+24+j;
        bpart{t2,3} = date_ahm(date_1(i,1)+j : date_1(i,1)+j+24);
        t2 =t2+1;
    end
end
