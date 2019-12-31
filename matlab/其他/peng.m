%% ���ݶ���
clc
clear
% ��������ͳһΪ2004��4��5��0��0ʱ0��
date_ah = csvread("D:\date_code\date\A_Hour.csv",4,1);
date_ad = csvread("D:\date_code\date\A_Day.csv",2,1);
date_ahm =mean(date_ah(1:74000, 2:3), 2); % �����м���ֵ�뼫Сֵ֮���ƽ��
date_ahtest = mean(date_ah(74001:end,2:3));

%% ƽ��
[THR,SORH,KEEPAPP,CRIT]=ddencmp('den','wp',date_ahm);
xx=wdencmp('gbl',date_ahm,'db6',3,THR,SORH,KEEPAPP);
date_s=smoothdata(xx,'movmean',4); % ʹ��������ƽ����С�ٽ�λ�õļ�ֵ
% date_s = date_ahm;

%% ��ֵ��������ȡ
datemin_index(:, 1) = find(diff(sign(diff(date_s)))>0)+1;%min
datemin_index(:, 2) = 1;
datemax_index(:, 1) = find(diff(sign(diff(date_s)))<0)+1;
datemax_index(:, 2) = 2;
extremum_index = sortrows([datemin_index; datemax_index]);
c = diff(extremum_index(:, 1));  % ������鿴����ֵ����ķֲ�����

%% ���潫���뷨2  ��ֵ����ϲ�
%% ���ÿ������ľ�ֵ
date_index =extremum_index;
date_index(1,3) = mean(date_ahm(1:extremum_index(1,1)));
for i = 2:length(extremum_index(:,1))
    date_index(i,3) = mean(date_ahm(extremum_index(i-1,1):extremum_index(i,1)));
end

h_index = date_index;
t1 = 0;
for i = 2:length(date_index(:,1))
    if abs(date_index(i,3)-date_index(i-1,3))<0.0025
        h_index(i-t1,:) = [];
        t1 = t1+1;
    end
end

%% ɾ������
figure(1)
plot(date_s)
hold on
plot(h_index(:,1),date_s(h_index(:,1)),'r')
hold off

%% ����ֲ�����
clear tr
close all
tr =diff(date_index(:,3));
figure(1)
histogram(tr)
figure(2)
p = capaplot(tr, [-0.0035,0.0035]);



%%  ���潫���뷨һ
%% �ֶΣ����һЩ����

a = extremum_index;
new_index = extremum_index;
ti=0;
if a(1,1)-1<15
    new_index(1,:) = [];
    ti =ti+1;
end
for i = 2: length(a(:,1))-2
    if  a(i,1) - a(i-1,1) <15
        if date_s(a(i,1)) - date_s(a(i-1, 1)) <0.015 
            new_index(i-ti,:) = [];
            ti=ti+1;
        end
    end   
end

%% ɾ������
figure(1)
plot(date_s)
hold on
% a = 1:120:74000;
% plot(a,date_s(a),'o')
plot(new_index(:,1),date_s(new_index(:,1)),'r')
hold off

%% ����
t=0;
for i = 1:length(c)
    if c(i)>15
        t = t+1;
    end
end


%% ��ͼ����
subplot(2,1,1);
plot(date_ah(1:2400));
subplot(2,1,2);
hold on
plot(date_ad(1:100,2),'r');
plot(date_ad(1:100,3),'k');

%% emd�ֽⳢ��
date_test = date_ad(1:100,2);
emd(date_test)

%% ��һС��
