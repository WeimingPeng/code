%% ���ݶ���
clc
clear
% ��������
date_ah = csvread("D:\date_code\date\B_Hour.csv",1,1);
date_ahm =mean(date_ah(1:74000, 2:3), 2); % �����м���ֵ�뼫Сֵ֮���ƽ��
date_ahtest = mean(date_ah(74001:end,2:3),2);

%% ������ȡ
% ѵ��������ȡ
p_train1 = pick(date_ahm,date_ahm); % ����������ȡ
p_test1 = pick(date_ahtest,date_ahtest); % ������ȡ��Ϊֻ�������ݣ���һ���ͷ���һ�������֣������ڴ˴�

%% ����ѵ����
% ��һ��
% step1 ѵ������һ��
[THR,SORH,KEEPAPP,CRIT]=ddencmp('den','wp',date_ahm');
xx=wdencmp('gbl',date_ahm','db6',3,THR,SORH,KEEPAPP);
date_s=smoothdata(xx,'movmean',4); % ʹ��������ƽ����С�ٽ�λ�õļ�ֵ
date_f = diff(date_s);
[date_diff, re] = mapminmax(date_f,0,1);
t1 =randperm(size(p_train1,1));
train_date = p_train1(t1,:); % ѵ�������ᱻ����
p_train11 = train_date(:,5)';
for i = 1:size(p_train11,2)
    t2 = mapminmax('apply',p_train11{1,i}',re);
    p_train{1,i} =t2';
end
% step2 ���Լ���һ��
for i = 1:size(p_test1, 1)
    t2 = mapminmax('apply',p_test1{i,5}',re);
    p_test{1,i} =t2';
end
% ��дѵ�����
for i = 1:size(train_date,1)
    switch train_date{i, 6}
        case 1
            t_train{1,i} = [1;0];
        case 2
            t_train{1,i} = [0;1];
    end
end


%% ����������
net = newrb(p_train,t_train, 0.2, 0.75);

%% Ԥ��
t_pr =sim(net,p_test);

%% Ч������
%�����ת��Ϊ���
for i =1:size(t_pr,2)
    p_out(i,1) = find(t_pr{1,i} ==max(t_pr{1,i}));
end
 % ��ȷ��
 t_right = 0;
 for i =1:size(t_pr,2)
    if p_out(i,1)-p_test1{i,6}==0
        t_right = t_right+1;
    end
end
right_ridio = t_right/size(t_pr,2);

