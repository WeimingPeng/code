%% ��peng�е�����extremum����
%�����ݰ�ƽ����ĵ㻮�ֳɶ�
date = extremum_index;
h_part = cell(length(date(1,:)),2);
h_part{1,1} = date(1,:);
h_part{1,2} = date_ahm(1:date(1,1));
for i = 2:length(date(:,1))
    h_part{i,1} = date(i,:);
    h_part{i,2} = date_ahm(date(i-1,1):date(i,1));
end

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
        



