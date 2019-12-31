function p_train =pick(date_a1,date_ahm)
    %% ����������������ȡ����������ǹ�һ��������ݣ�ԭʼ����
    % �����k1��k2���ֱ������ǰȡ�����ݺ����ȡ�����ݵĳ���
    k1 = 12;
    k2 = 5;
    %% ƽ��
    [THR,SORH,KEEPAPP,CRIT]=ddencmp('den','wp',date_a1');
    xx=wdencmp('gbl',date_a1','db6',3,THR,SORH,KEEPAPP);
    date_s=smoothdata(xx,'movmean',4); % ʹ��������ƽ����С�ٽ�λ�õļ�ֵ

    %% ��ֵ��������ȡ
    datemin_index(:, 1) = find(diff(sign(diff(date_s)))>0)+1; %min
    datemin_index(:, 2) = 1;
    datemax_index(:, 1) = find(diff(sign(diff(date_s)))<0)+1;
    datemax_index(:, 2) = 2;
    extremum_index = sortrows([datemin_index; datemax_index]);
    % c = diff(extremum_index(:, 1));  % ������鿴����ֵ����ķֲ�����

    %% ��ĳһ���ݳ�����ȡ���ݡ�Ŀǰ��λһ��
    % ���潫���ٽ��������ݽ��д��� ;������ȡ����Ϊ�ٽֵ����������
    date1 = extremum_index(3:end,:);
    apart = cell(length(date1(:,1)),4);
    for j = 1:length(date1(:,1))
        apart{j,1} = date1(j,2);
        apart{j,2} = date1(j,1);
        apart{j,3} = date_ahm(date1(j,1)+1:date1(j,1)+k2);
        apart{j,4} = date_ahm(date1(j,1)-k1:date1(j,1));
        apart{j,5} = 1;
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
        for j = 1:(date_1(i,2)-date_1(i,1)-k1-k2+1)
            bpart{t2,1} = date_1(i,3); % ��Ǵ���1���������䣬��Ǵ���2Ϊ�½�����
            bpart{t2,2} = date_1(i,1)+k1+j;
            bpart{t2,3} = date_ahm(date_1(i,1)+j+k1+1 : date_1(i,1)+j+k1+k2);
            bpart{t2,4} = date_ahm(date_1(i,1)+j : date_1(i,1)+j+k1);
            bpart{t2,5} =2;
            t2 =t2+1;
        end
    end
    %% ��������ȡ�����������ݺϲ�
    tt = size(apart,1)+size(bpart,1);
    date_part = cell(tt, 5);
    date_part(1:size(apart,1),:) =apart;
    date_part(size(apart,1)+1:end,:) =bpart;
    
    %% ����ȡ���ݴ���
    p_train = cell(size(date_part,1),6);
    p_train(:,1:4) =date_part(:,1:4);
    for i = 1:size(date_part,1)
        a = smoothdata(date_part{i,4},'movmean',3);
        p_train{i,5} = diff(a);
    end
    p_train(:,6) = date_part(:,5);
    
