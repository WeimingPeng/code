function p_train =pick(date_a1,date_ahm)
    %% 本函数用于数据提取，传入参数是归一化后的数据，原始数据
    % 下面的k1和k2将分别控制向前取的数据和向后取的数据的长度
    k1 = 12;
    k2 = 5;
    %% 平滑
    [THR,SORH,KEEPAPP,CRIT]=ddencmp('den','wp',date_a1');
    xx=wdencmp('gbl',date_a1','db6',3,THR,SORH,KEEPAPP);
    date_s=smoothdata(xx,'movmean',4); % 使用了两次平滑减小临近位置的极值

    %% 极值查找与提取
    datemin_index(:, 1) = find(diff(sign(diff(date_s)))>0)+1; %min
    datemin_index(:, 2) = 1;
    datemax_index(:, 1) = find(diff(sign(diff(date_s)))<0)+1;
    datemax_index(:, 2) = 2;
    extremum_index = sortrows([datemin_index; datemax_index]);
    % c = diff(extremum_index(:, 1));  % 在这里查看极大值间隔的分布特性

    %% 按某一数据长度提取数据。目前定位一天
    % 下面将对临界收敛数据进行处理 ;即所提取数据为临街点后收敛数据
    date1 = extremum_index(3:end,:);
    apart = cell(length(date1(:,1)),4);
    for j = 1:length(date1(:,1))
        apart{j,1} = date1(j,2);
        apart{j,2} = date1(j,1);
        apart{j,3} = date_ahm(date1(j,1)+1:date1(j,1)+k2);
        apart{j,4} = date_ahm(date1(j,1)-k1:date1(j,1));
        apart{j,5} = 1;
    end
    % 下面将提取发散数据
    % step1 筛选出连续发散区间，滚动取值
    % 提取连续变化说我数据
    t1 = 1;
    for i = 2:length(date1(:,1))
        if date1(i,1) - date1(i-1,1)>30
            date_1(t1, :) =[date1(i-1), date1(i,1), date1(i-1,2)];
            t1 = t1+1;
        end
    end
    % 提取目标数据
    t2 =1;
    for i = 1: length(date_1(:,1))
        for j = 1:(date_1(i,2)-date_1(i,1)-k1-k2+1)
            bpart{t2,1} = date_1(i,3); % 标记代号1是上升区间，标记代号2为下降区间
            bpart{t2,2} = date_1(i,1)+k1+j;
            bpart{t2,3} = date_ahm(date_1(i,1)+j+k1+1 : date_1(i,1)+j+k1+k2);
            bpart{t2,4} = date_ahm(date_1(i,1)+j : date_1(i,1)+j+k1);
            bpart{t2,5} =2;
            t2 =t2+1;
        end
    end
    %% 将上面提取的两部分数据合并
    tt = size(apart,1)+size(bpart,1);
    date_part = cell(tt, 5);
    date_part(1:size(apart,1),:) =apart;
    date_part(size(apart,1)+1:end,:) =bpart;
    
    %% 对提取数据处理
    p_train = cell(size(date_part,1),6);
    p_train(:,1:4) =date_part(:,1:4);
    for i = 1:size(date_part,1)
        a = smoothdata(date_part{i,4},'movmean',3);
        p_train{i,5} = diff(a);
    end
    p_train(:,6) = date_part(:,5);
    
