%% 本脚本需要配合test3文件使用
%% 每段长度分布
c = diff(extre_index(:, 1));  % 在这里查看极大值间隔
histogram(c)

%% 测试删除效果
t =extre_index;
%t = a_index;
plot(date_s);
hold on
plot(t(:,1),date_s(t(:,1)),'r');
hold off
%% 斜率编码测试
a = date_part{2,1};
a1 = 1:length(a);
%cc =a1'*0.1;
%cc1 = ones(length(a),1);
b = regress(a,[a1'*0.1,ones(length(a),1)]);
plot(a1,a);
