% 本代码需要在运行datapick后运行
%% 合并正例与负例
% 正例
for i = 2:size(postive_feat,1)
    postive_set(i-1,1:8) = [postive_feat(i-1,1:4),postive_feat(i,1:4)];
    postive_set(i-1,9) = postive_feat(i,5);
end
% 负例
for i = 1:size(p7_out,1)
    p1 = find(extra_index1<p7_out(i,5));
    p2 = p1(end);
    negative_set(i,1:8) = [postive_feat(p2,1:4),p7_out(i,1:4)];
    negative_set(i,9) = p7_out(i,5);
end
negative_set(:,10) = 1; %负例标签为1
postive_set(:,10) = 2; %正例标签为2
set_data = sortrows([postive_set; negative_set],9);
train_set = set_data(set_data(:,9) <= 75000,:);
test_set = set_data(set_data(:,9) > 75000,:);

%% 随机森林分类训练
p_train = train_set(:,1:8);
t_train = train_set(:,10);
p_test = test_set(:,1:8);
t_test = test_set(:,10);
% 随机森林代码为网络下载来源于
% http://code.google.com/p/randomforest-matlab/downloads/detail?name=Windows-Precompiled-RF_MexStandalone-v0.02-.zip&can=2&q=
model = classRF_train(p_train,t_train);

% 效果测试
[T_sim,votes] = classRF_predict(p_test,model);