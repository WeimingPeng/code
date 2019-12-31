%% 本文件处理接test3、
%% P1 特征数据排序处理,本处响应类为5类
p1 = cell2mat(date_feature);
p2 = cell2mat(date_feature1);
p3 = [p1;p2];
date_feat = sortrows(p3,7);

%% P2 特征数据处理，响应类为两类
date_feat1 = ones(size(date_feat,1),7);
date_feat1(:,1:6) = date_feat(:,1:6);
for i = 1:size(date_feat,1)
	if date_feat(i,8) == 5
		date_feat1(i,7) = 2;
	end
end

%% P3 划分训练和测试与归一化
% 划分数据
p1 = floor(size(date_feat1,1)*0.8);
p2 = date_feat1(1:p1,:); % 训练数据 
p3 = date_feat1(p1+1:end,:); % 测试数据
% 归一化
[p4,re] = mapminmax(p2(:,1:6)',0,1);
p_train = [p4',p2(:,7)];
p5 = mapminmax('apply',p3(:,1:6)',re);
p_test = p5';
% 这里的模型使用导出为结构的形式，如有需要可选择为导出为训练函数。具体看查看帮助 

%% P4  验证训练结果
yfit = trainedModel.predictFcn(p_test);
% 总体数据正确率
t_right = 0;
for i = 1:length(yfit)
	if yfit(i)-p3(i,7) == 0
		t_right = t_right+1;
	end
end
right_rate = t_right/length(yfit);
right_rate

% 突变点预测正确率
t_sum = 0;
t_right1 = 0;
for i = 1:length(yfit)
	if p3(i,7) == 1
		t_sum = t_sum+1;
		if yfit(i) == 1
			t_right1 = t_right1+1;
		end
	end
end
right_rate1 = t_right1/t_sum;
right_rate1
