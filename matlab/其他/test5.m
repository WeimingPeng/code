%% ���ļ������test3��
%% P1 ��������������,������Ӧ��Ϊ5��
p1 = cell2mat(date_feature);
p2 = cell2mat(date_feature1);
p3 = [p1;p2];
date_feat = sortrows(p3,7);

%% P2 �������ݴ�����Ӧ��Ϊ����
date_feat1 = ones(size(date_feat,1),7);
date_feat1(:,1:6) = date_feat(:,1:6);
for i = 1:size(date_feat,1)
	if date_feat(i,8) == 5
		date_feat1(i,7) = 2;
	end
end

%% P3 ����ѵ���Ͳ������һ��
% ��������
p1 = floor(size(date_feat1,1)*0.8);
p2 = date_feat1(1:p1,:); % ѵ������ 
p3 = date_feat1(p1+1:end,:); % ��������
% ��һ��
[p4,re] = mapminmax(p2(:,1:6)',0,1);
p_train = [p4',p2(:,7)];
p5 = mapminmax('apply',p3(:,1:6)',re);
p_test = p5';
% �����ģ��ʹ�õ���Ϊ�ṹ����ʽ��������Ҫ��ѡ��Ϊ����Ϊѵ�����������忴�鿴���� 

%% P4  ��֤ѵ�����
yfit = trainedModel.predictFcn(p_test);
% ����������ȷ��
t_right = 0;
for i = 1:length(yfit)
	if yfit(i)-p3(i,7) == 0
		t_right = t_right+1;
	end
end
right_rate = t_right/length(yfit);
right_rate

% ͻ���Ԥ����ȷ��
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
