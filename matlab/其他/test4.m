%% ���ű���Ҫ���test3�ļ�ʹ��
%% ÿ�γ��ȷֲ�
c = diff(extre_index(:, 1));  % ������鿴����ֵ���
histogram(c)

%% ����ɾ��Ч��
t =extre_index;
%t = a_index;
plot(date_s);
hold on
plot(t(:,1),date_s(t(:,1)),'r');
hold off
%% б�ʱ������
a = date_part{2,1};
a1 = 1:length(a);
%cc =a1'*0.1;
%cc1 = ones(length(a),1);
b = regress(a,[a1'*0.1,ones(length(a),1)]);
plot(a1,a);
