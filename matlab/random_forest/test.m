%% 本小节用于验证data_predict所产生结果正确率
%负例标签为1、正例标签为2
r1 = 0;t1=0;r=0;t=0;
for i = 1:length(t_test)
    if t_test(i) == 1
        t = t+1; 
        if T_sim(i) == t_test(i) 
            r = r+1;
        end
    end
    if t_test(i) == 2
        t1 = t1+1; 
        if T_sim(i) == t_test(i) 
            r1 = r1+1;
        end
    end
end
right_ratio1 = r/t;
right_ratio2 = r1/t1;
right_ratio = (r+r1)/(t+t1);
% 结果输出
disp(['负例预测正确率：' num2str(right_ratio1)])
disp(['正例预测正确率：' num2str(right_ratio2)])
disp(['总预测正确率：' num2str(right_ratio)])
r11 = 0;t11=0;r0=0;t0=0;
for i = 1:length(t_test)
    if T_sim(i) == 1
        t0 = t0+1; 
        if T_sim(i) == t_test(i) 
            r0 = r0+1;
        end
    end
    if T_sim(i) == 2
        t11 = t11+1; 
        if T_sim(i) == t_test(i) 
            r11 = r11+1;
        end
    end
end
right_ratio11 = r0/t0;
right_ratio22 = r11/t11;
right_ratio0 = (r0+r1)/(t0+t1);
% 结果输出
disp(['预测为负例的正确率：' num2str(right_ratio11)])
disp(['预测为正例的正确率：' num2str(right_ratio22)])
disp(['总预测正确率：' num2str(right_ratio0)])
%% 绘图
p1 = 1;
p2 = 1;
test_i = find(set_data(:,9) > 75000);
test_index = set_data(test_i,9);
for i = 1:length(test_index)
    if T_sim(i) == t_test(i)
        index_1(p1) = test_index(i);
        p1 = p1+1;
    else
        index_2(p2) = test_index(i);
        p2 = p2+1;
    end
end
p3 = 1; p4 =1;p5 = 1; p6 =1;
for i = 1:length(test_index)
    if t_test(i) == 1
        if T_sim(i) == t_test(i)
            index_3(p3) = test_index(i);
            p3 = p3+1;
        else
            index_4(p4) = test_index(i);
            p4 = p4+1;      
        end
    elseif  t_test(i) == 2
        if T_sim(i) == t_test(i)
            index_5(p5) = test_index(i);
            p5 = p5+1;
        else
            index_6(p6) = test_index(i);
            p6 = p6+1;      
        end
    end
end
% 平滑总体效果
figure(1)
plot(data)
hold on
plot(index_1,data(index_1),'o')
plot(index_2,data(index_2),'r*')
hold off
% 未平滑数据效果
figure(2)
plot(data_in)
hold on
plot(index_1,data_in(index_1),'ro')
plot(index_2,data_in(index_2),'b*')
hold off
%标记为负例的效果
figure(3)
plot(data_in)
hold on
plot(index_3,data_in(index_3),'bo')
plot(index_4,data_in(index_4),'r*')
hold off
%标记为正例的效果
figure(4)
plot(data_in)
hold on
plot(index_5,data_in(index_5),'ro')
plot(index_6,data_in(index_6),'b*')
hold off