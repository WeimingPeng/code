%% ����������֤
i=1;
ccc =0;
while i<length(extra_index)-3
    if extra_index(i+3) - extra_index(i)<4
        ccc = ccc+1;
        c_mark(ccc) = i;
    end
    i=i+1;
end

%% ɸѡ�����֤
plot(data,'b')
hold on
plot(extra_index1, data(extra_index1),'r')
hold off

%% ���䳤����֤
l = diff(extra_index1);
