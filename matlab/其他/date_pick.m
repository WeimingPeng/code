%% ��������
clc
clear
% ����������ȡ�������Сʱ���ݽ�����ȡ

date_1 = csvread('D:\date_code\date\B_Hour.csv',1,1);
date_2 = date_1(:, 2:3)'; % ��ȡ���ֵ����Сֵ����
date_3 = date_2';
date_0 = (mean(date_2))';%ȡСʱ�����ֵ��Сֵƽ��
 r = floor(date_0*100)/100;

%% ƽ��
[THR,SORH,KEEPAPP,CRIT]=ddencmp('den','wp',date_bhm);
xx=wdencmp('gbl',date_bhm,'db6',3,THR,SORH,KEEPAPP);
date_s=smoothdata(xx,'movmean',4); %ʹ������ƽ����С�ٽ�λ�õļ�ֵ

%% ��ֵ��������ȡ
datemin_index(:, 1) = find(diff(sign(diff(date_s)))>0)+1;%min
datemin_index(:, 2) = 1;
datemax_index(:, 1) = find(diff(sign(diff(date_s)))<0)+1;
datemax_index(:, 2) = 2;
extremum_index = sortrows([datemin_index; datemax_index]);
 c = diff(extremum_index(:, 1)); % ������鿴����ֵ����ķֲ�����

%%
dd=[];
dd((1:length(b)))=b';
dd((length(b)+1):(length(b)+length(c)))=c';
d=dd';
[e,f]=sort(d);%e��f�ϲ���һ��������һ���ж�e��ֵΪ����ֵ���Ǽ�Сֵ
test_p(1,1)=1;test_p(1,2)=2;test_p(1,3)=x(1);
test_p(2:length(e)+1,1)=e;
for i=1:length(e)
    if f(i)<=(length(b))
        test_p(i+1,2)=1;%min
    else
        test_p(i+1,2)=2;
    end
end
test_p(2:length(e)+1,3)=x(test_p(2:length(e)+1,1));
y=diff(test_p(:,3));%�����õ���һ�������ȥ����Ҫ�ļ�ֵ�㣬����һ�������⣬��ϸ������ʼ�����paper1
ya=abs(y);
s=2;
yb=[];yb(1,:)=test_p(1,:);
for i=1:length(ya)
    if ya(i)>0.008 % ������ֵ���Ӧ�ô������ֵ  
        yb(s,:)=test_p(i+1,:);
        s=s+1;
    end
end
yc=[];
yd=1;
for i=1:length(yb(:,1))-1
    if yb(i,2)==2
        if yb(i+1,2)==1
            yc(yd,:)=yb(i,:);
             yd=yd+1;
        end
    elseif yb(i,2)==1
       if yb(i+1,2)==2
            yc(yd,:)=yb(i,:);
             yd=yd+1;
        end
    end
end %yc Ϊ������ȡ�Ľ��
%% ��ֵ�����ݼ�������ȡ
fss=1;
fsr=1;
sd=[];
t=1;
q=1;
ws=[];
test_ss=[];
fsd=length(test_p(:,1));
for k=1:length(test_p(:,1))-1
     if test_p(k,2)==1
         ws(q,1)=test_p(k,1);
         ws(q,2)=1;
         q=q+1;
         ye=sign(x(test_p(k,1):end)-x(test_p(k,1)));
        qq=find(diff(ye)<0)+test_p(k,1);
        for g=1:length(qq)
            for i=1:(fsd-1)
                j=test_p(i,1);
                u=test_p(i+1,1);
                if  (j<qq(g))&&(qq(g)<u)
                       if ((qq(g)-j)>12)&&(qq(g)+6<93338)&&((u-qq(g))>5)
                           %�˴�������ȡ�м����Ҫ����ȡ��Ҫ��
                           if test_p(i,2)==2%max��min�м�
                                test_ss(t,1)=qq(g);
                                ws(q,1)=qq(g);
                                ws(q,2)=3;
                                test_ss(t,2)=3;
                                test_ss(t,3:15)=a(qq(g)-6:qq(g)+6);
                                t=t+1;
                                q=q+1;
                           end
                       end
                end
            end
        end
     elseif test_p(k,2)==2
         ws(q,1)=test_p(k,1);
         ws(q,2)=2;
         q=q+1;
         ye=sign(x(test_p(k,1):end)-x(test_p(k,1)));
        qq=find(diff(ye)>0)+test_p(k,1);
        for g=1:length(qq)
            for i=1:(fsd-1)
                j=test_p(i,1);
                u=test_p(i+1,1);
                if  (j<qq(g))&&(qq(g)<u)
                    fss=fss+1;
                       if ((qq(g)-j)>12)&&(qq(g)+6<93338)&&((u-qq(g))>5)
                           if test_p(i,2)==1%min��max�м�
                                test_ss(t,1)=qq(g);
                                ws(q,1)=qq(g);
                                ws(q,2)=4;
                                test_ss(t,2)=4;
                                test_ss(t,3:15)=a(qq(g)-6:qq(g)+6);
                                t=t+1;
                                q=q+1;
                           end
                      end
                end
            end
        end
     end
end
test_w=[];
test_w(1:length(yc),:)=yc(:,1:2);
test_w((length(yc)+1):(length(yc)+length(test_ss(:,1))),:)=test_ss(:,1:2);
test_ws=sortrows(test_w,1);%��ָ����һ������
test_wss=unique(test_ws,'rows') ;%ȥ���ظ���
xa=[];xb=[];xc=1;xd=1;
for i=1:length(test_wss(:,1))
    if test_wss(i,2)==3
        xa(xc)=test_wss(i,1);
        xc=xc+1;
    elseif test_wss(i,2)==4
        xb(xd)=test_wss(i,1);
        xd=xd+1;
    end
end %test_wssΪ������ݣ����ݴ�ȡ��ʱ��������Ϊdate_index


    
