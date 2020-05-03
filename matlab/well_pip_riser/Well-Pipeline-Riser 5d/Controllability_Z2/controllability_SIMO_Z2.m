clc
clear all
warning ('off')
Current_Folder = cd('../WPR_model');
%% **** System Inputs ****
par = v2_new_5d_parameters(); %All parameters defined in one file
z1 = 0.1;
z2 = 0.1;
global u0;
u0 = [z1;z2;0;0];
%% ******** Initializing *************
disp('Initializing...')
tic
[x0,y0,par] =  v2_new_5d_initialize(u0,par);
toc
disp('Initialization completed.');
disp('#########################');
%% ******** Linearization *************
[A,B,C,D]=v2_new_5d_linmodl_num(x0,u0,par);%numerical
sys =ss(A,B,C,D);
set(sys,'inputname',{'Z1','Z2','d1','d2'},'outputname',{'P_bh','P_wh','W_in','P_in','P_rb','DP_r','P_t','Q_out','W_out','Rho_t','Alpha_L'},'statename',{'m_Gw','m_Gp','m_Lp','m_Gr','m_Lr'});
%% ****Scaling****
%y output( P1,P2,W)
%d disturbance(wG_in, wL_in)
%u input (Z)

[De,Dd,Du] = scaling(z1,z2);

Bs1 = B(:,1:2)*Du;
Bs2 = B(:,3:4)*Dd;
Bs = [Bs1 Bs2];

Cs = (De^-1)*C;

Ds1 = (De^-1)*D(:,1:2)*Du;
Ds2 = (De^-1)*D(:,3:4)*Dd;
Ds = [Ds1 Ds2];

cd(Current_Folder);

sys =(ss(A,Bs,Cs,Ds));
set(sys,'inputname',{'Z1','Z2','d1','d2'},'outputname',{'P_bh','P_wh','W_in','P_in','P_rb','DP_r','P_t','Q_out','W_out','Rho_t','Alpha_L'},'statename',{'m_Gw','m_Gp','m_Lp','m_Gr','m_Lr'});
sys = minreal(sys);
Gu = tf(sys({'P_bh','P_wh','W_in','P_in','P_rb','DP_r','P_t','Q_out','W_out','Rho_t','Alpha_L'},'Z2'));
Gd = tf(sys({'P_bh','P_wh','W_in','P_in','P_rb','DP_r','P_t','Q_out','W_out','Rho_t','Alpha_L'},{'d1','d2'}));
%% ****Tranfer function and steady state values ****
G1 =(Gu('P_bh','Z2')); G2 =(Gu('P_wh','Z2')); G3 =(Gu('W_in','Z2')); G4 =(Gu('P_in','Z2'));
G5 =(Gu('P_rb','Z2')); G6 =(Gu('DP_r','Z2')); G7 =(Gu('P_t','Z2')); G8 =(Gu('Q_out','Z2'));
G9 =(Gu('W_out','Z2')); G10 =(Gu('Rho_t','Z2')); G11 =(Gu('Alpha_L','Z2'));

%Have to have the disturbance matrix Gd for each output,the gas and liquid flow rate into the riser are going
%to be analysed as two separately disturbances.
Gd11=(Gd('P_bh','d1')); Gd12=(Gd('P_bh','d2'));
Gd21=(Gd('P_wh','d1')); Gd22=(Gd('P_wh','d2'));
Gd31=(Gd('W_in','d1'));  Gd32=(Gd('W_in','d2'));
Gd41=(Gd('P_in','d1'));  Gd42=(Gd('P_in','d2'));
Gd51=(Gd('P_rb','d1'));  Gd52=(Gd('P_rb','d2'));
Gd61=(Gd('DP_r','d1'));  Gd62=(Gd('DP_r','d2'));
Gd71=(Gd('P_t','d1'));  Gd72=(Gd('P_t','d2'));
Gd81=(Gd('Q_out','d1'));  Gd82=(Gd('Q_out','d2'));
Gd91=(Gd('W_out','d1'));  Gd92=(Gd('W_out','d2'));
Gd101=(Gd('Rho_t','d1'));  Gd102=(Gd('Rho_t','d2'));
Gd111=(Gd('Alpha_L','d1'));  Gd112=(Gd('Alpha_L','d2'));

%Both G1, G2 and G3 have the same poles.
disp('poles:')
p=pole(G1);
disp(p)
disp('zeros of P_bh:')
z1=zero(G1);
disp(z1)
disp('zeros of P_wh:')
z2=zero(G2);
disp(z2)
disp('zeros of W_in:')
z3=zero(G3);
disp(z3)

disp('zeros of P_in:')
z4=zero(G4);
disp(z4)

disp('zeros of P_rb:')
z5=zero(G5);
disp(z5)

disp('zeros of DP_r:')
z6=zero(G6);
disp(z6)

disp('zeros of P_t:')
z7=zero(G7);
disp(z7)

disp('zeros of Q_out:')
z8=zero(G8);
disp(z8)

disp('zeros of W_out:')
z9=zero(G9);
disp(z9)

disp('zeros of Rho_t:')
z10=zero(G10);
disp(z10)

disp('zeros of Alpha_L:')
z11=zero(G11);
disp(z11)

% Stable phase, minimum phase and minimum stable phase of all TF, Functions
% are made to performed this.
[Gs1 Gm1 Gms1]=gmsi(G1);
[Gs2 Gm2 Gms2]=gmsi(G2);
[Gs3 Gm3 Gms3]=gmsi(G3);
[Gs4 Gm4 Gms4]=gmsi(G4);
[Gs5 Gm5 Gms5]=gmsi(G5);
[Gs6 Gm6 Gms6]=gmsi(G6);
[Gs7 Gm7 Gms7]=gmsi(G7);
[Gs8 Gm8 Gms8]=gmsi(G8);
[Gs9 Gm9 Gms9]=gmsi(G9);
[Gs10 Gm10 Gms10]=gmsi(G10);
[Gs11 Gm11 Gms11]=gmsi(G11);

[Gds11 Gdm11 Gdms11]=gmsdij(Gd11);
[Gds12 Gdm12 Gdms12]=gmsdij(Gd12);

[Gds21 Gdm21 Gdms21]=gmsdij(Gd21);
[Gds22 Gdm22 Gdms22]=gmsdij(Gd22);

[Gds31 Gdm31 Gdms31]=gmsdij(Gd31);
[Gds32 Gdm32 Gdms32]=gmsdij(Gd32);

[Gds41 Gdm41 Gdms41]=gmsdij(Gd41);
[Gds42 Gdm42 Gdms42]=gmsdij(Gd42);

[Gds51 Gdm51 Gdms51]=gmsdij(Gd51);
[Gds52 Gdm52 Gdms52]=gmsdij(Gd52);

[Gds61 Gdm61 Gdms61]=gmsdij(Gd61);
[Gds62 Gdm62 Gdms62]=gmsdij(Gd62);

[Gds71 Gdm71 Gdms71]=gmsdij(Gd71);
[Gds72 Gdm72 Gdms72]=gmsdij(Gd72);

[Gds81 Gdm81 Gdms81]=gmsdij(Gd81);
[Gds82 Gdm82 Gdms82]=gmsdij(Gd82);

[Gds91 Gdm91 Gdms91]=gmsdij(Gd91);
[Gds92 Gdm92 Gdms92]=gmsdij(Gd92);

[Gds101 Gdm101 Gdms101]=gmsdij(Gd101);
[Gds102 Gdm102 Gdms102]=gmsdij(Gd102);

[Gds111 Gdm111 Gdms111]=gmsdij(Gd111);
[Gds112 Gdm112 Gdms112]=gmsdij(Gd112);
%% ******** Controllability analysis*************
%*** Sensitivity and complementery sensitivity peak***
%Sensitivity peak for combination of P1 and P2.
sys47 = sys({'P_in','P_t'},'Z2');   
G47 = tf(sys47);
[p, z47]=pzmap(G47);
RHPp47=p(find(p>0));
RHPz47=z47(find(z47>0));

if (isempty(RHPp47) || isempty(RHPz47))
    Msmin47=1;
else
    np47=length(RHPp47); nz47=length(RHPz47);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp47 = zeros(l,np47);
    Yz47 = zeros(l,nz47);
    for i=1:np47
        Yp47(:,i) = Cs([4 7],:)*V(:,i)/norm(Cs([4 7],:)*V(:,i)); %Pole direction
    end
    for j=1:nz47
        [U,S,V]=svd(evalfr(G47,RHPz47(j))); Yz47(:,j)=U(:,end); %zero direction
    end
    Qz47 = zeros(nz47,nz47); Qp47 = zeros(np47,np47); Qzp47 = zeros(nz47,np47); 
    for i=1:nz47
        for j=1:nz47
            Qz47(i,j) = ctranspose(Yz47(:,i))*Yz47(:,j)/(RHPz47(i)+conj(RHPz47(j)));
        end
    end
    for i=1:np47
        for j=1:np47
            Qp47(i,j) = ctranspose(Yp47(:,i))*Yp47(:,j)/(conj(RHPp47(i))+RHPp47(j));
        end
    end
    for i=1:nz47
        for j=1:np47
            Qzp47(i,j) = ctranspose(Yz47(:,i))*Yp47(:,j)/(RHPz47(i)-RHPp47(j));
        end
    end
%     Qp1=(Yp'*Yp).*(1./(diag(RHPp1')*ones(np1)+ones(np1)*diag(RHPp1)));
%     Qz1=(Yz'*Yz).*(1./(diag(RHPz1)*ones(nz1)+ones(nz1)*diag(RHPz1')));
%     Qzp1=(Yz'*Yp).*(1./(diag(RHPz1)*ones(nz1,np1)-ones(nz1,np1)*diag(RHPp1)));
    Msmin47=sqrt(1+norm(sqrtm(inv(Qz47))*Qzp47*sqrtm(inv(Qp47)))^2);
end

%Sensitivity peak for combination of P_t and W_out
sys79 = sys({'P_t','W_out'},'Z2');
G79 = tf(sys79);
[p, z79]=pzmap(G79);
RHPp79=p(find(p>0));
RHPz79=z79(find(z79>0));

if (isempty(RHPp79) || isempty(RHPz79))
    Msmin79=1;
else
    np79=length(RHPp79); nz79=length(RHPz79);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp79 = zeros(l,np79);
    Yz79 = zeros(l,nz79);
    for i=1:np79
        Yp79(:,i) = Cs([7 9],:)*V(:,i)/norm(Cs([7 9],:)*V(:,i)); %Pole direction
    end
    for j=1:nz79
        [U,S,V]=svd(evalfr(G79,RHPz79(j))); Yz79(:,j)=U(:,end); %zero direction
    end
    Qz79 = zeros(nz79,nz79); Qp79 = zeros(np79,np79); Qzp79 = zeros(nz79,np79); 
    for i=1:nz79
        for j=1:nz79
            Qz79(i,j) = ctranspose(Yz79(:,i))*Yz79(:,j)/(RHPz79(i)+conj(RHPz79(j)));
        end
    end
    for i=1:np79
        for j=1:np79
            Qp79(i,j) = ctranspose(Yp79(:,i))*Yp79(:,j)/(conj(RHPp79(i))+RHPp79(j));
        end
    end
    for i=1:nz79
        for j=1:np79
            Qzp79(i,j) = ctranspose(Yz79(:,i))*Yp79(:,j)/(RHPz79(i)-RHPp79(j));
        end
    end
%     Qp1=(Yp'*Yp).*(1./(diag(RHPp1')*ones(np1)+ones(np1)*diag(RHPp1)));
%     Qz1=(Yz'*Yz).*(1./(diag(RHPz1)*ones(nz1)+ones(nz1)*diag(RHPz1')));
%     Qzp1=(Yz'*Yp).*(1./(diag(RHPz1)*ones(nz1,np1)-ones(nz1,np1)*diag(RHPp1)));
    Msmin79=sqrt(1+norm(sqrtm(inv(Qz79))*Qzp79*sqrtm(inv(Qp79)))^2);
end




%Sensitivity peak for combination of P_bh and W_out
sys19 = sys({'P_bh','W_out'},'Z2');
G19 = tf(sys19);
[p, z19]=pzmap(G19);
RHPp19=p(find(p>0));
RHPz19=z19(find(z19>0));

if (isempty(RHPp19) || isempty(RHPz19))
    Msmin19=1;
else
    np19=length(RHPp19); nz19=length(RHPz19);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp19 = zeros(l,np19);
    Yz19 = zeros(l,nz19);
    for i=1:np19
        Yp19(:,i) = Cs([1 9],:)*V(:,i)/norm(Cs([1 9],:)*V(:,i)); %Pole direction
    end
    for j=1:nz19
        [U,S,V]=svd(evalfr(G19,RHPz19(j))); Yz19(:,j)=U(:,end); %zero direction
    end
    Qz19 = zeros(nz19,nz19); Qp19 = zeros(np19,np19); Qzp19 = zeros(nz19,np19); 
    for i=1:nz19
        for j=1:nz19
            Qz19(i,j) = ctranspose(Yz19(:,i))*Yz19(:,j)/(RHPz19(i)+conj(RHPz19(j)));
        end
    end
    for i=1:np19
        for j=1:np19
            Qp19(i,j) = ctranspose(Yp19(:,i))*Yp19(:,j)/(conj(RHPp19(i))+RHPp19(j));
        end
    end
    for i=1:nz19
        for j=1:np19
            Qzp19(i,j) = ctranspose(Yz19(:,i))*Yp19(:,j)/(RHPz19(i)-RHPp19(j));
        end
    end
%     Qp1=(Yp'*Yp).*(1./(diag(RHPp1')*ones(np1)+ones(np1)*diag(RHPp1)));
%     Qz1=(Yz'*Yz).*(1./(diag(RHPz1)*ones(nz1)+ones(nz1)*diag(RHPz1')));
%     Qzp1=(Yz'*Yp).*(1./(diag(RHPz1)*ones(nz1,np1)-ones(nz1,np1)*diag(RHPp1)));
    Msmin19=sqrt(1+norm(sqrtm(inv(Qz19))*Qzp19*sqrtm(inv(Qp19)))^2);
end




%Sensitivity peak for combination of P_in and W_out
sys49 = sys({'P_in','W_out'},'Z2');
G49 = tf(sys49);
[p, z49]=pzmap(G49);
RHPp49=p(find(p>0));
RHPz49=z49(find(z49>0));

if (isempty(RHPp49) || isempty(RHPz49))
    Msmin49=1;
else
    np49=length(RHPp49); nz49=length(RHPz49);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp49 = zeros(l,np49);
    Yz49 = zeros(l,nz49);
    for i=1:np49
        Yp49(:,i) = Cs([4 9],:)*V(:,i)/norm(Cs([4 9],:)*V(:,i)); %Pole direction
    end
    for j=1:nz49
        [U,S,V]=svd(evalfr(G49,RHPz49(j))); Yz49(:,j)=U(:,end); %zero direction
    end
    Qz49 = zeros(nz49,nz49); Qp49 = zeros(np49,np49); Qzp49 = zeros(nz49,np49); 
    for i=1:nz49
        for j=1:nz49
            Qz49(i,j) = ctranspose(Yz49(:,i))*Yz49(:,j)/(RHPz49(i)+conj(RHPz49(j)));
        end
    end
    for i=1:np49
        for j=1:np49
            Qp49(i,j) = ctranspose(Yp49(:,i))*Yp49(:,j)/(conj(RHPp49(i))+RHPp49(j));
        end
    end
    for i=1:nz49
        for j=1:np49
            Qzp49(i,j) = ctranspose(Yz49(:,i))*Yp49(:,j)/(RHPz49(i)-RHPp49(j));
        end
    end
%     Qp1=(Yp'*Yp).*(1./(diag(RHPp1')*ones(np1)+ones(np1)*diag(RHPp1)));
%     Qz1=(Yz'*Yz).*(1./(diag(RHPz1)*ones(nz1)+ones(nz1)*diag(RHPz1')));
%     Qzp1=(Yz'*Yp).*(1./(diag(RHPz1)*ones(nz1,np1)-ones(nz1,np1)*diag(RHPp1)));
    Msmin49=sqrt(1+norm(sqrtm(inv(Qz49))*Qzp49*sqrtm(inv(Qp49)))^2);
end





%Sensitivity peak for combination of P1 and Rho.
sys710 = sys({'P_t','Rho_t'},'Z2');
G710 = tf(sys710);
[p, z710]=pzmap(G710);
RHPp710=p(find(p>0));
RHPz710=z710(find(z710>0));

if (isempty(RHPp710) || isempty(RHPz710))
    Msmin710=1;
else
    np710=length(RHPp710); nz710=length(RHPz710);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp710 = zeros(l,np710);
    Yz710 = zeros(l,nz710);
    for i=1:np710
        Yp710(:,i) = Cs([7 10],:)*V(:,i)/norm(Cs([7 10],:)*V(:,i)); %Pole direction
    end
    for j=1:nz710
        [U,S,V]=svd(evalfr(G710,RHPz710(j))); Yz710(:,j)=U(:,end); %zero direction
    end
    Qz710 = zeros(nz710,nz710); Qp710 = zeros(np710,np710); Qzp710 = zeros(nz710,np710); 
    for i=1:nz710
        for j=1:nz710
            Qz710(i,j) = ctranspose(Yz710(:,i))*Yz710(:,j)/(RHPz710(i)+conj(RHPz710(j)));
        end
    end
    for i=1:np710
        for j=1:np710
            Qp710(i,j) = ctranspose(Yp710(:,i))*Yp710(:,j)/(conj(RHPp710(i))+RHPp710(j));
        end
    end
    for i=1:nz710
        for j=1:np710
            Qzp710(i,j) = ctranspose(Yz710(:,i))*Yp710(:,j)/(RHPz710(i)-RHPp710(j));
        end
    end
%     Qp1=(Yp'*Yp).*(1./(diag(RHPp1')*ones(np1)+ones(np1)*diag(RHPp1)));
%     Qz1=(Yz'*Yz).*(1./(diag(RHPz1)*ones(nz1)+ones(nz1)*diag(RHPz1')));
%     Qzp1=(Yz'*Yp).*(1./(diag(RHPz1)*ones(nz1,np1)-ones(nz1,np1)*diag(RHPp1)));
    Msmin710=sqrt(1+norm(sqrtm(inv(Qz710))*Qzp710*sqrtm(inv(Qp710)))^2);
end

%% *** Bounds on KS ****
%To find the bound for KS the mirror image of the antistable part of the
%G`s are used. If a zero is very close to be negativ, this zero will be
%implemented in the antistable part.

%Lower bound on KS for the output P1 & P2
%Split G into stable and anti stable part, by using a matlab
%function. 
[GSu47, GNSu47]=gnsi(tf(sys47));

tf_mirror=(GNSu47)';
h47=hsvd(tf_mirror);
KSmin47=1/min(h47);

%Lower bound on KS for the output P_t and W_out
%Split G into stable and anti stable part, by using a matlab
%function. 
[GSu79, GNSu79]=gnsi(tf(sys79));

tf_mirror=(GNSu79)';
h79=hsvd(tf_mirror);
KSmin79=1/min(h79);

%Lower bound on KS for the output P_bh and W_out
%Split G into stable and anti stable part, by using a matlab
%function. 
[GSu19, GNSu19]=gnsi(G19);

tf_mirror=(GNSu19)';
h19=hsvd(tf_mirror);
KSmin19=1/min(h19);

%Lower bound on KS for the output P_in and W_out
%Split G into stable and anti stable part, by using a matlab
%function. 
[GSu49, GNSu49]=gnsi(tf(sys49));

tf_mirror=(GNSu49)';
h49=hsvd(tf_mirror);
KSmin49=1/min(h49);

%Lower bound on KS for the output P1 & Rho
%Split G into stable and anti stable part, by using a matlab
%function. 
[GSu710, GNSu710]=gnsi(tf(sys710));

tf_mirror=(GNSu710)';
h710=hsvd(tf_mirror);
KSmin710=1/min(h710);

%% *** Bounds on KSGd **** form 6.26 page 780
%To find the bound for KSGd the mirror image of the antistable part of the
%1/(Gd,ms)*G`s are used. If a zero is very close to be negativ, this zero will be
%implemented in the antistable part.
s=tf('s');
%% Lower bound on KSGd1 (first disturbance) for the output P1, P2 and d1
Gdms_inv = [1/Gdms41 1/Gdms71];
temp = minreal(Gdms_inv*G47);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu471 = gnsdij(temp);
[GSu471, GNSu471]=gnsi(temp);

tf_mirror=(GNSu471)';
h471=hsvd(tf_mirror);
KSGdmin471=1/min(h471);


%% Lower bound on KSGd1 (first disturbance) for the output P1, P2 and d2
Gdms_inv = [1/Gdms42 1/Gdms72];
temp = minreal(Gdms_inv*G47);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu472 = gnsdij(temp);
[GSu472, GNSu472]=gnsi(temp);

tf_mirror=(GNSu472)';
h472=hsvd(tf_mirror);
KSGdmin472=1/min(h472);

%% Lower bound on KS for the outputs P_t and W_out and disturbance d1
%Split G into stable and anti stable part, by using a matlab
%function. 
Gdms_inv = [1/Gdms71  1/Gdms91];
temp = minreal(Gdms_inv*G79);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu791 = gnsdij(temp);
[GSu791, GNSu791]=gnsi(temp);

tf_mirror=(GNSu791)';
h791=hsvd(tf_mirror);
KSGdmin791=1/min(h791);
%% Lower bound on KS for the outputs P_t and W_out and disturbance d2
%Split G into stable and anti stable part, by using a matlab
%function. 
Gdms_inv = [1/Gdms72  1/Gdms92];
temp = minreal(Gdms_inv*G79);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu792 = gnsdij(temp);
[GSu792, GNSu792]=gnsi(temp);

tf_mirror=(GNSu792)';
h792=hsvd(tf_mirror);
KSGdmin792=1/min(h792)

%% Lower bound on KSGd1 for the outputs P_bh and W_out and disturbance d1
%Split G into stable and anti stable part, by using a matlab
%function. 
Gdms_inv = [1/Gdms11  1/Gdms91];
temp = minreal(Gdms_inv*G19);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu191 = gnsdij(temp);
[GSu191, GNSu191]=gnsi(temp);

tf_mirror=(GNSu191)';
h191=hsvd(tf_mirror);
KSGdmin191=1/min(h191);

%% Lower bound on KSGd2 for the outputs P_bh and W_out and disturbance d2
%Split G into stable and anti stable part, by using a matlab
%function. 
Gdms_inv = [1/Gdms12  1/Gdms92];
temp = minreal(Gdms_inv*G19);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu192 = gnsdij(temp);
[GSu192, GNSu192]=gnsi(temp);

tf_mirror=(GNSu192)';
h192=hsvd(tf_mirror);
KSGdmin192=1/min(h192);


%% Lower bound on KSGd1 for the outputs P_in and W_out and disturbance d1
%Split G into stable and anti stable part, by using a matlab
%function. 
Gdms_inv = [1/Gdms41  1/Gdms91];
temp = minreal(Gdms_inv*G49);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu491 = gnsdij(temp);
[GSu491, GNSu491]=gnsi(temp);

tf_mirror=(GNSu491)';
h491=hsvd(tf_mirror);
KSGdmin491=1/min(h491);

%% Lower bound on KSGd2 for the outputs P_in and W_out and disturbance d2
%Split G into stable and anti stable part, by using a matlab
%function. 
Gdms_inv = [1/Gdms42  1/Gdms92];
temp = minreal(Gdms_inv*G49);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu492 = gnsdij(temp);
[GSu492, GNSu492]=gnsi(temp);

tf_mirror=(GNSu492)';
h492=hsvd(tf_mirror);
KSGdmin492=1/min(h492);
%% Lower bound on KS for the outputs P1 & Rho and disturbance d1
%Split G into stable and anti stable part, by using a matlab
%function. 
Gdms_inv = [1/Gdms71  1/Gdms101];
temp = minreal(Gdms_inv*G710);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu7101 = gnsdij(temp);
[GSu7101, GNSu7101]=gnsi(temp);

tf_mirror=(GNSu7101)';
h7101=hsvd(tf_mirror);
KSGdmin7101=1/min(h7101);
%% Lower bound on KS for the outputs P1 & Rho and disturbance d2
%Split G into stable and anti stable part, by using a matlab
%function. 
Gdms_inv = [1/Gdms72  1/Gdms102];
temp = minreal(Gdms_inv*G710);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu7102 = gnsdij(temp);
[GSu7102, GNSu7102]=gnsi(temp);

tf_mirror=(GNSu7102)';
h7102=hsvd(tf_mirror);
KSGdmin7102=1/min(h7102);
%% *** Bounds on KSGd **** form 6.27 page 780
% having only one input, input pole direction is always 1

%% Lower bound on KSGd471 (first disturbance) for the output P1, P2 and d1
Gs47_inv = [1/Gs4 1/Gs7];
Gdms471 = [Gdms41;Gdms71];
RHPp47=p(find(p>0));
np47 = length(RHPp47);
KSGd_min = zeros(1,np47);
for i = np47
    KSGd_min(i) = norm(evalfr(Gs47_inv,RHPp47(i))*evalfr(Gdms471,RHPp47(i)));
end
KSGd471_min = max(KSGd_min);

%% Lower bound on KSGd472 (first disturbance) for the output P1, P2 and d2
Gs47_inv = [1/Gs4 1/Gs7];
Gdms472 = [Gdms42;Gdms72];
RHPp47=p(find(p>0));
np47 = length(RHPp47);
KSGd_min = zeros(1,np47);
for i = np47
    KSGd_min(i) = norm(evalfr(Gs47_inv,RHPp47(i))*evalfr(Gdms472,RHPp47(i)));
end
KSGd472_min = max(KSGd_min);

%% Lower bound on KSGd791 (first disturbance) for the output P_t and W_out and d1
Gs79_inv = [1/Gs7 1/Gs9];
Gdms791 = [Gdms71;Gdms91];
RHPp79=p(find(p>0));
np79 = length(RHPp79);
KSGd_min = zeros(1,np79);
for i = np79
    KSGd_min(i) = norm(evalfr(Gs79_inv,RHPp79(i))*evalfr(Gdms791,RHPp79(i)));
end
KSGd791_min = max(KSGd_min);
%% Lower bound on KSGd791 (first disturbance) for the output P_t and W_out and d2
Gs79_inv = [1/Gs7 1/Gs9];
Gdms792 = [Gdms72;Gdms92];
RHPp79=p(find(p>0));
np79 = length(RHPp79);
KSGd_min = zeros(1,np79);
for i = np79
    KSGd_min(i) = norm(evalfr(Gs79_inv,RHPp79(i))*evalfr(Gdms792,RHPp79(i)));
end
KSGd792_min = max(KSGd_min);

%% Lower bound on KSGd191 (first disturbance) for the output P_bh and W_out and d1
Gs19_inv = [1/Gs1 1/Gs9];
Gdms191 = [Gdms11;Gdms91];
RHPp19=p(find(p>0));
np19 = length(RHPp19);
KSGd_min = zeros(1,np19);
for i = np19
    KSGd_min(i) = norm(evalfr(Gs19_inv,RHPp19(i))*evalfr(Gdms191,RHPp19(i)));
end
KSGd191_min = max(KSGd_min);
%% Lower bound on KSGd191 (second disturbance) for the output P_bh and W_out and d2
Gs19_inv = [1/Gs1 1/Gs9];
Gdms192 = [Gdms12;Gdms92];
RHPp19=p(find(p>0));
np19 = length(RHPp19);
KSGd_min = zeros(1,np19);
for i = np19
    KSGd_min(i) = norm(evalfr(Gs19_inv,RHPp19(i))*evalfr(Gdms192,RHPp19(i)));
end
KSGd192_min = max(KSGd_min);

%% Lower bound on KSGd491 (first disturbance) for the output P_in and W_out and d1
Gs49_inv = [1/Gs4 1/Gs9];
Gdms491 = [Gdms41;Gdms91];
RHPp49=p(find(p>0));
np49 = length(RHPp49);
KSGd_min = zeros(1,np49);
for i = np49
    KSGd_min(i) = norm(evalfr(Gs49_inv,RHPp49(i))*evalfr(Gdms491,RHPp49(i)));
end
KSGd491_min = max(KSGd_min);

%% Lower bound on KSGd491 (second disturbance) for the output P_in and W_out and d2
Gs49_inv = [1/Gs4 1/Gs9];
Gdms492 = [Gdms42;Gdms92];
RHPp49=p(find(p>0));
np49 = length(RHPp49);
KSGd_min = zeros(1,np49);
for i = np49
    KSGd_min(i) = norm(evalfr(Gs49_inv,RHPp49(i))*evalfr(Gdms492,RHPp49(i)));
end
KSGd492_min = max(KSGd_min);

%% Lower bound on KSGd7101 (first disturbance) for the output P1, Rho and d1
Gs710_inv = [1/Gs7 1/Gs10];
Gdms7101 = [Gdms71;Gdms101];
RHPp710=p(find(p>0));
np710 = length(RHPp710);
KSGd_min = zeros(1,np710);
for i = np710
    KSGd_min(i) = norm(evalfr(Gs710_inv,RHPp710(i))*evalfr(Gdms7101,RHPp710(i)));
end
KSGd7101_min = max(KSGd_min);
%% Lower bound on KSGd7101 (first disturbance) for the output P1, Rho and d2
Gs710_inv = [1/Gs7 1/Gs10];
Gdms7102 = [Gdms72;Gdms102];
RHPp710=p(find(p>0));
np710 = length(RHPp710);
KSGd_min = zeros(1,np710);
for i = np710
    KSGd_min(i) = norm(evalfr(Gs710_inv,RHPp710(i))*evalfr(Gdms7102,RHPp710(i)));
end
KSGd7102_min = max(KSGd_min);
%% Bounds on SG47 (P1 and P2).
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = G21ms(s)
sys47 = sys({'P_in','P_t'},'Z2');
z47 = zero(sys47);
l=2; % Number of outputs
RHPp47=p(find(p>0));
RHPz47=z47(find(z47>0));

if (isempty(RHPp47) || isempty(RHPz47))
     GammaSG47_min=0;
else
    np47=length(RHPp47); nz47=length(RHPz47);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp47 = zeros(l,np47);
    Yz47 = zeros(l,nz47);
    for i=1:np47
        Yp47(:,i) = Cs([4 7],:)*V(:,i)/norm(Cs([4 7],:)*V(:,i)); %Pole direction
    end
    for j=1:nz47
        [U,S,V]=svd(evalfr(sys47,RHPz47(j)));
        Yz47(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    W2 = [Gms4; Gms7];
    Qz1 = zeros(nz47,nz47); Qp1 = zeros(np47,np47); Qzp1 = zeros(nz47,np47); 
    Qz2 = zeros(nz47,nz47); Qp2 = zeros(np47,np47); Qzp2 = zeros(nz47,np47); 
    for i=1:nz47
        for j=1:nz47
            W1_zi = W1; %evalfr(W1,RHPz1(i));
            W1_zj = W1; %evalfr(W1,RHPz1(j));
            Qz1(i,j) = ctranspose(Yz47(:,i))*(1/W1_zi)*ctranspose(inv(W1_zj))*Yz47(:,j)/(RHPz47(i)+conj(RHPz47(j)));
            
            W2_zi = evalfr(W2,RHPz47(i));
            W2_zj = evalfr(W2,RHPz47(j));
            Qz2(i,j) = ctranspose(Yz47(:,i))*W2_zi*ctranspose(W2_zj)*Yz47(:,j)/(RHPz47(i)+conj(RHPz47(j)));
        end
    end
    for i=1:np47
        for j=1:np47
            W1_pi = W1; %evalfr(W1,RHPp1(i));
            W1_pj = W1; %evalfr(W1,RHPp1(j));
            Qp1(i,j) = ctranspose(Yp47(:,i))*ctranspose(W1)*W1_pj*Yp47(:,j)/(conj(RHPp47(i))+RHPp47(j));
            
            W2_pi = evalfr(W2,RHPp47(i));
            W2_pj = evalfr(W2,RHPp47(j));
            Qp2(i,j) = ctranspose(Yp47(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp47(:,j)/(conj(RHPp47(i))+RHPp47(j));
        end
    end
    for i=1:nz47
        for j=1:np47
            W1_zi = W1; %evalfr(W1,RHPz1(i));
            W1_pj = W1; %evalfr(W1,RHPp1(j));
            Qzp1(i,j) = ctranspose(Yz47(:,i))*inv(W1_zi)*W1_pj*Yp47(:,j)/(RHPz47(i)-RHPp47(j));
            
            W2_zi = evalfr(W2,RHPz47(i));
            W2_pj = evalfr(W2,RHPp47(j));
            Qzp2(i,j) = ctranspose(Yz47(:,i))*W2_zi*inv(W2_pj)*Yp47(:,j)/(RHPz47(i)-RHPp47(j));
        end
    end
    temp47 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSG47_min = sqrt(max(abs(eig(temp47))));

   
end

%% Bounds on SG79 (P_t and W_out)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = G79ms(s)
sys79 = sys({'P_t','W_out'},'Z2');
l=2; % Number of outputs
RHPp79=p(find(p>0));
z79 = zero(sys79);
RHPz79=z79(find(z79>0));

if (isempty(RHPp79) || isempty(RHPz79))
    GammaSG79_min=0;
else
    np79=length(RHPp79); nz79=length(RHPz79);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp79 = zeros(l,np79);
    Yz79 = zeros(l,nz79);
    for i=1:np79
        Yp79(:,i) = Cs([7 9],:)*V(:,i)/norm(Cs([7 9],:)*V(:,i)); %Pole direction
    end
    for j=1:nz79
        [U,S,V]=svd(evalfr(sys79,RHPz79(j))); 
        Yz79(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gms7;Gms9];

    Qz1 = zeros(nz79,nz79); Qp1 = zeros(np79,np79); Qzp1 = zeros(nz79,np79); 
    Qz2 = zeros(nz79,nz79); Qp2 = zeros(np79,np79); Qzp2 = zeros(nz79,np79); 
    for i=1:nz79
        for j=1:nz79
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz79(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz79(:,j)/(RHPz79(i)+conj(RHPz79(j)));
            
            W2_zi = evalfr(W2,RHPz79(i));
            W2_zj = evalfr(W2,RHPz79(j));
            Qz2(i,j) = ctranspose(Yz79(:,i))*W2_zi*ctranspose(W2_zj)*Yz79(:,j)/(RHPz79(i)+conj(RHPz79(j)));
        end
    end
    for i=1:np79
        for j=1:np79
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp79(:,i))*ctranspose(W1)*W1_pj*Yp79(:,j)/(conj(RHPp79(i))+RHPp79(j));
            
            W2_pi = evalfr(W2,RHPp79(i));
            W2_pj = evalfr(W2,RHPp79(j));
            Qp2(i,j) = ctranspose(Yp79(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp79(:,j)/(conj(RHPp79(i))+RHPp79(j));
        end
    end
    for i=1:nz79
        for j=1:np79
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz79(:,i))*inv(W1_zi)*W1_pj*Yp79(:,j)/(RHPz79(i)-RHPp79(j));
            
            W2_zi = evalfr(W2,RHPz79(i));
            W2_pj = evalfr(W2,RHPp79(j));
            Qzp2(i,j) = ctranspose(Yz79(:,i))*W2_zi*inv(W2_pj)*Yp79(:,j)/(RHPz79(i)-RHPp79(j));
        end
    end
    temp79 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSG79_min = sqrt(max(eig(temp79)));
end

%% Bounds on SG19 (P_bh and W_out)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = G19ms(s)
sys19 = sys({'P_bh','W_out'},'Z2');
l=2; % Number of outputs
RHPp19=p(find(p>0));
z19 = zero(sys19);
RHPz19=z19(find(z19>0));

if (isempty(RHPp19) || isempty(RHPz19))
    GammaSG19_min=0;
else
    np19=length(RHPp19); nz19=length(RHPz19);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp19 = zeros(l,np19);
    Yz19 = zeros(l,nz19);
    for i=1:np19
        Yp19(:,i) = Cs([1 9],:)*V(:,i)/norm(Cs([1 9],:)*V(:,i)); %Pole direction
    end
    for j=1:nz19
        [U,S,V]=svd(evalfr(sys19,RHPz19(j))); 
        Yz19(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gms1;Gms9];

    Qz1 = zeros(nz19,nz19); Qp1 = zeros(np19,np19); Qzp1 = zeros(nz19,np19); 
    Qz2 = zeros(nz19,nz19); Qp2 = zeros(np19,np19); Qzp2 = zeros(nz19,np19); 
    for i=1:nz19
        for j=1:nz19
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz19(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz19(:,j)/(RHPz19(i)+conj(RHPz19(j)));
            
            W2_zi = evalfr(W2,RHPz19(i));
            W2_zj = evalfr(W2,RHPz19(j));
            Qz2(i,j) = ctranspose(Yz19(:,i))*W2_zi*ctranspose(W2_zj)*Yz19(:,j)/(RHPz19(i)+conj(RHPz19(j)));
        end
    end
    for i=1:np19
        for j=1:np19
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp19(:,i))*ctranspose(W1)*W1_pj*Yp19(:,j)/(conj(RHPp19(i))+RHPp19(j));
            
            W2_pi = evalfr(W2,RHPp19(i));
            W2_pj = evalfr(W2,RHPp19(j));
            Qp2(i,j) = ctranspose(Yp19(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp19(:,j)/(conj(RHPp19(i))+RHPp19(j));
        end
    end
    for i=1:nz19
        for j=1:np19
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz19(:,i))*inv(W1_zi)*W1_pj*Yp19(:,j)/(RHPz19(i)-RHPp19(j));
            
            W2_zi = evalfr(W2,RHPz19(i));
            W2_pj = evalfr(W2,RHPp19(j));
            Qzp2(i,j) = ctranspose(Yz19(:,i))*W2_zi*inv(W2_pj)*Yp19(:,j)/(RHPz19(i)-RHPp19(j));
        end
    end
    temp19 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSG19_min = sqrt(max(eig(temp19)));
end

%% Bounds on SG49 (P_in and W_out)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = G49ms(s)
sys49 = sys({'P_in','W_out'},'Z2');
l=2; % Number of outputs
RHPp49=p(find(p>0));
z49 = zero(sys49);
RHPz49=z49(find(z49>0));

if (isempty(RHPp49) || isempty(RHPz49))
    GammaSG49_min=0;
else
    np49=length(RHPp49); nz49=length(RHPz49);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp49 = zeros(l,np49);
    Yz49 = zeros(l,nz49);
    for i=1:np49
        Yp49(:,i) = Cs([4 9],:)*V(:,i)/norm(Cs([4 9],:)*V(:,i)); %Pole direction
    end
    for j=1:nz49
        [U,S,V]=svd(evalfr(sys49,RHPz49(j))); 
        Yz49(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gms4;Gms9];

    Qz1 = zeros(nz49,nz49); Qp1 = zeros(np49,np49); Qzp1 = zeros(nz49,np49); 
    Qz2 = zeros(nz49,nz49); Qp2 = zeros(np49,np49); Qzp2 = zeros(nz49,np49); 
    for i=1:nz49
        for j=1:nz49
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz49(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz49(:,j)/(RHPz49(i)+conj(RHPz49(j)));
            
            W2_zi = evalfr(W2,RHPz49(i));
            W2_zj = evalfr(W2,RHPz49(j));
            Qz2(i,j) = ctranspose(Yz49(:,i))*W2_zi*ctranspose(W2_zj)*Yz49(:,j)/(RHPz49(i)+conj(RHPz49(j)));
        end
    end
    for i=1:np49
        for j=1:np49
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp49(:,i))*ctranspose(W1)*W1_pj*Yp49(:,j)/(conj(RHPp49(i))+RHPp49(j));
            
            W2_pi = evalfr(W2,RHPp49(i));
            W2_pj = evalfr(W2,RHPp49(j));
            Qp2(i,j) = ctranspose(Yp49(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp49(:,j)/(conj(RHPp49(i))+RHPp49(j));
        end
    end
    for i=1:nz49
        for j=1:np49
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz49(:,i))*inv(W1_zi)*W1_pj*Yp49(:,j)/(RHPz49(i)-RHPp49(j));
            
            W2_zi = evalfr(W2,RHPz49(i));
            W2_pj = evalfr(W2,RHPp49(j));
            Qzp2(i,j) = ctranspose(Yz49(:,i))*W2_zi*inv(W2_pj)*Yp49(:,j)/(RHPz49(i)-RHPp49(j));
        end
    end
    temp49 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSG49_min = sqrt(max(eig(temp49)));
end



%% Bounds on SG710 (P1 and Rho)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = G710ms(s)
sys710 = sys({'P_t','Rho_t'},'Z2');
l=2; % Number of outputs
RHPp710=p(find(p>0));
z710 = zero(sys710);
RHPz710=z710(find(z710>0));

if (isempty(RHPp710) || isempty(RHPz710))
    GammaSG710_min=0;
else
    np710=length(RHPp710); nz710=length(RHPz710);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp710 = zeros(l,np710);
    Yz710 = zeros(l,nz710);
    for i=1:np710
        Yp710(:,i) = Cs([7 10],:)*V(:,i)/norm(Cs([7 10],:)*V(:,i)); %Pole direction
    end
    for j=1:nz710
        [U,S,V]=svd(evalfr(sys710,RHPz710(j))); 
        Yz710(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gms7;Gms10];

    Qz1 = zeros(nz710,nz710); Qp1 = zeros(np710,np710); Qzp1 = zeros(nz710,np710); 
    Qz2 = zeros(nz710,nz710); Qp2 = zeros(np710,np710); Qzp2 = zeros(nz710,np710); 
    for i=1:nz710
        for j=1:nz710
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz710(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz710(:,j)/(RHPz710(i)+conj(RHPz710(j)));
            
            W2_zi = evalfr(W2,RHPz710(i));
            W2_zj = evalfr(W2,RHPz710(j));
            Qz2(i,j) = ctranspose(Yz710(:,i))*W2_zi*ctranspose(W2_zj)*Yz710(:,j)/(RHPz710(i)+conj(RHPz710(j)));
        end
    end
    for i=1:np710
        for j=1:np710
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp710(:,i))*ctranspose(W1)*W1_pj*Yp710(:,j)/(conj(RHPp710(i))+RHPp710(j));
            
            W2_pi = evalfr(W2,RHPp710(i));
            W2_pj = evalfr(W2,RHPp710(j));
            Qp2(i,j) = ctranspose(Yp710(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp710(:,j)/(conj(RHPp710(i))+RHPp710(j));
        end
    end
    for i=1:nz710
        for j=1:np710
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz710(:,i))*inv(W1_zi)*W1_pj*Yp710(:,j)/(RHPz710(i)-RHPp710(j));
            
            W2_zi = evalfr(W2,RHPz710(i));
            W2_pj = evalfr(W2,RHPp710(j));
            Qzp2(i,j) = ctranspose(Yz710(:,i))*W2_zi*inv(W2_pj)*Yp710(:,j)/(RHPz710(i)-RHPp710(j));
        end
    end
    temp710 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSG710_min = sqrt(max(eig(temp710)));
end

%% Bounds on SGd471 (P1,P2,d1).
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 =[Gdms11;Gdms21]
sys47 = sys({'P_in','P_t'},'Z2');
z47 = zero(sys47);
l=2; % Number of outputs
RHPp47=p(find(p>0));
RHPz47=z47(find(z47>0));

if (isempty(RHPp47) || isempty(RHPz47))
     GammaSGd471_min=0;
else
    np47=length(RHPp47); nz47=length(RHPz47);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp47 = zeros(l,np47);
    Yz47 = zeros(l,nz47);
    for i=1:np47
        Yp47(:,i) = Cs([4 7],:)*V(:,i)/norm(Cs([4 7],:)*V(:,i)); %Pole direction
    end
    for j=1:nz47
        [U,S,V]=svd(evalfr(sys47,RHPz47(j)));
        Yz47(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    W2 = [Gdms41;Gdms71];
    Qz1 = zeros(nz1,nz47); Qp1 = zeros(np1,np47); Qzp1 = zeros(nz1,np47); 
    Qz2 = zeros(nz1,nz47); Qp2 = zeros(np1,np47); Qzp2 = zeros(nz1,np47); 
    for i=1:nz47
        for j=1:nz47
            W1_zi = W1; %evalfr(W1,RHPz1(i));
            W1_zj = W1; %evalfr(W1,RHPz1(j));
            Qz1(i,j) = ctranspose(Yz47(:,i))*(1/W1_zi)*ctranspose(inv(W1_zj))*Yz47(:,j)/(RHPz47(i)+conj(RHPz47(j)));
            
            W2_zi = evalfr(W2,RHPz47(i));
            W2_zj = evalfr(W2,RHPz47(j));
            Qz2(i,j) = ctranspose(Yz47(:,i))*W2_zi*ctranspose(W2_zj)*Yz47(:,j)/(RHPz47(i)+conj(RHPz47(j)));
        end
    end
    for i=1:np47
        for j=1:np47
            W1_pi = W1; %evalfr(W1,RHPp1(i));
            W1_pj = W1; %evalfr(W1,RHPp1(j));
            Qp1(i,j) = ctranspose(Yp47(:,i))*ctranspose(W1)*W1_pj*Yp47(:,j)/(conj(RHPp47(i))+RHPp47(j));
            
            W2_pi = evalfr(W2,RHPp47(i));
            W2_pj = evalfr(W2,RHPp47(j));
            Qp2(i,j) = ctranspose(Yp47(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp47(:,j)/(conj(RHPp47(i))+RHPp47(j));
        end
    end
    for i=1:nz47
        for j=1:np47
            W1_zi = W1; %evalfr(W1,RHPz1(i));
            W1_pj = W1; %evalfr(W1,RHPp1(j));
            Qzp1(i,j) = ctranspose(Yz47(:,i))*inv(W1_zi)*W1_pj*Yp47(:,j)/(RHPz47(i)-RHPp47(j));
            
            W2_zi = evalfr(W2,RHPz47(i));
            W2_pj = evalfr(W2,RHPp47(j));
            Qzp2(i,j) = ctranspose(Yz47(:,i))*W2_zi*inv(W2_pj)*Yp47(:,j)/(RHPz47(i)-RHPp47(j));
        end
    end
    temp47 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd471_min = sqrt(max(abs(eig(temp47))));
end
%% Bounds on SGd472 (P1,P2,d2).
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 =[Gdms47;Gdms22]
sys47 = sys({'P_in','P_t'},'Z2');
z47 = zero(sys47);
l=2; % Number of outputs
RHPp47=p(find(p>0));
RHPz47=z47(find(z47>0));

if (isempty(RHPp47) || isempty(RHPz47))
     GammaSGd472_min=0;
else
    np47=length(RHPp47); nz47=length(RHPz47);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp47 = zeros(l,np47);
    Yz47 = zeros(l,nz47);
    for i=1:np47
        Yp47(:,i) = Cs([4 7],:)*V(:,i)/norm(Cs([4 7],:)*V(:,i)); %Pole direction
    end
    for j=1:nz47
        [U,S,V]=svd(evalfr(sys47,RHPz47(j)));
        Yz47(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    W2 = [Gdms42;Gdms72];
    Qz1 = zeros(nz47,nz47); Qp1 = zeros(np47,np47); Qzp1 = zeros(nz47,np47); 
    Qz2 = zeros(nz47,nz47); Qp2 = zeros(np47,np47); Qzp2 = zeros(nz47,np47); 
    for i=1:nz47
        for j=1:nz47
            W1_zi = W1; %evalfr(W1,RHPz1(i));
            W1_zj = W1; %evalfr(W1,RHPz1(j));
            Qz1(i,j) = ctranspose(Yz47(:,i))*(1/W1_zi)*ctranspose(inv(W1_zj))*Yz47(:,j)/(RHPz47(i)+conj(RHPz47(j)));
            
            W2_zi = evalfr(W2,RHPz47(i));
            W2_zj = evalfr(W2,RHPz47(j));
            Qz2(i,j) = ctranspose(Yz47(:,i))*W2_zi*ctranspose(W2_zj)*Yz47(:,j)/(RHPz47(i)+conj(RHPz47(j)));
        end
    end
    for i=1:np47
        for j=1:np47
            W1_pi = W1; %evalfr(W1,RHPp1(i));
            W1_pj = W1; %evalfr(W1,RHPp1(j));
            Qp1(i,j) = ctranspose(Yp47(:,i))*ctranspose(W1)*W1_pj*Yp47(:,j)/(conj(RHPp47(i))+RHPp47(j));
            
            W2_pi = evalfr(W2,RHPp47(i));
            W2_pj = evalfr(W2,RHPp47(j));
            Qp2(i,j) = ctranspose(Yp47(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp47(:,j)/(conj(RHPp47(i))+RHPp47(j));
        end
    end
    for i=1:nz47
        for j=1:np47
            W1_zi = W1; %evalfr(W1,RHPz1(i));
            W1_pj = W1; %evalfr(W1,RHPp1(j));
            Qzp1(i,j) = ctranspose(Yz47(:,i))*inv(W1_zi)*W1_pj*Yp47(:,j)/(RHPz47(i)-RHPp47(j));
            
            W2_zi = evalfr(W2,RHPz47(i));
            W2_pj = evalfr(W2,RHPp47(j));
            Qzp2(i,j) = ctranspose(Yz47(:,i))*W2_zi*inv(W2_pj)*Yp47(:,j)/(RHPz47(i)-RHPp47(j));
        end
    end
    temp47 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd472_min = sqrt(max(abs(eig(temp47))));
end
%% Bounds on SGd791 (P_t and W_out,d1)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms21;Gdms31]
sys79 = sys({'P_t','W_out'},'Z2');
l=2; % Number of outputs
RHPp79=p(find(p>0));
z79 = zero(sys79);
RHPz79=z79(find(z79>0));

if (isempty(RHPp79) || isempty(RHPz79))
    GammaSGd791_min=0;
else
    np79=length(RHPp79); nz79=length(RHPz79);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp79 = zeros(l,np79);
    Yz79 = zeros(l,nz79);
    for i=1:np79
        Yp79(:,i) = Cs([7 9],:)*V(:,i)/norm(Cs([7 9],:)*V(:,i)); %Pole direction
    end
    for j=1:nz79
        [U,S,V]=svd(evalfr(sys79,RHPz79(j))); 
        Yz79(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms71;Gdms91];

    Qz1 = zeros(nz79,nz79); Qp1 = zeros(np79,np79); Qzp1 = zeros(nz79,np79); 
    Qz2 = zeros(nz79,nz79); Qp2 = zeros(np79,np79); Qzp2 = zeros(nz79,np79); 
    for i=1:nz79
        for j=1:nz79
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz79(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz79(:,j)/(RHPz79(i)+conj(RHPz79(j)));
            
            W2_zi = evalfr(W2,RHPz79(i));
            W2_zj = evalfr(W2,RHPz79(j));
            Qz2(i,j) = ctranspose(Yz79(:,i))*W2_zi*ctranspose(W2_zj)*Yz79(:,j)/(RHPz79(i)+conj(RHPz79(j)));
        end
    end
    for i=1:np79
        for j=1:np79
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp79(:,i))*ctranspose(W1)*W1_pj*Yp79(:,j)/(conj(RHPp79(i))+RHPp79(j));
            
            W2_pi = evalfr(W2,RHPp79(i));
            W2_pj = evalfr(W2,RHPp79(j));
            Qp2(i,j) = ctranspose(Yp79(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp79(:,j)/(conj(RHPp79(i))+RHPp79(j));
        end
    end
    for i=1:nz79
        for j=1:np79
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz79(:,i))*inv(W1_zi)*W1_pj*Yp79(:,j)/(RHPz79(i)-RHPp79(j));
            
            W2_zi = evalfr(W2,RHPz79(i));
            W2_pj = evalfr(W2,RHPp79(j));
            Qzp2(i,j) = ctranspose(Yz79(:,i))*W2_zi*inv(W2_pj)*Yp79(:,j)/(RHPz79(i)-RHPp79(j));
        end
    end
    temp79 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd791_min = sqrt(max(eig(temp79)));
end

%% Bounds on SGd791 (P_t and W_out,d2)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms22;Gdms32]
sys79 = sys({'P_t','W_out'},'Z2');
l=2; % Number of outputs
RHPp79=p(find(p>0));
z79 = zero(sys79);
RHPz79=z79(find(z79>0));

if (isempty(RHPp79) || isempty(RHPz79))
    GammaSGd792_min=0;
else
    np79=length(RHPp79); nz79=length(RHPz79);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp79 = zeros(l,np79);
    Yz79 = zeros(l,nz79);
    for i=1:np79
        Yp79(:,i) = Cs([7 9],:)*V(:,i)/norm(Cs([7 9],:)*V(:,i)); %Pole direction
    end
    for j=1:nz79
        [U,S,V]=svd(evalfr(sys79,RHPz79(j))); 
        Yz79(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms72;Gdms92];

    Qz1 = zeros(nz79,nz79); Qp1 = zeros(np79,np79); Qzp1 = zeros(nz79,np79); 
    Qz2 = zeros(nz79,nz79); Qp2 = zeros(np79,np79); Qzp2 = zeros(nz79,np79); 
    for i=1:nz79
        for j=1:nz79
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz79(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz79(:,j)/(RHPz79(i)+conj(RHPz79(j)));
            
            W2_zi = evalfr(W2,RHPz79(i));
            W2_zj = evalfr(W2,RHPz79(j));
            Qz2(i,j) = ctranspose(Yz79(:,i))*W2_zi*ctranspose(W2_zj)*Yz79(:,j)/(RHPz79(i)+conj(RHPz79(j)));
        end
    end
    for i=1:np79
        for j=1:np79
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp79(:,i))*ctranspose(W1)*W1_pj*Yp79(:,j)/(conj(RHPp79(i))+RHPp79(j));
            
            W2_pi = evalfr(W2,RHPp79(i));
            W2_pj = evalfr(W2,RHPp79(j));
            Qp2(i,j) = ctranspose(Yp79(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp79(:,j)/(conj(RHPp79(i))+RHPp79(j));
        end
    end
    for i=1:nz79
        for j=1:np79
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz79(:,i))*inv(W1_zi)*W1_pj*Yp79(:,j)/(RHPz79(i)-RHPp79(j));
            
            W2_zi = evalfr(W2,RHPz79(i));
            W2_pj = evalfr(W2,RHPp79(j));
            Qzp2(i,j) = ctranspose(Yz79(:,i))*W2_zi*inv(W2_pj)*Yp79(:,j)/(RHPz79(i)-RHPp79(j));
        end
    end
    temp79 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd792_min = sqrt(max(eig(temp79)));
end


%% Bounds on SGd191 (P_bh and W_out,d1)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms21;Gdms31]
sys19 = sys({'P_bh','W_out'},'Z2');
l=2; % Number of outputs
RHPp19=p(find(p>0));
z19 = zero(sys19);
RHPz19=z19(find(z19>0));

if (isempty(RHPp19) || isempty(RHPz19))
    GammaSGd191_min=0;
else
    np19=length(RHPp19); nz19=length(RHPz19);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp19 = zeros(l,np19);
    Yz19 = zeros(l,nz19);
    for i=1:np19
        Yp19(:,i) = Cs([1 9],:)*V(:,i)/norm(Cs([1 9],:)*V(:,i)); %Pole direction
    end
    for j=1:nz19
        [U,S,V]=svd(evalfr(sys19,RHPz19(j))); 
        Yz19(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms11;Gdms91];

    Qz1 = zeros(nz19,nz19); Qp1 = zeros(np19,np19); Qzp1 = zeros(nz19,np19); 
    Qz2 = zeros(nz19,nz19); Qp2 = zeros(np19,np19); Qzp2 = zeros(nz19,np19); 
    for i=1:nz19
        for j=1:nz19
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz19(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz19(:,j)/(RHPz19(i)+conj(RHPz19(j)));
            
            W2_zi = evalfr(W2,RHPz19(i));
            W2_zj = evalfr(W2,RHPz19(j));
            Qz2(i,j) = ctranspose(Yz19(:,i))*W2_zi*ctranspose(W2_zj)*Yz19(:,j)/(RHPz19(i)+conj(RHPz19(j)));
        end
    end
    for i=1:np19
        for j=1:np19
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp19(:,i))*ctranspose(W1)*W1_pj*Yp19(:,j)/(conj(RHPp19(i))+RHPp19(j));
            
            W2_pi = evalfr(W2,RHPp19(i));
            W2_pj = evalfr(W2,RHPp19(j));
            Qp2(i,j) = ctranspose(Yp19(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp19(:,j)/(conj(RHPp19(i))+RHPp19(j));
        end
    end
    for i=1:nz19
        for j=1:np19
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz19(:,i))*inv(W1_zi)*W1_pj*Yp19(:,j)/(RHPz19(i)-RHPp19(j));
            
            W2_zi = evalfr(W2,RHPz19(i));
            W2_pj = evalfr(W2,RHPp19(j));
            Qzp2(i,j) = ctranspose(Yz19(:,i))*W2_zi*inv(W2_pj)*Yp19(:,j)/(RHPz19(i)-RHPp19(j));
        end
    end
    temp19 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd191_min = sqrt(max(eig(temp19)));
end



%% Bounds on SGd192 (P_bh and W_out,d2)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms22;Gdms32]
sys19 = sys({'P_bh','W_out'},'Z2');
l=2; % Number of outputs
RHPp19=p(find(p>0));
z19 = zero(sys19);
RHPz19=z19(find(z19>0));

if (isempty(RHPp19) || isempty(RHPz19))
    GammaSGd192_min=0;
else
    np19=length(RHPp19); nz19=length(RHPz19);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp19 = zeros(l,np19);
    Yz19 = zeros(l,nz19);
    for i=1:np19
        Yp19(:,i) = Cs([1 9],:)*V(:,i)/norm(Cs([1 9],:)*V(:,i)); %Pole direction
    end
    for j=1:nz19
        [U,S,V]=svd(evalfr(sys19,RHPz19(j))); 
        Yz19(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms12;Gdms92];

    Qz1 = zeros(nz19,nz19); Qp1 = zeros(np19,np19); Qzp1 = zeros(nz19,np19); 
    Qz2 = zeros(nz19,nz19); Qp2 = zeros(np19,np19); Qzp2 = zeros(nz19,np19); 
    for i=1:nz19
        for j=1:nz19
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz19(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz19(:,j)/(RHPz19(i)+conj(RHPz19(j)));
            
            W2_zi = evalfr(W2,RHPz19(i));
            W2_zj = evalfr(W2,RHPz19(j));
            Qz2(i,j) = ctranspose(Yz19(:,i))*W2_zi*ctranspose(W2_zj)*Yz19(:,j)/(RHPz19(i)+conj(RHPz19(j)));
        end
    end
    for i=1:np19
        for j=1:np19
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp19(:,i))*ctranspose(W1)*W1_pj*Yp19(:,j)/(conj(RHPp19(i))+RHPp19(j));
            
            W2_pi = evalfr(W2,RHPp19(i));
            W2_pj = evalfr(W2,RHPp19(j));
            Qp2(i,j) = ctranspose(Yp19(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp19(:,j)/(conj(RHPp19(i))+RHPp19(j));
        end
    end
    for i=1:nz19
        for j=1:np19
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz19(:,i))*inv(W1_zi)*W1_pj*Yp19(:,j)/(RHPz19(i)-RHPp19(j));
            
            W2_zi = evalfr(W2,RHPz19(i));
            W2_pj = evalfr(W2,RHPp19(j));
            Qzp2(i,j) = ctranspose(Yz19(:,i))*W2_zi*inv(W2_pj)*Yp19(:,j)/(RHPz19(i)-RHPp19(j));
        end
    end
    temp19 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd192_min = sqrt(max(eig(temp19)));
end


%% Bounds on SGd491 (P_in and W_out,d1)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms21;Gdms31]
sys49 = sys({'P_in','W_out'},'Z2');
l=2; % Number of outputs
RHPp49=p(find(p>0));
z49 = zero(sys49);
RHPz49=z49(find(z49>0));

if (isempty(RHPp49) || isempty(RHPz49))
    GammaSGd491_min=0;
else
    np49=length(RHPp49); nz49=length(RHPz49);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp49 = zeros(l,np49);
    Yz49 = zeros(l,nz49);
    for i=1:np49
        Yp49(:,i) = Cs([4 9],:)*V(:,i)/norm(Cs([4 9],:)*V(:,i)); %Pole direction
    end
    for j=1:nz49
        [U,S,V]=svd(evalfr(sys49,RHPz49(j))); 
        Yz49(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms41;Gdms91];

    Qz1 = zeros(nz49,nz49); Qp1 = zeros(np49,np49); Qzp1 = zeros(nz49,np49); 
    Qz2 = zeros(nz49,nz49); Qp2 = zeros(np49,np49); Qzp2 = zeros(nz49,np49); 
    for i=1:nz49
        for j=1:nz49
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz49(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz49(:,j)/(RHPz49(i)+conj(RHPz49(j)));
            
            W2_zi = evalfr(W2,RHPz49(i));
            W2_zj = evalfr(W2,RHPz49(j));
            Qz2(i,j) = ctranspose(Yz49(:,i))*W2_zi*ctranspose(W2_zj)*Yz49(:,j)/(RHPz49(i)+conj(RHPz49(j)));
        end
    end
    for i=1:np49
        for j=1:np49
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp49(:,i))*ctranspose(W1)*W1_pj*Yp49(:,j)/(conj(RHPp49(i))+RHPp49(j));
            
            W2_pi = evalfr(W2,RHPp49(i));
            W2_pj = evalfr(W2,RHPp49(j));
            Qp2(i,j) = ctranspose(Yp49(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp49(:,j)/(conj(RHPp49(i))+RHPp49(j));
        end
    end
    for i=1:nz49
        for j=1:np49
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz49(:,i))*inv(W1_zi)*W1_pj*Yp49(:,j)/(RHPz49(i)-RHPp49(j));
            
            W2_zi = evalfr(W2,RHPz49(i));
            W2_pj = evalfr(W2,RHPp49(j));
            Qzp2(i,j) = ctranspose(Yz49(:,i))*W2_zi*inv(W2_pj)*Yp49(:,j)/(RHPz49(i)-RHPp49(j));
        end
    end
    temp49 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd491_min = sqrt(max(eig(temp49)));
end



%% Bounds on SGd492 (P_in and W_out,d2)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms22;Gdms32]
sys49 = sys({'P_in','W_out'},'Z2');
l=2; % Number of outputs
RHPp49=p(find(p>0));
z49 = zero(sys49);
RHPz49=z49(find(z49>0));

if (isempty(RHPp49) || isempty(RHPz49))
    GammaSGd492_min=0;
else
    np49=length(RHPp49); nz49=length(RHPz49);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp49 = zeros(l,np49);
    Yz49 = zeros(l,nz49);
    for i=1:np49
        Yp49(:,i) = Cs([4 9],:)*V(:,i)/norm(Cs([4 9],:)*V(:,i)); %Pole direction
    end
    for j=1:nz49
        [U,S,V]=svd(evalfr(sys49,RHPz49(j))); 
        Yz49(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms42;Gdms92];

    Qz1 = zeros(nz49,nz49); Qp1 = zeros(np49,np49); Qzp1 = zeros(nz49,np49); 
    Qz2 = zeros(nz49,nz49); Qp2 = zeros(np49,np49); Qzp2 = zeros(nz49,np49); 
    for i=1:nz49
        for j=1:nz49
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz49(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz49(:,j)/(RHPz49(i)+conj(RHPz49(j)));
            
            W2_zi = evalfr(W2,RHPz49(i));
            W2_zj = evalfr(W2,RHPz49(j));
            Qz2(i,j) = ctranspose(Yz49(:,i))*W2_zi*ctranspose(W2_zj)*Yz49(:,j)/(RHPz49(i)+conj(RHPz49(j)));
        end
    end
    for i=1:np49
        for j=1:np49
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp49(:,i))*ctranspose(W1)*W1_pj*Yp49(:,j)/(conj(RHPp49(i))+RHPp49(j));
            
            W2_pi = evalfr(W2,RHPp49(i));
            W2_pj = evalfr(W2,RHPp49(j));
            Qp2(i,j) = ctranspose(Yp49(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp49(:,j)/(conj(RHPp49(i))+RHPp49(j));
        end
    end
    for i=1:nz49
        for j=1:np49
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz49(:,i))*inv(W1_zi)*W1_pj*Yp49(:,j)/(RHPz49(i)-RHPp49(j));
            
            W2_zi = evalfr(W2,RHPz49(i));
            W2_pj = evalfr(W2,RHPp49(j));
            Qzp2(i,j) = ctranspose(Yz49(:,i))*W2_zi*inv(W2_pj)*Yp49(:,j)/(RHPz49(i)-RHPp49(j));
        end
    end
    temp49 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd492_min = sqrt(max(eig(temp49)));
end




%% Bounds on SGd7101 (P1,Rho,d1)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms21;Gdms31]
sys710 = sys({'P_t','Rho_t'},'Z2');
l=2; % Number of outputs
RHPp710=p(find(p>0));
z710 = zero(sys710);
RHPz710=z710(find(z710>0));

if (isempty(RHPp710) || isempty(RHPz710))
    GammaSGd7101_min=0;
else
    np710=length(RHPp710); nz710=length(RHPz710);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp710 = zeros(l,np710);
    Yz710 = zeros(l,nz710);
    for i=1:np710
        Yp710(:,i) = Cs([7 10],:)*V(:,i)/norm(Cs([7 10],:)*V(:,i)); %Pole direction
    end
    for j=1:nz710
        [U,S,V]=svd(evalfr(sys710,RHPz710(j))); 
        Yz710(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms71;Gdms101];

    Qz1 = zeros(nz710,nz710); Qp1 = zeros(np710,np710); Qzp1 = zeros(nz710,np710); 
    Qz2 = zeros(nz710,nz710); Qp2 = zeros(np710,np710); Qzp2 = zeros(nz710,np710); 
    for i=1:nz710
        for j=1:nz710
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz710(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz710(:,j)/(RHPz710(i)+conj(RHPz710(j)));
            
            W2_zi = evalfr(W2,RHPz710(i));
            W2_zj = evalfr(W2,RHPz710(j));
            Qz2(i,j) = ctranspose(Yz710(:,i))*W2_zi*ctranspose(W2_zj)*Yz710(:,j)/(RHPz710(i)+conj(RHPz710(j)));
        end
    end
    for i=1:np710
        for j=1:np710
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp710(:,i))*ctranspose(W1)*W1_pj*Yp710(:,j)/(conj(RHPp710(i))+RHPp710(j));
            
            W2_pi = evalfr(W2,RHPp710(i));
            W2_pj = evalfr(W2,RHPp710(j));
            Qp2(i,j) = ctranspose(Yp710(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp710(:,j)/(conj(RHPp710(i))+RHPp710(j));
        end
    end
    for i=1:nz710
        for j=1:np710
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz710(:,i))*inv(W1_zi)*W1_pj*Yp710(:,j)/(RHPz710(i)-RHPp710(j));
            
            W2_zi = evalfr(W2,RHPz710(i));
            W2_pj = evalfr(W2,RHPp710(j));
            Qzp2(i,j) = ctranspose(Yz710(:,i))*W2_zi*inv(W2_pj)*Yp710(:,j)/(RHPz710(i)-RHPp710(j));
        end
    end
    temp710 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd7101_min = sqrt(max(eig(temp710)));
end
%% Bounds on SGd7101 (P1,Rho,d2)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms22;Gdms32]
sys710 = sys({'P_t','Rho_t'},'Z2');
l=2; % Number of outputs
RHPp710=p(find(p>0));
z710 = zero(sys710);
RHPz710=z710(find(z710>0));

if (isempty(RHPp710) || isempty(RHPz710))
    GammaSGd7102_min=0;
else
    np710=length(RHPp710); nz710=length(RHPz710);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp710 = zeros(l,np710);
    Yz710 = zeros(l,nz710);
    for i=1:np710
        Yp710(:,i) = Cs([7 10],:)*V(:,i)/norm(Cs([7 10],:)*V(:,i)); %Pole direction
    end
    for j=1:nz710
        [U,S,V]=svd(evalfr(sys710,RHPz710(j))); 
        Yz710(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms72;Gdms102];

    Qz1 = zeros(nz710,nz710); Qp1 = zeros(np710,np710); Qzp1 = zeros(nz710,np710); 
    Qz2 = zeros(nz710,nz710); Qp2 = zeros(np710,np710); Qzp2 = zeros(nz710,np710); 
    for i=1:nz710
        for j=1:nz710
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz710(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz710(:,j)/(RHPz710(i)+conj(RHPz710(j)));
            
            W2_zi = evalfr(W2,RHPz710(i));
            W2_zj = evalfr(W2,RHPz710(j));
            Qz2(i,j) = ctranspose(Yz710(:,i))*W2_zi*ctranspose(W2_zj)*Yz710(:,j)/(RHPz710(i)+conj(RHPz710(j)));
        end
    end
    for i=1:np710
        for j=1:np710
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp710(:,i))*ctranspose(W1)*W1_pj*Yp710(:,j)/(conj(RHPp710(i))+RHPp710(j));
            
            W2_pi = evalfr(W2,RHPp710(i));
            W2_pj = evalfr(W2,RHPp710(j));
            Qp2(i,j) = ctranspose(Yp710(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp710(:,j)/(conj(RHPp710(i))+RHPp710(j));
        end
    end
    for i=1:nz710
        for j=1:np710
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz710(:,i))*inv(W1_zi)*W1_pj*Yp710(:,j)/(RHPz710(i)-RHPp710(j));
            
            W2_zi = evalfr(W2,RHPz710(i));
            W2_pj = evalfr(W2,RHPp710(j));
            Qzp2(i,j) = ctranspose(Yz710(:,i))*W2_zi*inv(W2_pj)*Yp710(:,j)/(RHPz710(i)-RHPp710(j));
        end
    end
    temp710 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd7102_min = sqrt(max(eig(temp710)));
end
%%
disp('Bounds on S and T');
disp(['The lowest achivable peak for S47 and T47 is: ', num2str(Msmin47)]);
disp(['The lowest achivable peak for S19 and T19 is: ', num2str(Msmin19)]);
disp(['The lowest achivable peak for S49 and T49 is: ', num2str(Msmin49)]);
disp(['The lowest achivable peak for S79 and T79 is: ', num2str(Msmin79)]);
disp(['The lowest achivable peak for S710 and T710 is: ', num2str(Msmin710)]);
disp('##########################################################');
disp('Bounds on KS');
disp(['The lowest achivable peak for KS47 is: ', num2str(KSmin47)]);
disp(['The lowest achivable peak for KS19 is: ', num2str(KSmin19)]);
disp(['The lowest achivable peak for KS49 is: ', num2str(KSmin49)]);
disp(['The lowest achivable peak for KS79 is: ', num2str(KSmin79)]);
disp(['The lowest achivable peak for KS710 is: ', num2str(KSmin710)]);
disp('##########################################################');
disp('Bounds on KSGd form equation 6.26 (tight)');
disp(['The lowest achivable peak for KSGd471 (P1,P2,d1) is: ', num2str(KSGdmin471)]);
disp(['The lowest achivable peak for KSGd472 (P1,P2,d2) is: ', num2str(KSGdmin472)]);
disp(['The lowest achivable peak for KSGd191 (P1,Q,d1) is: ', num2str(KSGdmin191)]);
disp(['The lowest achivable peak for KSGd192 (P1,Q,d2) is: ', num2str(KSGdmin192)]);
disp(['The lowest achivable peak for KSGd491 (P1,Q,d1) is: ', num2str(KSGdmin491)]);
disp(['The lowest achivable peak for KSGd492 (P1,Q,d2) is: ', num2str(KSGdmin492)]);
disp(['The lowest achivable peak for KSGd791 (P2,Q,d1) is: ', num2str(KSGdmin791)]);
disp(['The lowest achivable peak for KSGd792 (P2,Q,d2) is: ', num2str(KSGdmin792)]);
disp(['The lowest achivable peak for KSGd7101 (P2,Rho,d1) is: ', num2str(KSGdmin7101)]);
disp(['The lowest achivable peak for KSGd7102 (P2,Rho,d2) is: ', num2str(KSGdmin7102)]);

disp('##########################################################');
disp('Bounds on KSGd form equation 6.27 (not tight)');
disp(['The lowest achivable peak for KSGd471 (P1,P2,d1) is: ', num2str(KSGd471_min)]);
disp(['The lowest achivable peak for KSGd472 (P1,P2,d2) is: ', num2str(KSGd472_min)]);
disp(['The lowest achivable peak for KSGd191 (P1,Q,d1) is: ', num2str(KSGd191_min)]);
disp(['The lowest achivable peak for KSGd192 (P1,Q,d2) is: ', num2str(KSGd192_min)]);
disp(['The lowest achivable peak for KSGd491 (P1,Q,d1) is: ', num2str(KSGd491_min)]);
disp(['The lowest achivable peak for KSGd492 (P1,Q,d2) is: ', num2str(KSGd492_min)]);
disp(['The lowest achivable peak for KSGd791 (P2,Q,d1) is: ', num2str(KSGd791_min)]);
disp(['The lowest achivable peak for KSGd792 (P2,Q,d2) is: ', num2str(KSGd792_min)]);
disp(['The lowest achivable peak for KSGd7101 (P2,Rho,d1) is: ', num2str(KSGd7101_min)]);
disp(['The lowest achivable peak for KSGd7102 (P2,Rho,d2) is: ', num2str(KSGd7102_min)]);
disp('##########################################################');
disp('Bounds on SG');
disp(['The lowest achivable peak for SG47 (P1,P2) is: ', num2str(GammaSG47_min)]);
disp(['The lowest achivable peak for SG19 (P1,Q) is: ', num2str(GammaSG19_min)]);
disp(['The lowest achivable peak for SG49 (P1,Q) is: ', num2str(GammaSG49_min)]);
disp(['The lowest achivable peak for SG79 (P1,Q) is: ', num2str(GammaSG79_min)]);
disp(['The lowest achivable peak for SG710 (P2,Rho) is: ', num2str(GammaSG710_min)]);
disp('##########################################################');
disp('Bounds on SGd');
disp(['The lowest achivable peak for SGd471 (P1,P2,d1) is: ',num2str(GammaSGd471_min )]);
disp(['The lowest achivable peak for SGd472 (P1,P2,d2) is: ',num2str(GammaSGd472_min )]);
disp(['The lowest achivable peak for SGd191 (P1,Q,d1) is: ',num2str(GammaSGd191_min )]);
disp(['The lowest achivable peak for SGd192 (P1,Q,d2) is: ',num2str(GammaSGd192_min )]);
disp(['The lowest achivable peak for SGd491 (P1,Q,d1) is: ',num2str(GammaSGd491_min )]);
disp(['The lowest achivable peak for SGd492 (P1,Q,d2) is: ',num2str(GammaSGd492_min )]);
disp(['The lowest achivable peak for SGd791 (P2,Q,d1) is: ',num2str(GammaSGd791_min )]);
disp(['The lowest achivable peak for SGd792 (P2,Q,d2) is: ',num2str(GammaSGd792_min )]);
disp(['The lowest achivable peak for SGd7101 (P2,Rho,d1) is: ',num2str(GammaSGd7101_min )]);
disp(['The lowest achivable peak for SGd7102 (P2,Rho,d2) is: ',num2str(GammaSGd7102_min )]);

