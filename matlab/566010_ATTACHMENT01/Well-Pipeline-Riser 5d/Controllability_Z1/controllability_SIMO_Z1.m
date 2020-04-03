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

cd(Current_Folder);

Ds1 = (De^-1)*D(:,1:2)*Du;
Ds2 = (De^-1)*D(:,3:4)*Dd;
Ds = [Ds1 Ds2];
sys =(ss(A,Bs,Cs,Ds));
set(sys,'inputname',{'Z1','Z2','d1','d2'},'outputname',{'P_bh','P_wh','W_in','P_in','P_rb','DP_r','P_t','Q_out','W_out','Rho_t','Alpha_L'},'statename',{'m_Gw','m_Gp','m_Lp','m_Gr','m_Lr'});
sys = minreal(sys);
Gu = tf(sys({'P_bh','P_wh','W_in','P_in','P_rb','DP_r','P_t','Q_out','W_out','Rho_t','Alpha_L'},'Z1'));
Gd = tf(sys({'P_bh','P_wh','W_in','P_in','P_rb','DP_r','P_t','Q_out','W_out','Rho_t','Alpha_L'},{'d1','d2'}));
%% ****Tranfer function and steady state values ****
G1 =(Gu('P_bh','Z1')); G2 =(Gu('P_wh','Z1')); G3 =(Gu('W_in','Z1')); G4 =(Gu('P_in','Z1'));
G5 =(Gu('P_rb','Z1')); G6 =(Gu('DP_r','Z1')); G7 =(Gu('P_t','Z1')); G8 =(Gu('Q_out','Z1'));
G9 =(Gu('W_out','Z1')); G10 =(Gu('Rho_t','Z1')); G11 =(Gu('Alpha_L','Z1'));

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
%Sensitivity peak for combination of P_bh and P_wh.
l=2; % Number of outputs
sys12 = sys({'P_bh','P_wh'},'Z1');
G12 = tf(sys12);
[p, z12]=pzmap(G12);
RHPp12=p(find(p>0));
RHPz12=z12(find(z12>0));

if (isempty(RHPp12) || isempty(RHPz12))
    Msmin12=1;
else
    np12=length(RHPp12); nz12=length(RHPz12);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp12 = zeros(l,np12);
    Yz12 = zeros(l,nz12);
    for i=1:np12
        Yp12(:,i) = Cs([1 2],:)*V(:,i)/norm(Cs([1 2],:)*V(:,i)); %Pole direction
    end
    for j=1:nz12
        [U,S,V]=svd(evalfr(G12,RHPz12(j))); Yz12(:,j)=U(:,end); %zero direction
    end
    Qz12 = zeros(nz12,nz12); Qp12 = zeros(np12,np12); Qzp12 = zeros(nz12,np12); 
    for i=1:nz12
        for j=1:nz12
            Qz12(i,j) = ctranspose(Yz12(:,i))*Yz12(:,j)/(RHPz12(i)+conj(RHPz12(j)));
        end
    end
    for i=1:np12
        for j=1:np12
            Qp12(i,j) = ctranspose(Yp12(:,i))*Yp12(:,j)/(conj(RHPp12(i))+RHPp12(j));
        end
    end
    for i=1:nz12
        for j=1:np12
            Qzp12(i,j) = ctranspose(Yz12(:,i))*Yp12(:,j)/(RHPz12(i)-RHPp12(j));
        end
    end
%     Qp1=(Yp'*Yp).*(1./(diag(RHPp1')*ones(np1)+ones(np1)*diag(RHPp1)));
%     Qz1=(Yz'*Yz).*(1./(diag(RHPz1)*ones(nz1)+ones(nz1)*diag(RHPz1')));
%     Qzp1=(Yz'*Yp).*(1./(diag(RHPz1)*ones(nz1,np1)-ones(nz1,np1)*diag(RHPp1)));
    Msmin12=sqrt(1+norm(sqrtm(inv(Qz12))*Qzp12*sqrtm(inv(Qp12)))^2);
end

%% Sensitivity peak for combination of P_bh and W_in.
l=2; % Number of outputs
sys13 = sys({'P_bh','W_in'},'Z1');
G13 = tf(sys13);
[p, z13]=pzmap(G13);
RHPp13=p(find(p>0));
RHPz13=z13(find(z13>0));

if (isempty(RHPp13) || isempty(RHPz13))
    Msmin13=1;
else
    np13=length(RHPp13); nz13=length(RHPz13);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp13 = zeros(l,np13);
    Yz13 = zeros(l,nz13);
    for i=1:np13
        Yp13(:,i) = Cs([1 3],:)*V(:,i)/norm(Cs([1 3],:)*V(:,i)); %Pole direction
    end
    for j=1:nz13
        [U,S,V]=svd(evalfr(G13,RHPz13(j))); Yz13(:,j)=U(:,end); %zero direction
    end
    Qz13 = zeros(nz13,nz13); Qp13 = zeros(np13,np13); Qzp13 = zeros(nz13,np13); 
    for i=1:nz13
        for j=1:nz13
            Qz13(i,j) = ctranspose(Yz13(:,i))*Yz13(:,j)/(RHPz13(i)+conj(RHPz13(j)));
        end
    end
    for i=1:np13
        for j=1:np13
            Qp13(i,j) = ctranspose(Yp13(:,i))*Yp13(:,j)/(conj(RHPp13(i))+RHPp13(j));
        end
    end
    for i=1:nz13
        for j=1:np13
            Qzp13(i,j) = ctranspose(Yz13(:,i))*Yp13(:,j)/(RHPz13(i)-RHPp13(j));
        end
    end
%     Qp1=(Yp'*Yp).*(1./(diag(RHPp1')*ones(np1)+ones(np1)*diag(RHPp1)));
%     Qz1=(Yz'*Yz).*(1./(diag(RHPz1)*ones(nz1)+ones(nz1)*diag(RHPz1')));
%     Qzp1=(Yz'*Yp).*(1./(diag(RHPz1)*ones(nz1,np1)-ones(nz1,np1)*diag(RHPp1)));
    Msmin13=sqrt(1+norm(sqrtm(inv(Qz13))*Qzp13*sqrtm(inv(Qp13)))^2);
end

%% Sensitivity peak for combination of P_bh and P_in.
l=2; % Number of outputs
sys14 = sys({'P_bh','P_in'},'Z1');
G14 = tf(sys14);
[p, z14]=pzmap(G14);
RHPp14=p(find(p>0));
RHPz14=z14(find(z14>0));

if (isempty(RHPp14) || isempty(RHPz14))
    Msmin14=1;
else
    np14=length(RHPp14); nz14=length(RHPz14);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp14 = zeros(l,np14);
    Yz14 = zeros(l,nz14);
    for i=1:np14
        Yp14(:,i) = Cs([1 4],:)*V(:,i)/norm(Cs([1 4],:)*V(:,i)); %Pole direction
    end
    for j=1:nz14
        [U,S,V]=svd(evalfr(G14,RHPz14(j))); Yz14(:,j)=U(:,end); %zero direction
    end
    Qz14 = zeros(nz14,nz14); Qp14 = zeros(np14,np14); Qzp14 = zeros(nz14,np14); 
    for i=1:nz14
        for j=1:nz14
            Qz14(i,j) = ctranspose(Yz14(:,i))*Yz14(:,j)/(RHPz14(i)+conj(RHPz14(j)));
        end
    end
    for i=1:np14
        for j=1:np14
            Qp14(i,j) = ctranspose(Yp14(:,i))*Yp14(:,j)/(conj(RHPp14(i))+RHPp14(j));
        end
    end
    for i=1:nz14
        for j=1:np14
            Qzp14(i,j) = ctranspose(Yz14(:,i))*Yp14(:,j)/(RHPz14(i)-RHPp14(j));
        end
    end
%     Qp1=(Yp'*Yp).*(1./(diag(RHPp1')*ones(np1)+ones(np1)*diag(RHPp1)));
%     Qz1=(Yz'*Yz).*(1./(diag(RHPz1)*ones(nz1)+ones(nz1)*diag(RHPz1')));
%     Qzp1=(Yz'*Yp).*(1./(diag(RHPz1)*ones(nz1,np1)-ones(nz1,np1)*diag(RHPp1)));
    Msmin14=sqrt(1+norm(sqrtm(inv(Qz14))*Qzp14*sqrtm(inv(Qp14)))^2);
end

%% Sensitivity peak for combination of P_wh and W_in.
l=2; % Number of outputs
sys23 = sys({'P_wh','W_in'},'Z1');
G23 = tf(sys23);
[p, z23]=pzmap(G23);
RHPp23=p(find(p>0));
RHPz23=z23(find(z23>0));

if (isempty(RHPp23) || isempty(RHPz23))
    Msmin23=1;
else
    np23=length(RHPp23); nz23=length(RHPz23);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp23 = zeros(l,np23);
    Yz23 = zeros(l,nz23);
    for i=1:np23
        Yp23(:,i) = Cs([2 3],:)*V(:,i)/norm(Cs([2 3],:)*V(:,i)); %Pole direction
    end
    for j=1:nz23
        [U,S,V]=svd(evalfr(G23,RHPz23(j))); Yz23(:,j)=U(:,end); %zero direction
    end
    Qz23 = zeros(nz23,nz23); Qp23 = zeros(np23,np23); Qzp23 = zeros(nz23,np23); 
    for i=1:nz23
        for j=1:nz23
            Qz23(i,j) = ctranspose(Yz23(:,i))*Yz23(:,j)/(RHPz23(i)+conj(RHPz23(j)));
        end
    end
    for i=1:np23
        for j=1:np23
            Qp23(i,j) = ctranspose(Yp23(:,i))*Yp23(:,j)/(conj(RHPp23(i))+RHPp23(j));
        end
    end
    for i=1:nz23
        for j=1:np23
            Qzp23(i,j) = ctranspose(Yz23(:,i))*Yp23(:,j)/(RHPz23(i)-RHPp23(j));
        end
    end
%     Qp1=(Yp'*Yp).*(1./(diag(RHPp1')*ones(np1)+ones(np1)*diag(RHPp1)));
%     Qz1=(Yz'*Yz).*(1./(diag(RHPz1)*ones(nz1)+ones(nz1)*diag(RHPz1')));
%     Qzp1=(Yz'*Yp).*(1./(diag(RHPz1)*ones(nz1,np1)-ones(nz1,np1)*diag(RHPp1)));
    Msmin23=sqrt(1+norm(sqrtm(inv(Qz23))*Qzp23*sqrtm(inv(Qp23)))^2);
end

%% Sensitivity peak for combination of P_wh and P_in.
l=2; % Number of outputs
sys24 = sys({'P_wh','P_in'},'Z1');
G24 = tf(sys24);
[p, z24]=pzmap(G24);
RHPp24=p(find(p>0));
RHPz24=z24(find(z24>0));

if (isempty(RHPp24) || isempty(RHPz24))
    Msmin24=1;
else
    np24=length(RHPp24); nz24=length(RHPz24);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp24 = zeros(l,np24);
    Yz24 = zeros(l,nz24);
    for i=1:np24
        Yp24(:,i) = Cs([2 4],:)*V(:,i)/norm(Cs([2 4],:)*V(:,i)); %Pole direction
    end
    for j=1:nz24
        [U,S,V]=svd(evalfr(G24,RHPz24(j))); Yz24(:,j)=U(:,end); %zero direction
    end
    Qz24 = zeros(nz24,nz24); Qp24 = zeros(np24,np24); Qzp24 = zeros(nz24,np24); 
    for i=1:nz24
        for j=1:nz24
            Qz24(i,j) = ctranspose(Yz24(:,i))*Yz24(:,j)/(RHPz24(i)+conj(RHPz24(j)));
        end
    end
    for i=1:np24
        for j=1:np24
            Qp24(i,j) = ctranspose(Yp24(:,i))*Yp24(:,j)/(conj(RHPp24(i))+RHPp24(j));
        end
    end
    for i=1:nz24
        for j=1:np24
            Qzp24(i,j) = ctranspose(Yz24(:,i))*Yp24(:,j)/(RHPz24(i)-RHPp24(j));
        end
    end
%     Qp1=(Yp'*Yp).*(1./(diag(RHPp1')*ones(np1)+ones(np1)*diag(RHPp1)));
%     Qz1=(Yz'*Yz).*(1./(diag(RHPz1)*ones(nz1)+ones(nz1)*diag(RHPz1')));
%     Qzp1=(Yz'*Yp).*(1./(diag(RHPz1)*ones(nz1,np1)-ones(nz1,np1)*diag(RHPp1)));
    Msmin24=sqrt(1+norm(sqrtm(inv(Qz24))*Qzp24*sqrtm(inv(Qp24)))^2);
end
%% Sensitivity peak for combination of W_in and P_in.
l=2; % Number of outputs
sys34 = sys({'W_in','P_in'},'Z1');
G34 = tf(sys34);
[p, z34]=pzmap(G34);
RHPp34=p(find(p>0));
RHPz34=z34(find(z34>0));

if (isempty(RHPp34) || isempty(RHPz34))
    Msmin34=1;
else
    np34=length(RHPp34); nz34=length(RHPz34);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp34 = zeros(l,np34);
    Yz34 = zeros(l,nz34);
    for i=1:np34
        Yp34(:,i) = Cs([3 4],:)*V(:,i)/norm(Cs([3 4],:)*V(:,i)); %Pole direction
    end
    for j=1:nz34
        [U,S,V]=svd(evalfr(G34,RHPz34(j))); Yz34(:,j)=U(:,end); %zero direction
    end
    Qz34 = zeros(nz34,nz34); Qp34 = zeros(np34,np34); Qzp34 = zeros(nz34,np34); 
    for i=1:nz34
        for j=1:nz34
            Qz34(i,j) = ctranspose(Yz34(:,i))*Yz34(:,j)/(RHPz34(i)+conj(RHPz34(j)));
        end
    end
    for i=1:np34
        for j=1:np34
            Qp34(i,j) = ctranspose(Yp34(:,i))*Yp34(:,j)/(conj(RHPp34(i))+RHPp34(j));
        end
    end
    for i=1:nz34
        for j=1:np34
            Qzp34(i,j) = ctranspose(Yz34(:,i))*Yp34(:,j)/(RHPz34(i)-RHPp34(j));
        end
    end
%     Qp1=(Yp'*Yp).*(1./(diag(RHPp1')*ones(np1)+ones(np1)*diag(RHPp1)));
%     Qz1=(Yz'*Yz).*(1./(diag(RHPz1)*ones(nz1)+ones(nz1)*diag(RHPz1')));
%     Qzp1=(Yz'*Yp).*(1./(diag(RHPz1)*ones(nz1,np1)-ones(nz1,np1)*diag(RHPp1)));
    Msmin34=sqrt(1+norm(sqrtm(inv(Qz34))*Qzp34*sqrtm(inv(Qp34)))^2);
end

%% *** Bounds on KS ****
%To find the bound for KS the mirror image of the antistable part of the
%G`s are used. If a zero is very close to be negativ, this zero will be
%implemented in the antistable part.

%Lower bound on KS for the output P_bh & P_wh
%Split G into stable and anti stable part, by using a matlab
%function. 
[GSu12, GNSu12]=gnsi(tf(sys12));

tf_mirror=(GNSu12)';
h12=hsvd(tf_mirror);
KSmin12=1/min(h12);

%Lower bound on KS for the output P_bh & W_in
%Split G into stable and anti stable part, by using a matlab
%function. 
[GSu13, GNSu13]=gnsi(tf(sys13));

tf_mirror=(GNSu13)';
h13=hsvd(tf_mirror);
KSmin13=1/min(h13);

%Lower bound on KS for the output P_bh & P_in
%Split G into stable and anti stable part, by using a matlab
%function. 
[GSu14, GNSu14]=gnsi(tf(sys14));

tf_mirror=(GNSu14)';
h14=hsvd(tf_mirror);
KSmin14=1/min(h14);

%Lower bound on KS for the output P_wh & W_in
%Split G into stable and anti stable part, by using a matlab
%function. 
[GSu23, GNSu23]=gnsi(tf(sys23));

tf_mirror=(GNSu23)';
h23=hsvd(tf_mirror);
KSmin23=1/min(h23);

%Lower bound on KS for the output P_wh & P_in
%Split G into stable and anti stable part, by using a matlab
%function. 
[GSu24, GNSu24]=gnsi(tf(sys24));

tf_mirror=(GNSu24)';
h24=hsvd(tf_mirror);
KSmin24=1/min(h24);

%Lower bound on KS for the output W_in & P_in
%Split G into stable and anti stable part, by using a matlab
%function. 
[GSu34, GNSu34]=gnsi(tf(sys34));

tf_mirror=(GNSu34)';
h34=hsvd(tf_mirror);
KSmin34=1/min(h34);

%% *** Bounds on KSGd **** form 6.26 page 780
%To find the bound for KSGd the mirror image of the antistable part of the
%1/(Gd,ms)*G`s are used. If a zero is very close to be negativ, this zero will be
%implemented in the antistable part.
s=tf('s');
%% Lower bound on KSGd1 (first disturbance) for the output P_bh, P_wh and d1
Gdms_inv = [1/Gdms11 1/Gdms21];
temp = minreal(Gdms_inv*G12);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu471 = gnsdij(temp);
[GSu121, GNSu121]=gnsi(temp);

tf_mirror=(GNSu121)';
h121=hsvd(tf_mirror);
KSGdmin121=1/min(h121);

%% Lower bound on KSGd1 (first disturbance) for the output P_bh, W_in and d1
Gdms_inv = [1/Gdms11 1/Gdms31];
temp = minreal(Gdms_inv*G13);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu471 = gnsdij(temp);
[GSu131, GNSu131]=gnsi(temp);

tf_mirror=(GNSu131)';
h131=hsvd(tf_mirror);
KSGdmin131=1/min(h131);

%% Lower bound on KSGd1 (first disturbance) for the output P_bh, P_in and d1
Gdms_inv = [1/Gdms11 1/Gdms41];
temp = minreal(Gdms_inv*G14);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu471 = gnsdij(temp);
[GSu141, GNSu141]=gnsi(temp);

tf_mirror=(GNSu141)';
h141=hsvd(tf_mirror);
KSGdmin141=1/min(h141);

%% Lower bound on KSGd1 (first disturbance) for the output P_wh, W_in and d1
Gdms_inv = [1/Gdms21 1/Gdms31];
temp = minreal(Gdms_inv*G23);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu471 = gnsdij(temp);
[GSu231, GNSu231]=gnsi(temp);

tf_mirror=(GNSu231)';
h231=hsvd(tf_mirror);
KSGdmin231=1/min(h231);

%% Lower bound on KSGd1 (first disturbance) for the output P_wh, P_in and d1
Gdms_inv = [1/Gdms21 1/Gdms41];
temp = minreal(Gdms_inv*G24);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu471 = gnsdij(temp);
[GSu241, GNSu241]=gnsi(temp);

tf_mirror=(GNSu241)';
h241=hsvd(tf_mirror);
KSGdmin241=1/min(h241);

%% Lower bound on KSGd1 (first disturbance) for the output P_wh, P_in and d1
Gdms_inv = [1/Gdms31 1/Gdms41];
temp = minreal(Gdms_inv*G34);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu471 = gnsdij(temp);
[GSu341, GNSu341]=gnsi(temp);

tf_mirror=(GNSu341)';
h341=hsvd(tf_mirror);
KSGdmin341=1/min(h341);

%% Lower bound on KSGd2 (second disturbance) for the output P_bh, P_wh and d2
Gdms_inv = [1/Gdms12 1/Gdms22];
temp = minreal(Gdms_inv*G12);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu122 = gnsdij(temp);
[GSu122, GNSu122]=gnsi(temp);

tf_mirror=(GNSu122)';
h122=hsvd(tf_mirror);
KSGdmin122=1/min(h122);

%% Lower bound on KSGd2 (second disturbance) for the output P_bh, W_in and d2
Gdms_inv = [1/Gdms12 1/Gdms32];
temp = minreal(Gdms_inv*G13);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu132 = gnsdij(temp);
[GSu132, GNSu132]=gnsi(temp);

tf_mirror=(GNSu132)';
h132=hsvd(tf_mirror);
KSGdmin132=1/min(h132);

%% Lower bound on KSGd2 (second disturbance) for the output P_bh, P_in and d2
Gdms_inv = [1/Gdms12 1/Gdms42];
temp = minreal(Gdms_inv*G14);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu142 = gnsdij(temp);
[GSu142, GNSu142]=gnsi(temp);

tf_mirror=(GNSu142)';
h142=hsvd(tf_mirror);
KSGdmin142=1/min(h142);

%% Lower bound on KSGd2 (second disturbance) for the output P_wh, W_in and d2
Gdms_inv = [1/Gdms22 1/Gdms32];
temp = minreal(Gdms_inv*G23);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu232 = gnsdij(temp);
[GSu232, GNSu232]=gnsi(temp);

tf_mirror=(GNSu232)';
h232=hsvd(tf_mirror);
KSGdmin232=1/min(h232);

%% Lower bound on KSGd2 (second disturbance) for the output P_wh, P_in and d2
Gdms_inv = [1/Gdms22 1/Gdms42];
temp = minreal(Gdms_inv*G24);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu242 = gnsdij(temp);
[GSu242, GNSu242]=gnsi(temp);

tf_mirror=(GNSu242)';
h242=hsvd(tf_mirror);
KSGdmin242=1/min(h242);

%% Lower bound on KSGd2 (second disturbance) for the output W_in, P_in and d2
Gdms_inv = [1/Gdms32 1/Gdms42];
temp = minreal(Gdms_inv*G34);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu342 = gnsdij(temp);
[GSu342, GNSu342]=gnsi(temp);

tf_mirror=(GNSu342)';
h342=hsvd(tf_mirror);
KSGdmin342=1/min(h342);

%% Bounds on SG12 (P_bh, P_wh).
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = G21ms(s)
sys12 = sys({'P_bh','P_wh'},'Z1');
z12 = zero(sys12);
l=2; % Number of outputs
RHPp12=p(find(p>0));
RHPz12=z12(find(z12>0));

if (isempty(RHPp12) || isempty(RHPz12))
     GammaSG12_min=0;
else
    np12=length(RHPp12); nz12=length(RHPz12);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp12 = zeros(l,np12);
    Yz12 = zeros(l,nz12);
    for i=1:np12
        Yp12(:,i) = Cs([1 2],:)*V(:,i)/norm(Cs([1 2],:)*V(:,i)); %Pole direction
    end
    for j=1:nz12
        [U,S,V]=svd(evalfr(sys12,RHPz12(j)));
        Yz12(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    W2 = [Gms1; Gms2];
    Qz1 = zeros(nz12,nz12); Qp1 = zeros(np12,np12); Qzp1 = zeros(nz12,np12); 
    Qz2 = zeros(nz12,nz12); Qp2 = zeros(np12,np12); Qzp2 = zeros(nz12,np12); 
    for i=1:nz12
        for j=1:nz12
            W1_zi = W1; %evalfr(W1,RHPz1(i));
            W1_zj = W1; %evalfr(W1,RHPz1(j));
            Qz1(i,j) = ctranspose(Yz12(:,i))*(1/W1_zi)*ctranspose(inv(W1_zj))*Yz12(:,j)/(RHPz12(i)+conj(RHPz12(j)));
            
            W2_zi = evalfr(W2,RHPz12(i));
            W2_zj = evalfr(W2,RHPz12(j));
            Qz2(i,j) = ctranspose(Yz12(:,i))*W2_zi*ctranspose(W2_zj)*Yz12(:,j)/(RHPz12(i)+conj(RHPz12(j)));
        end
    end
    for i=1:np12
        for j=1:np12
            W1_pi = W1; %evalfr(W1,RHPp1(i));
            W1_pj = W1; %evalfr(W1,RHPp1(j));
            Qp1(i,j) = ctranspose(Yp12(:,i))*ctranspose(W1)*W1_pj*Yp12(:,j)/(conj(RHPp12(i))+RHPp12(j));
            
            W2_pi = evalfr(W2,RHPp12(i));
            W2_pj = evalfr(W2,RHPp12(j));
            Qp2(i,j) = ctranspose(Yp12(:,i))*ctranspose(pinv(W2_pi))*pinv(W2_pj)*Yp12(:,j)/(conj(RHPp12(i))+RHPp12(j));
        end
    end
    for i=1:nz12
        for j=1:np12
            W1_zi = W1; %evalfr(W1,RHPz1(i));
            W1_pj = W1; %evalfr(W1,RHPp1(j));
            Qzp1(i,j) = ctranspose(Yz12(:,i))*inv(W1_zi)*W1_pj*Yp12(:,j)/(RHPz12(i)-RHPp12(j));
            
            W2_zi = evalfr(W2,RHPz12(i));
            W2_pj = evalfr(W2,RHPp12(j));
            Qzp2(i,j) = ctranspose(Yz12(:,i))*W2_zi*pinv(W2_pj)*Yp12(:,j)/(RHPz12(i)-RHPp12(j));
        end
    end
    temp12 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSG12_min = sqrt(max(abs(eig(temp12))));
end

%% Bounds on SG13 (P_bh and W_in)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = G13ms(s)
sys13 = sys({'P_bh','W_in'},'Z1');
l=2; % Number of outputs
RHPp13=p(find(p>0));
z13 = zero(sys13);
RHPz13=z13(find(z13>0));

if (isempty(RHPp13) || isempty(RHPz13))
    GammaSG13_min=0;
else
    np13=length(RHPp13); nz13=length(RHPz13);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp13 = zeros(l,np13);
    Yz13 = zeros(l,nz13);
    for i=1:np13
        Yp13(:,i) = Cs([1 3],:)*V(:,i)/norm(Cs([1 3],:)*V(:,i)); %Pole direction
    end
    for j=1:nz13
        [U,S,V]=svd(evalfr(sys13,RHPz13(j))); 
        Yz13(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gms1;Gms3];

    Qz1 = zeros(nz13,nz13); Qp1 = zeros(np13,np13); Qzp1 = zeros(nz13,np13); 
    Qz2 = zeros(nz13,nz13); Qp2 = zeros(np13,np13); Qzp2 = zeros(nz13,np13); 
    for i=1:nz13
        for j=1:nz13
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz13(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz13(:,j)/(RHPz13(i)+conj(RHPz13(j)));
            
            W2_zi = evalfr(W2,RHPz13(i));
            W2_zj = evalfr(W2,RHPz13(j));
            Qz2(i,j) = ctranspose(Yz13(:,i))*W2_zi*ctranspose(W2_zj)*Yz13(:,j)/(RHPz13(i)+conj(RHPz13(j)));
        end
    end
    for i=1:np13
        for j=1:np13
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp13(:,i))*ctranspose(W1)*W1_pj*Yp13(:,j)/(conj(RHPp13(i))+RHPp13(j));
            
            W2_pi = evalfr(W2,RHPp13(i));
            W2_pj = evalfr(W2,RHPp13(j));
            Qp2(i,j) = ctranspose(Yp13(:,i))*ctranspose(pinv(W2_pi))*pinv(W2_pj)*Yp13(:,j)/(conj(RHPp13(i))+RHPp13(j));
        end
    end
    for i=1:nz13
        for j=1:np13
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz13(:,i))*inv(W1_zi)*W1_pj*Yp13(:,j)/(RHPz13(i)-RHPp13(j));
            
            W2_zi = evalfr(W2,RHPz13(i));
            W2_pj = evalfr(W2,RHPp13(j));
            Qzp2(i,j) = ctranspose(Yz13(:,i))*W2_zi*pinv(W2_pj)*Yp13(:,j)/(RHPz13(i)-RHPp13(j));
        end
    end
    temp13 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSG13_min = real(sqrt(max(eig(temp13))));
end

%% Bounds on SG14 (P_bh and P_in)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = G14ms(s)
sys14 = sys({'P_bh','P_in'},'Z1');
l=2; % Number of outputs
RHPp14=p(find(p>0));
z14 = zero(sys14);
RHPz14=z14(find(z14>0));

if (isempty(RHPp14) || isempty(RHPz14))
    GammaSG14_min=0;
else
    np14=length(RHPp14); nz14=length(RHPz14);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp14 = zeros(l,np14);
    Yz14 = zeros(l,nz14);
    for i=1:np14
        Yp14(:,i) = Cs([1 4],:)*V(:,i)/norm(Cs([1 4],:)*V(:,i)); %Pole direction
    end
    for j=1:nz14
        [U,S,V]=svd(evalfr(sys14,RHPz14(j))); 
        Yz14(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gms1;Gms4];

    Qz1 = zeros(nz14,nz14); Qp1 = zeros(np14,np14); Qzp1 = zeros(nz14,np14); 
    Qz2 = zeros(nz14,nz14); Qp2 = zeros(np14,np14); Qzp2 = zeros(nz14,np14); 
    for i=1:nz14
        for j=1:nz14
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz14(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz14(:,j)/(RHPz14(i)+conj(RHPz14(j)));
            
            W2_zi = evalfr(W2,RHPz14(i));
            W2_zj = evalfr(W2,RHPz14(j));
            Qz2(i,j) = ctranspose(Yz14(:,i))*W2_zi*ctranspose(W2_zj)*Yz14(:,j)/(RHPz14(i)+conj(RHPz14(j)));
        end
    end
    for i=1:np14
        for j=1:np14
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp14(:,i))*ctranspose(W1)*W1_pj*Yp14(:,j)/(conj(RHPp14(i))+RHPp14(j));
            
            W2_pi = evalfr(W2,RHPp14(i));
            W2_pj = evalfr(W2,RHPp14(j));
            Qp2(i,j) = ctranspose(Yp14(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp14(:,j)/(conj(RHPp14(i))+RHPp14(j));
        end
    end
    for i=1:nz14
        for j=1:np14
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz14(:,i))*inv(W1_zi)*W1_pj*Yp14(:,j)/(RHPz14(i)-RHPp14(j));
            
            W2_zi = evalfr(W2,RHPz14(i));
            W2_pj = evalfr(W2,RHPp14(j));
            Qzp2(i,j) = ctranspose(Yz14(:,i))*W2_zi*inv(W2_pj)*Yp14(:,j)/(RHPz14(i)-RHPp14(j));
        end
    end
    temp14 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSG14_min = sqrt(max(eig(temp14)));
end

%% Bounds on SG23 (P_wh and W_in)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = G23ms(s)
sys23 = sys({'P_wh','W_in'},'Z1');
l=2; % Number of outputs
RHPp23=p(find(p>0));
z23 = zero(sys23);
RHPz23=z23(find(z23>0));

if (isempty(RHPp23) || isempty(RHPz23))
    GammaSG23_min=0;
else
    np23=length(RHPp23); nz23=length(RHPz23);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp23 = zeros(l,np23);
    Yz23 = zeros(l,nz23);
    for i=1:np23
        Yp23(:,i) = Cs([2 3],:)*V(:,i)/norm(Cs([2 3],:)*V(:,i)); %Pole direction
    end
    for j=1:nz23
        [U,S,V]=svd(evalfr(sys23,RHPz23(j))); 
        Yz23(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gms2;Gms3];

    Qz1 = zeros(nz23,nz23); Qp1 = zeros(np23,np23); Qzp1 = zeros(nz23,np23); 
    Qz2 = zeros(nz23,nz23); Qp2 = zeros(np23,np23); Qzp2 = zeros(nz23,np23); 
    for i=1:nz23
        for j=1:nz23
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz23(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz23(:,j)/(RHPz23(i)+conj(RHPz23(j)));
            
            W2_zi = evalfr(W2,RHPz23(i));
            W2_zj = evalfr(W2,RHPz23(j));
            Qz2(i,j) = ctranspose(Yz23(:,i))*W2_zi*ctranspose(W2_zj)*Yz23(:,j)/(RHPz23(i)+conj(RHPz23(j)));
        end
    end
    for i=1:np23
        for j=1:np23
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp23(:,i))*ctranspose(W1)*W1_pj*Yp23(:,j)/(conj(RHPp23(i))+RHPp23(j));
            
            W2_pi = evalfr(W2,RHPp23(i));
            W2_pj = evalfr(W2,RHPp23(j));
            Qp2(i,j) = ctranspose(Yp23(:,i))*ctranspose(pinv(W2_pi))*pinv(W2_pj)*Yp23(:,j)/(conj(RHPp23(i))+RHPp23(j));
        end
    end
    for i=1:nz23
        for j=1:np23
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz23(:,i))*inv(W1_zi)*W1_pj*Yp23(:,j)/(RHPz23(i)-RHPp23(j));
            
            W2_zi = evalfr(W2,RHPz23(i));
            W2_pj = evalfr(W2,RHPp23(j));
            Qzp2(i,j) = ctranspose(Yz23(:,i))*W2_zi*pinv(W2_pj)*Yp23(:,j)/(RHPz23(i)-RHPp23(j));
        end
    end
    temp23 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSG23_min = real(sqrt(max(eig(temp23))));
end

%% Bounds on SG24 (P_wh and P_in)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = G24ms(s)
sys24 = sys({'P_wh','P_in'},'Z1');
l=2; % Number of outputs
RHPp24=p(find(p>0));
z24 = zero(sys24);
RHPz24=z24(find(z24>0));

if (isempty(RHPp24) || isempty(RHPz24))
    GammaSG24_min=0;
else
    np24=length(RHPp24); nz24=length(RHPz24);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp24 = zeros(l,np24);
    Yz24 = zeros(l,nz24);
    for i=1:np24
        Yp24(:,i) = Cs([2 4],:)*V(:,i)/norm(Cs([2 4],:)*V(:,i)); %Pole direction
    end
    for j=1:nz24
        [U,S,V]=svd(evalfr(sys24,RHPz24(j))); 
        Yz24(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gms2;Gms4];

    Qz1 = zeros(nz24,nz24); Qp1 = zeros(np24,np24); Qzp1 = zeros(nz24,np24); 
    Qz2 = zeros(nz24,nz24); Qp2 = zeros(np24,np24); Qzp2 = zeros(nz24,np24); 
    for i=1:nz24
        for j=1:nz24
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz24(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz24(:,j)/(RHPz24(i)+conj(RHPz24(j)));
            
            W2_zi = evalfr(W2,RHPz24(i));
            W2_zj = evalfr(W2,RHPz24(j));
            Qz2(i,j) = ctranspose(Yz24(:,i))*W2_zi*ctranspose(W2_zj)*Yz24(:,j)/(RHPz24(i)+conj(RHPz24(j)));
        end
    end
    for i=1:np24
        for j=1:np24
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp24(:,i))*ctranspose(W1)*W1_pj*Yp24(:,j)/(conj(RHPp24(i))+RHPp24(j));
            
            W2_pi = evalfr(W2,RHPp24(i));
            W2_pj = evalfr(W2,RHPp24(j));
            Qp2(i,j) = ctranspose(Yp24(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp24(:,j)/(conj(RHPp24(i))+RHPp24(j));
        end
    end
    for i=1:nz24
        for j=1:np24
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz24(:,i))*inv(W1_zi)*W1_pj*Yp24(:,j)/(RHPz24(i)-RHPp24(j));
            
            W2_zi = evalfr(W2,RHPz24(i));
            W2_pj = evalfr(W2,RHPp24(j));
            Qzp2(i,j) = ctranspose(Yz24(:,i))*W2_zi*inv(W2_pj)*Yp24(:,j)/(RHPz24(i)-RHPp24(j));
        end
    end
    temp24 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSG24_min = sqrt(max(eig(temp24)));
end

%% Bounds on SG34 (W_in and P_in)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = G34ms(s)
sys34 = sys({'W_in','P_in'},'Z1');
l=2; % Number of outputs
RHPp34=p(find(p>0));
z34 = zero(sys34);
RHPz34=z34(find(z34>0));

if (isempty(RHPp34) || isempty(RHPz34))
    GammaSG34_min=0;
else
    np34=length(RHPp34); nz34=length(RHPz34);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp34 = zeros(l,np34);
    Yz34 = zeros(l,nz34);
    for i=1:np34
        Yp34(:,i) = Cs([3 4],:)*V(:,i)/norm(Cs([3 4],:)*V(:,i)); %Pole direction
    end
    for j=1:nz34
        [U,S,V]=svd(evalfr(sys34,RHPz34(j))); 
        Yz34(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gms3;Gms4];

    Qz1 = zeros(nz34,nz34); Qp1 = zeros(np34,np34); Qzp1 = zeros(nz34,np34); 
    Qz2 = zeros(nz34,nz34); Qp2 = zeros(np34,np34); Qzp2 = zeros(nz34,np34); 
    for i=1:nz34
        for j=1:nz34
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz34(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz34(:,j)/(RHPz34(i)+conj(RHPz34(j)));
            
            W2_zi = evalfr(W2,RHPz34(i));
            W2_zj = evalfr(W2,RHPz34(j));
            Qz2(i,j) = ctranspose(Yz34(:,i))*W2_zi*ctranspose(W2_zj)*Yz34(:,j)/(RHPz34(i)+conj(RHPz34(j)));
        end
    end
    for i=1:np34
        for j=1:np34
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp34(:,i))*ctranspose(W1)*W1_pj*Yp34(:,j)/(conj(RHPp34(i))+RHPp34(j));
            
            W2_pi = evalfr(W2,RHPp34(i));
            W2_pj = evalfr(W2,RHPp34(j));
            Qp2(i,j) = ctranspose(Yp34(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp34(:,j)/(conj(RHPp34(i))+RHPp34(j));
        end
    end
    for i=1:nz34
        for j=1:np34
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz34(:,i))*inv(W1_zi)*W1_pj*Yp34(:,j)/(RHPz34(i)-RHPp34(j));
            
            W2_zi = evalfr(W2,RHPz34(i));
            W2_pj = evalfr(W2,RHPp34(j));
            Qzp2(i,j) = ctranspose(Yz34(:,i))*W2_zi*inv(W2_pj)*Yp34(:,j)/(RHPz34(i)-RHPp34(j));
        end
    end
    temp34 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSG34_min = sqrt(max(eig(temp34)));
end

%% Bounds on SGd121 (P_bh,P_wh,d1)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms21;Gdms31]
sys12 = sys({'P_bh','P_wh'},'Z1');
l=2; % Number of outputs
RHPp12=p(find(p>0));
z12 = zero(sys12);
RHPz12=z12(find(z12>0));

if (isempty(RHPp12) || isempty(RHPz12))
    GammaSGd121_min=0;
else
    np12=length(RHPp12); nz12=length(RHPz12);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp12 = zeros(l,np12);
    Yz12 = zeros(l,nz12);
    for i=1:np12
        Yp12(:,i) = Cs([1 2],:)*V(:,i)/norm(Cs([1 2],:)*V(:,i)); %Pole direction
    end
    for j=1:nz12
        [U,S,V]=svd(evalfr(sys12,RHPz12(j))); 
        Yz12(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms11;Gdms21];

    Qz1 = zeros(nz12,nz12); Qp1 = zeros(np12,np12); Qzp1 = zeros(nz12,np12); 
    Qz2 = zeros(nz12,nz12); Qp2 = zeros(np12,np12); Qzp2 = zeros(nz12,np12); 
    for i=1:nz12
        for j=1:nz12
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz12(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz12(:,j)/(RHPz12(i)+conj(RHPz12(j)));
            
            W2_zi = evalfr(W2,RHPz12(i));
            W2_zj = evalfr(W2,RHPz12(j));
            Qz2(i,j) = ctranspose(Yz12(:,i))*W2_zi*ctranspose(W2_zj)*Yz12(:,j)/(RHPz12(i)+conj(RHPz12(j)));
        end
    end
    for i=1:np12
        for j=1:np12
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp12(:,i))*ctranspose(W1)*W1_pj*Yp12(:,j)/(conj(RHPp12(i))+RHPp12(j));
            
            W2_pi = evalfr(W2,RHPp12(i));
            W2_pj = evalfr(W2,RHPp12(j));
            Qp2(i,j) = ctranspose(Yp12(:,i))*ctranspose(pinv(W2_pi))*pinv(W2_pj)*Yp12(:,j)/(conj(RHPp12(i))+RHPp12(j));
        end
    end
    for i=1:nz12
        for j=1:np12
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz12(:,i))*inv(W1_zi)*W1_pj*Yp12(:,j)/(RHPz12(i)-RHPp12(j));
            
            W2_zi = evalfr(W2,RHPz12(i));
            W2_pj = evalfr(W2,RHPp12(j));
            Qzp2(i,j) = ctranspose(Yz12(:,i))*W2_zi*pinv(W2_pj)*Yp12(:,j)/(RHPz12(i)-RHPp12(j));
        end
    end
    temp12 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd121_min = real(sqrt(max(eig(temp12))));
end

%% Bounds on SGd131 (P_bh,W_in,d1)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms21;Gdms31]
sys13 = sys({'P_bh','W_in'},'Z1');
l=2; % Number of outputs
RHPp13=p(find(p>0));
z13 = zero(sys13);
RHPz13=z13(find(z13>0));

if (isempty(RHPp13) || isempty(RHPz13))
    GammaSGd131_min=0;
else
    np13=length(RHPp13); nz13=length(RHPz13);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp13 = zeros(l,np13);
    Yz13 = zeros(l,nz13);
    for i=1:np13
        Yp13(:,i) = Cs([1 3],:)*V(:,i)/norm(Cs([1 3],:)*V(:,i)); %Pole direction
    end
    for j=1:nz13
        [U,S,V]=svd(evalfr(sys13,RHPz13(j))); 
        Yz13(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms11;Gdms31];

    Qz1 = zeros(nz13,nz13); Qp1 = zeros(np13,np13); Qzp1 = zeros(nz13,np13); 
    Qz2 = zeros(nz13,nz13); Qp2 = zeros(np13,np13); Qzp2 = zeros(nz13,np13); 
    for i=1:nz13
        for j=1:nz13
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz13(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz13(:,j)/(RHPz13(i)+conj(RHPz13(j)));
            
            W2_zi = evalfr(W2,RHPz13(i));
            W2_zj = evalfr(W2,RHPz13(j));
            Qz2(i,j) = ctranspose(Yz13(:,i))*W2_zi*ctranspose(W2_zj)*Yz13(:,j)/(RHPz13(i)+conj(RHPz13(j)));
        end
    end
    for i=1:np13
        for j=1:np13
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp13(:,i))*ctranspose(W1)*W1_pj*Yp13(:,j)/(conj(RHPp13(i))+RHPp13(j));
            
            W2_pi = evalfr(W2,RHPp13(i));
            W2_pj = evalfr(W2,RHPp13(j));
            Qp2(i,j) = ctranspose(Yp13(:,i))*ctranspose(pinv(W2_pi))*pinv(W2_pj)*Yp13(:,j)/(conj(RHPp13(i))+RHPp13(j));
        end
    end
    for i=1:nz13
        for j=1:np13
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz13(:,i))*inv(W1_zi)*W1_pj*Yp13(:,j)/(RHPz13(i)-RHPp13(j));
            
            W2_zi = evalfr(W2,RHPz13(i));
            W2_pj = evalfr(W2,RHPp13(j));
            Qzp2(i,j) = ctranspose(Yz13(:,i))*W2_zi*pinv(W2_pj)*Yp13(:,j)/(RHPz13(i)-RHPp13(j));
        end
    end
    temp13 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd131_min = real(sqrt(max(eig(temp13))));
end

%% Bounds on SGd141 (P_bh,P_in,d1)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms21;Gdms31]
sys14 = sys({'P_bh','P_in'},'Z1');
l=2; % Number of outputs
RHPp14=p(find(p>0));
z14 = zero(sys14);
RHPz14=z14(find(z14>0));

if (isempty(RHPp14) || isempty(RHPz14))
    GammaSGd141_min=0;
else
    np14=length(RHPp14); nz14=length(RHPz14);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp14 = zeros(l,np14);
    Yz14 = zeros(l,nz14);
    for i=1:np14
        Yp14(:,i) = Cs([1 4],:)*V(:,i)/norm(Cs([1 4],:)*V(:,i)); %Pole direction
    end
    for j=1:nz14
        [U,S,V]=svd(evalfr(sys14,RHPz14(j))); 
        Yz14(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms11;Gdms41];

    Qz1 = zeros(nz14,nz14); Qp1 = zeros(np14,np14); Qzp1 = zeros(nz14,np14); 
    Qz2 = zeros(nz14,nz14); Qp2 = zeros(np14,np14); Qzp2 = zeros(nz14,np14); 
    for i=1:nz14
        for j=1:nz14
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz14(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz14(:,j)/(RHPz14(i)+conj(RHPz14(j)));
            
            W2_zi = evalfr(W2,RHPz14(i));
            W2_zj = evalfr(W2,RHPz14(j));
            Qz2(i,j) = ctranspose(Yz14(:,i))*W2_zi*ctranspose(W2_zj)*Yz14(:,j)/(RHPz14(i)+conj(RHPz14(j)));
        end
    end
    for i=1:np14
        for j=1:np14
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp14(:,i))*ctranspose(W1)*W1_pj*Yp14(:,j)/(conj(RHPp14(i))+RHPp14(j));
            
            W2_pi = evalfr(W2,RHPp14(i));
            W2_pj = evalfr(W2,RHPp14(j));
            Qp2(i,j) = ctranspose(Yp14(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp14(:,j)/(conj(RHPp14(i))+RHPp14(j));
        end
    end
    for i=1:nz14
        for j=1:np14
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz14(:,i))*inv(W1_zi)*W1_pj*Yp14(:,j)/(RHPz14(i)-RHPp14(j));
            
            W2_zi = evalfr(W2,RHPz14(i));
            W2_pj = evalfr(W2,RHPp14(j));
            Qzp2(i,j) = ctranspose(Yz14(:,i))*W2_zi*inv(W2_pj)*Yp14(:,j)/(RHPz14(i)-RHPp14(j));
        end
    end
    temp14 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd141_min = sqrt(max(eig(temp14)));
end

%% Bounds on SGd231 (P_wh,W_in,d1)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms21;Gdms31]
sys23 = sys({'P_wh','W_in'},'Z1');
l=2; % Number of outputs
RHPp23=p(find(p>0));
z23 = zero(sys23);
RHPz23=z23(find(z23>0));

if (isempty(RHPp23) || isempty(RHPz23))
    GammaSGd231_min=0;
else
    np23=length(RHPp23); nz23=length(RHPz23);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp23 = zeros(l,np23);
    Yz23 = zeros(l,nz23);
    for i=1:np23
        Yp23(:,i) = Cs([2 3],:)*V(:,i)/norm(Cs([2 3],:)*V(:,i)); %Pole direction
    end
    for j=1:nz23
        [U,S,V]=svd(evalfr(sys23,RHPz23(j))); 
        Yz23(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms21;Gdms31];

    Qz1 = zeros(nz23,nz23); Qp1 = zeros(np23,np23); Qzp1 = zeros(nz23,np23); 
    Qz2 = zeros(nz23,nz23); Qp2 = zeros(np23,np23); Qzp2 = zeros(nz23,np23); 
    for i=1:nz23
        for j=1:nz23
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz23(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz23(:,j)/(RHPz23(i)+conj(RHPz23(j)));
            
            W2_zi = evalfr(W2,RHPz23(i));
            W2_zj = evalfr(W2,RHPz23(j));
            Qz2(i,j) = ctranspose(Yz23(:,i))*W2_zi*ctranspose(W2_zj)*Yz23(:,j)/(RHPz23(i)+conj(RHPz23(j)));
        end
    end
    for i=1:np23
        for j=1:np23
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp23(:,i))*ctranspose(W1)*W1_pj*Yp23(:,j)/(conj(RHPp23(i))+RHPp23(j));
            
            W2_pi = evalfr(W2,RHPp23(i));
            W2_pj = evalfr(W2,RHPp23(j));
            Qp2(i,j) = ctranspose(Yp23(:,i))*ctranspose(pinv(W2_pi))*pinv(W2_pj)*Yp23(:,j)/(conj(RHPp23(i))+RHPp23(j));
        end
    end
    for i=1:nz23
        for j=1:np23
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz23(:,i))*inv(W1_zi)*W1_pj*Yp23(:,j)/(RHPz23(i)-RHPp23(j));
            
            W2_zi = evalfr(W2,RHPz23(i));
            W2_pj = evalfr(W2,RHPp23(j));
            Qzp2(i,j) = ctranspose(Yz23(:,i))*W2_zi*pinv(W2_pj)*Yp23(:,j)/(RHPz23(i)-RHPp23(j));
        end
    end
    temp23 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd231_min = real(sqrt(max(eig(temp23))));
end

%% Bounds on SGd241 (P_wh,P_in,d1)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms21;Gdms31]
sys24 = sys({'P_wh','P_in'},'Z1');
l=2; % Number of outputs
RHPp24=p(find(p>0));
z24 = zero(sys24);
RHPz24=z24(find(z24>0));

if (isempty(RHPp24) || isempty(RHPz24))
    GammaSGd241_min=0;
else
    np24=length(RHPp24); nz24=length(RHPz24);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp24 = zeros(l,np24);
    Yz24 = zeros(l,nz24);
    for i=1:np24
        Yp24(:,i) = Cs([2 4],:)*V(:,i)/norm(Cs([2 4],:)*V(:,i)); %Pole direction
    end
    for j=1:nz24
        [U,S,V]=svd(evalfr(sys24,RHPz24(j))); 
        Yz24(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms21;Gdms41];

    Qz1 = zeros(nz24,nz24); Qp1 = zeros(np24,np24); Qzp1 = zeros(nz24,np24); 
    Qz2 = zeros(nz24,nz24); Qp2 = zeros(np24,np24); Qzp2 = zeros(nz24,np24); 
    for i=1:nz24
        for j=1:nz24
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz24(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz24(:,j)/(RHPz24(i)+conj(RHPz24(j)));
            
            W2_zi = evalfr(W2,RHPz24(i));
            W2_zj = evalfr(W2,RHPz24(j));
            Qz2(i,j) = ctranspose(Yz24(:,i))*W2_zi*ctranspose(W2_zj)*Yz24(:,j)/(RHPz24(i)+conj(RHPz24(j)));
        end
    end
    for i=1:np24
        for j=1:np24
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp24(:,i))*ctranspose(W1)*W1_pj*Yp24(:,j)/(conj(RHPp24(i))+RHPp24(j));
            
            W2_pi = evalfr(W2,RHPp24(i));
            W2_pj = evalfr(W2,RHPp24(j));
            Qp2(i,j) = ctranspose(Yp24(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp24(:,j)/(conj(RHPp24(i))+RHPp24(j));
        end
    end
    for i=1:nz24
        for j=1:np24
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz24(:,i))*inv(W1_zi)*W1_pj*Yp24(:,j)/(RHPz24(i)-RHPp24(j));
            
            W2_zi = evalfr(W2,RHPz24(i));
            W2_pj = evalfr(W2,RHPp24(j));
            Qzp2(i,j) = ctranspose(Yz24(:,i))*W2_zi*inv(W2_pj)*Yp24(:,j)/(RHPz24(i)-RHPp24(j));
        end
    end
    temp24 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd241_min = sqrt(max(eig(temp24)));
end

%% Bounds on SGd341 (W_in,P_in,d1)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms21;Gdms31]
sys34 = sys({'W_in','P_in'},'Z1');
l=2; % Number of outputs
RHPp34=p(find(p>0));
z34 = zero(sys34);
RHPz34=z34(find(z34>0));

if (isempty(RHPp34) || isempty(RHPz34))
    GammaSGd341_min=0;
else
    np34=length(RHPp34); nz34=length(RHPz34);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp34 = zeros(l,np34);
    Yz34 = zeros(l,nz34);
    for i=1:np34
        Yp34(:,i) = Cs([3 4],:)*V(:,i)/norm(Cs([3 4],:)*V(:,i)); %Pole direction
    end
    for j=1:nz34
        [U,S,V]=svd(evalfr(sys34,RHPz34(j))); 
        Yz34(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms31;Gdms41];

    Qz1 = zeros(nz34,nz34); Qp1 = zeros(np34,np34); Qzp1 = zeros(nz34,np34); 
    Qz2 = zeros(nz34,nz34); Qp2 = zeros(np34,np34); Qzp2 = zeros(nz34,np34); 
    for i=1:nz34
        for j=1:nz34
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz34(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz34(:,j)/(RHPz34(i)+conj(RHPz34(j)));
            
            W2_zi = evalfr(W2,RHPz34(i));
            W2_zj = evalfr(W2,RHPz34(j));
            Qz2(i,j) = ctranspose(Yz34(:,i))*W2_zi*ctranspose(W2_zj)*Yz34(:,j)/(RHPz34(i)+conj(RHPz34(j)));
        end
    end
    for i=1:np34
        for j=1:np34
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp34(:,i))*ctranspose(W1)*W1_pj*Yp34(:,j)/(conj(RHPp34(i))+RHPp34(j));
            
            W2_pi = evalfr(W2,RHPp34(i));
            W2_pj = evalfr(W2,RHPp34(j));
            Qp2(i,j) = ctranspose(Yp34(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp34(:,j)/(conj(RHPp34(i))+RHPp34(j));
        end
    end
    for i=1:nz34
        for j=1:np34
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz34(:,i))*inv(W1_zi)*W1_pj*Yp34(:,j)/(RHPz34(i)-RHPp34(j));
            
            W2_zi = evalfr(W2,RHPz34(i));
            W2_pj = evalfr(W2,RHPp34(j));
            Qzp2(i,j) = ctranspose(Yz34(:,i))*W2_zi*inv(W2_pj)*Yp34(:,j)/(RHPz34(i)-RHPp34(j));
        end
    end
    temp34 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd341_min = sqrt(max(eig(temp34)));
end

%% Bounds on SGd122 (P_bh,P_wh,d2)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms22;Gdms32]
sys12 = sys({'P_bh','P_wh'},'Z1');
l=2; % Number of outputs
RHPp12=p(find(p>0));
z12 = zero(sys12);
RHPz12=z12(find(z12>0));

if (isempty(RHPp12) || isempty(RHPz12))
    GammaSGd122_min=0;
else
    np12=length(RHPp12); nz12=length(RHPz12);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp12 = zeros(l,np12);
    Yz12 = zeros(l,nz12);
    for i=1:np12
        Yp12(:,i) = Cs([1 2],:)*V(:,i)/norm(Cs([1 2],:)*V(:,i)); %Pole direction
    end
    for j=1:nz12
        [U,S,V]=svd(evalfr(sys12,RHPz12(j))); 
        Yz12(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms12;Gdms22];

    Qz1 = zeros(nz12,nz12); Qp1 = zeros(np12,np12); Qzp1 = zeros(nz12,np12); 
    Qz2 = zeros(nz12,nz12); Qp2 = zeros(np12,np12); Qzp2 = zeros(nz12,np12); 
    for i=1:nz12
        for j=1:nz12
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz12(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz12(:,j)/(RHPz12(i)+conj(RHPz12(j)));
            
            W2_zi = evalfr(W2,RHPz12(i));
            W2_zj = evalfr(W2,RHPz12(j));
            Qz2(i,j) = ctranspose(Yz12(:,i))*W2_zi*ctranspose(W2_zj)*Yz12(:,j)/(RHPz12(i)+conj(RHPz12(j)));
        end
    end
    for i=1:np12
        for j=1:np12
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp12(:,i))*ctranspose(W1)*W1_pj*Yp12(:,j)/(conj(RHPp12(i))+RHPp12(j));
            
            W2_pi = evalfr(W2,RHPp12(i));
            W2_pj = evalfr(W2,RHPp12(j));
            Qp2(i,j) = ctranspose(Yp12(:,i))*ctranspose(pinv(W2_pi))*pinv(W2_pj)*Yp12(:,j)/(conj(RHPp12(i))+RHPp12(j));
        end
    end
    for i=1:nz12
        for j=1:np12
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz12(:,i))*inv(W1_zi)*W1_pj*Yp12(:,j)/(RHPz12(i)-RHPp12(j));
            
            W2_zi = evalfr(W2,RHPz12(i));
            W2_pj = evalfr(W2,RHPp12(j));
            Qzp2(i,j) = ctranspose(Yz12(:,i))*W2_zi*pinv(W2_pj)*Yp12(:,j)/(RHPz12(i)-RHPp12(j));
        end
    end
    temp12 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd122_min = real(sqrt(max(eig(temp12))));
end

%% Bounds on SGd132 (P_bh,W_in,d2)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms22;Gdms32]
sys13 = sys({'P_bh','W_in'},'Z1');
l=2; % Number of outputs
RHPp13=p(find(p>0));
z13 = zero(sys13);
RHPz13=z13(find(z13>0));

if (isempty(RHPp13) || isempty(RHPz13))
    GammaSGd132_min=0;
else
    np13=length(RHPp13); nz13=length(RHPz13);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp13 = zeros(l,np13);
    Yz13 = zeros(l,nz13);
    for i=1:np13
        Yp13(:,i) = Cs([1 3],:)*V(:,i)/norm(Cs([1 3],:)*V(:,i)); %Pole direction
    end
    for j=1:nz13
        [U,S,V]=svd(evalfr(sys13,RHPz13(j))); 
        Yz13(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms12;Gdms32];

    Qz1 = zeros(nz13,nz13); Qp1 = zeros(np13,np13); Qzp1 = zeros(nz13,np13); 
    Qz2 = zeros(nz13,nz13); Qp2 = zeros(np13,np13); Qzp2 = zeros(nz13,np13); 
    for i=1:nz13
        for j=1:nz13
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz13(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz13(:,j)/(RHPz13(i)+conj(RHPz13(j)));
            
            W2_zi = evalfr(W2,RHPz13(i));
            W2_zj = evalfr(W2,RHPz13(j));
            Qz2(i,j) = ctranspose(Yz13(:,i))*W2_zi*ctranspose(W2_zj)*Yz13(:,j)/(RHPz13(i)+conj(RHPz13(j)));
        end
    end
    for i=1:np13
        for j=1:np13
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp13(:,i))*ctranspose(W1)*W1_pj*Yp13(:,j)/(conj(RHPp13(i))+RHPp13(j));
            
            W2_pi = evalfr(W2,RHPp13(i));
            W2_pj = evalfr(W2,RHPp13(j));
            Qp2(i,j) = ctranspose(Yp13(:,i))*ctranspose(pinv(W2_pi))*pinv(W2_pj)*Yp13(:,j)/(conj(RHPp13(i))+RHPp13(j));
        end
    end
    for i=1:nz13
        for j=1:np13
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz13(:,i))*inv(W1_zi)*W1_pj*Yp13(:,j)/(RHPz13(i)-RHPp13(j));
            
            W2_zi = evalfr(W2,RHPz13(i));
            W2_pj = evalfr(W2,RHPp13(j));
            Qzp2(i,j) = ctranspose(Yz13(:,i))*W2_zi*pinv(W2_pj)*Yp13(:,j)/(RHPz13(i)-RHPp13(j));
        end
    end
    temp13 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd132_min = real(sqrt(max(eig(temp13))));
end

%% Bounds on SGd141 (P_bh,P_in,d2)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms22;Gdms32]
sys14 = sys({'P_bh','P_in'},'Z1');
l=2; % Number of outputs
RHPp14=p(find(p>0));
z14 = zero(sys14);
RHPz14=z14(find(z14>0));

if (isempty(RHPp14) || isempty(RHPz14))
    GammaSGd142_min=0;
else
    np14=length(RHPp14); nz14=length(RHPz14);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp14 = zeros(l,np14);
    Yz14 = zeros(l,nz14);
    for i=1:np14
        Yp14(:,i) = Cs([1 4],:)*V(:,i)/norm(Cs([1 4],:)*V(:,i)); %Pole direction
    end
    for j=1:nz14
        [U,S,V]=svd(evalfr(sys14,RHPz14(j))); 
        Yz14(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms12;Gdms42];

    Qz1 = zeros(nz14,nz14); Qp1 = zeros(np14,np14); Qzp1 = zeros(nz14,np14); 
    Qz2 = zeros(nz14,nz14); Qp2 = zeros(np14,np14); Qzp2 = zeros(nz14,np14); 
    for i=1:nz14
        for j=1:nz14
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz14(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz14(:,j)/(RHPz14(i)+conj(RHPz14(j)));
            
            W2_zi = evalfr(W2,RHPz14(i));
            W2_zj = evalfr(W2,RHPz14(j));
            Qz2(i,j) = ctranspose(Yz14(:,i))*W2_zi*ctranspose(W2_zj)*Yz14(:,j)/(RHPz14(i)+conj(RHPz14(j)));
        end
    end
    for i=1:np14
        for j=1:np14
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp14(:,i))*ctranspose(W1)*W1_pj*Yp14(:,j)/(conj(RHPp14(i))+RHPp14(j));
            
            W2_pi = evalfr(W2,RHPp14(i));
            W2_pj = evalfr(W2,RHPp14(j));
            Qp2(i,j) = ctranspose(Yp14(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp14(:,j)/(conj(RHPp14(i))+RHPp14(j));
        end
    end
    for i=1:nz14
        for j=1:np14
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz14(:,i))*inv(W1_zi)*W1_pj*Yp14(:,j)/(RHPz14(i)-RHPp14(j));
            
            W2_zi = evalfr(W2,RHPz14(i));
            W2_pj = evalfr(W2,RHPp14(j));
            Qzp2(i,j) = ctranspose(Yz14(:,i))*W2_zi*inv(W2_pj)*Yp14(:,j)/(RHPz14(i)-RHPp14(j));
        end
    end
    temp14 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd142_min = sqrt(max(eig(temp14)));
end

%% Bounds on SGd231 (P_wh,W_in,d2)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms22;Gdms32]
sys23 = sys({'P_wh','W_in'},'Z1');
l=2; % Number of outputs
RHPp23=p(find(p>0));
z23 = zero(sys23);
RHPz23=z23(find(z23>0));

if (isempty(RHPp23) || isempty(RHPz23))
    GammaSGd232_min=0;
else
    np23=length(RHPp23); nz23=length(RHPz23);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp23 = zeros(l,np23);
    Yz23 = zeros(l,nz23);
    for i=1:np23
        Yp23(:,i) = Cs([2 3],:)*V(:,i)/norm(Cs([2 3],:)*V(:,i)); %Pole direction
    end
    for j=1:nz23
        [U,S,V]=svd(evalfr(sys23,RHPz23(j))); 
        Yz23(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms22;Gdms32];

    Qz1 = zeros(nz23,nz23); Qp1 = zeros(np23,np23); Qzp1 = zeros(nz23,np23); 
    Qz2 = zeros(nz23,nz23); Qp2 = zeros(np23,np23); Qzp2 = zeros(nz23,np23); 
    for i=1:nz23
        for j=1:nz23
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz23(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz23(:,j)/(RHPz23(i)+conj(RHPz23(j)));
            
            W2_zi = evalfr(W2,RHPz23(i));
            W2_zj = evalfr(W2,RHPz23(j));
            Qz2(i,j) = ctranspose(Yz23(:,i))*W2_zi*ctranspose(W2_zj)*Yz23(:,j)/(RHPz23(i)+conj(RHPz23(j)));
        end
    end
    for i=1:np23
        for j=1:np23
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp23(:,i))*ctranspose(W1)*W1_pj*Yp23(:,j)/(conj(RHPp23(i))+RHPp23(j));
            
            W2_pi = evalfr(W2,RHPp23(i));
            W2_pj = evalfr(W2,RHPp23(j));
            Qp2(i,j) = ctranspose(Yp23(:,i))*ctranspose(pinv(W2_pi))*pinv(W2_pj)*Yp23(:,j)/(conj(RHPp23(i))+RHPp23(j));
        end
    end
    for i=1:nz23
        for j=1:np23
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz23(:,i))*inv(W1_zi)*W1_pj*Yp23(:,j)/(RHPz23(i)-RHPp23(j));
            
            W2_zi = evalfr(W2,RHPz23(i));
            W2_pj = evalfr(W2,RHPp23(j));
            Qzp2(i,j) = ctranspose(Yz23(:,i))*W2_zi*pinv(W2_pj)*Yp23(:,j)/(RHPz23(i)-RHPp23(j));
        end
    end
    temp23 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd232_min = real(sqrt(max(eig(temp23))));
end

%% Bounds on SGd241 (P_wh,P_in,d2)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms22;Gdms32]
sys24 = sys({'P_wh','P_in'},'Z1');
l=2; % Number of outputs
RHPp24=p(find(p>0));
z24 = zero(sys24);
RHPz24=z24(find(z24>0));

if (isempty(RHPp24) || isempty(RHPz24))
    GammaSGd242_min=0;
else
    np24=length(RHPp24); nz24=length(RHPz24);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp24 = zeros(l,np24);
    Yz24 = zeros(l,nz24);
    for i=1:np24
        Yp24(:,i) = Cs([2 4],:)*V(:,i)/norm(Cs([2 4],:)*V(:,i)); %Pole direction
    end
    for j=1:nz24
        [U,S,V]=svd(evalfr(sys24,RHPz24(j))); 
        Yz24(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms22;Gdms42];

    Qz1 = zeros(nz24,nz24); Qp1 = zeros(np24,np24); Qzp1 = zeros(nz24,np24); 
    Qz2 = zeros(nz24,nz24); Qp2 = zeros(np24,np24); Qzp2 = zeros(nz24,np24); 
    for i=1:nz24
        for j=1:nz24
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz24(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz24(:,j)/(RHPz24(i)+conj(RHPz24(j)));
            
            W2_zi = evalfr(W2,RHPz24(i));
            W2_zj = evalfr(W2,RHPz24(j));
            Qz2(i,j) = ctranspose(Yz24(:,i))*W2_zi*ctranspose(W2_zj)*Yz24(:,j)/(RHPz24(i)+conj(RHPz24(j)));
        end
    end
    for i=1:np24
        for j=1:np24
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp24(:,i))*ctranspose(W1)*W1_pj*Yp24(:,j)/(conj(RHPp24(i))+RHPp24(j));
            
            W2_pi = evalfr(W2,RHPp24(i));
            W2_pj = evalfr(W2,RHPp24(j));
            Qp2(i,j) = ctranspose(Yp24(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp24(:,j)/(conj(RHPp24(i))+RHPp24(j));
        end
    end
    for i=1:nz24
        for j=1:np24
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz24(:,i))*inv(W1_zi)*W1_pj*Yp24(:,j)/(RHPz24(i)-RHPp24(j));
            
            W2_zi = evalfr(W2,RHPz24(i));
            W2_pj = evalfr(W2,RHPp24(j));
            Qzp2(i,j) = ctranspose(Yz24(:,i))*W2_zi*inv(W2_pj)*Yp24(:,j)/(RHPz24(i)-RHPp24(j));
        end
    end
    temp24 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd242_min = sqrt(max(eig(temp24)));
end

%% Bounds on SGd341 (W_in,P_in,d2)
% Page 226, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms22;Gdms32]
sys34 = sys({'W_in','P_in'},'Z1');
l=2; % Number of outputs
RHPp34=p(find(p>0));
z34 = zero(sys34);
RHPz34=z34(find(z34>0));

if (isempty(RHPp34) || isempty(RHPz34))
    GammaSGd342_min=0;
else
    np34=length(RHPp34); nz34=length(RHPz34);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp34 = zeros(l,np34);
    Yz34 = zeros(l,nz34);
    for i=1:np34
        Yp34(:,i) = Cs([3 4],:)*V(:,i)/norm(Cs([3 4],:)*V(:,i)); %Pole direction
    end
    for j=1:nz34
        [U,S,V]=svd(evalfr(sys34,RHPz34(j))); 
        Yz34(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms32;Gdms42];

    Qz1 = zeros(nz34,nz34); Qp1 = zeros(np34,np34); Qzp1 = zeros(nz34,np34); 
    Qz2 = zeros(nz34,nz34); Qp2 = zeros(np34,np34); Qzp2 = zeros(nz34,np34); 
    for i=1:nz34
        for j=1:nz34
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz34(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz34(:,j)/(RHPz34(i)+conj(RHPz34(j)));
            
            W2_zi = evalfr(W2,RHPz34(i));
            W2_zj = evalfr(W2,RHPz34(j));
            Qz2(i,j) = ctranspose(Yz34(:,i))*W2_zi*ctranspose(W2_zj)*Yz34(:,j)/(RHPz34(i)+conj(RHPz34(j)));
        end
    end
    for i=1:np34
        for j=1:np34
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp34(:,i))*ctranspose(W1)*W1_pj*Yp34(:,j)/(conj(RHPp34(i))+RHPp34(j));
            
            W2_pi = evalfr(W2,RHPp34(i));
            W2_pj = evalfr(W2,RHPp34(j));
            Qp2(i,j) = ctranspose(Yp34(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp34(:,j)/(conj(RHPp34(i))+RHPp34(j));
        end
    end
    for i=1:nz34
        for j=1:np34
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz34(:,i))*inv(W1_zi)*W1_pj*Yp34(:,j)/(RHPz34(i)-RHPp34(j));
            
            W2_zi = evalfr(W2,RHPz34(i));
            W2_pj = evalfr(W2,RHPp34(j));
            Qzp2(i,j) = ctranspose(Yz34(:,i))*W2_zi*inv(W2_pj)*Yp34(:,j)/(RHPz34(i)-RHPp34(j));
        end
    end
    temp34 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd342_min = sqrt(max(eig(temp34)));
end


%%
disp('Bounds on S and T');
disp(['The lowest achivable peak for S12 and T12 is: ', num2str(Msmin12)]);
disp(['The lowest achivable peak for S13 and T13 is: ', num2str(Msmin13)]);
disp(['The lowest achivable peak for S14 and T14 is: ', num2str(Msmin14)]);
disp(['The lowest achivable peak for S23 and T23 is: ', num2str(Msmin23)]);
disp(['The lowest achivable peak for S24 and T24 is: ', num2str(Msmin24)]);
disp(['The lowest achivable peak for S34 and T34 is: ', num2str(Msmin34)]);
disp('##########################################################');
disp('Bounds on KS');
disp(['The lowest achivable peak for KS12 is: ', num2str(KSmin12)]);
disp(['The lowest achivable peak for KS14 is: ', num2str(KSmin13)]);
disp(['The lowest achivable peak for KS14 is: ', num2str(KSmin14)]);
disp(['The lowest achivable peak for KS23 is: ', num2str(KSmin23)]);
disp(['The lowest achivable peak for KS24 is: ', num2str(KSmin24)]);
disp(['The lowest achivable peak for KS34 is: ', num2str(KSmin34)]);
disp('##########################################################');
disp('Bounds on KSGd form equation 6.26 (tight)');
disp(['The lowest achivable peak for KSGd121 (P_bh, P_wh, d1) is: ', num2str(KSGdmin121)]);
disp(['The lowest achivable peak for KSGd131 (P_bh, W_in, d1) is: ', num2str(KSGdmin131)]);
disp(['The lowest achivable peak for KSGd141 (P_bh, P_in, d1) is: ', num2str(KSGdmin141)]);
disp(['The lowest achivable peak for KSGd231 (P_wh, W_in, d1) is: ', num2str(KSGdmin231)]);
disp(['The lowest achivable peak for KSGd241 (P_wh, P_in, d1) is: ', num2str(KSGdmin241)]);
disp(['The lowest achivable peak for KSGd341 (W_in, P_in, d1) is: ', num2str(KSGdmin341)]);
disp('##########################################################');
disp(['The lowest achivable peak for KSGd122 (P_bh, P_wh, d2) is: ', num2str(KSGdmin122)]);
disp(['The lowest achivable peak for KSGd132 (P_bh, W_in, d2) is: ', num2str(KSGdmin132)]);
disp(['The lowest achivable peak for KSGd142 (P_bh, P_in, d2) is: ', num2str(KSGdmin142)]);
disp(['The lowest achivable peak for KSGd232 (P_wh, W_in, d2) is: ', num2str(KSGdmin232)]);
disp(['The lowest achivable peak for KSGd242 (P_wh, P_in, d2) is: ', num2str(KSGdmin242)]);
disp(['The lowest achivable peak for KSGd342 (W_in, P_in, d2) is: ', num2str(KSGdmin342)]);
disp('##########################################################');
disp('Bounds on SG');
disp(['The lowest achivable peak for SG12 (P_bh, P_wh) is: ', num2str(GammaSG12_min)]);
disp(['The lowest achivable peak for SG13 (P_bh, W_in) is: ', num2str(GammaSG13_min)]);
disp(['The lowest achivable peak for SG14 (P_bh, P_in) is: ', num2str(GammaSG14_min)]);
disp(['The lowest achivable peak for SG23 (P_wh, W_in) is: ', num2str(GammaSG23_min)]);
disp(['The lowest achivable peak for SG24 (P_wh, P_in) is: ', num2str(GammaSG24_min)]);
disp(['The lowest achivable peak for SG34 (W_in, P_in) is: ', num2str(GammaSG34_min)]);
disp('##########################################################');
disp('Bounds on SGd');
disp(['The lowest achivable peak for SGd121 (P_bh, P_wh, d1) is: ',num2str(GammaSGd121_min)]);
disp(['The lowest achivable peak for SGd131 (P_bh, W_in, d1) is: ',num2str(GammaSGd131_min)]);
disp(['The lowest achivable peak for SGd141 (P_bh, P_in, d1) is: ',num2str(GammaSGd141_min)]);
disp(['The lowest achivable peak for SGd231 (P_wh, W_in, d1) is: ',num2str(GammaSGd231_min)]);
disp(['The lowest achivable peak for SGd241 (P_wh, P_in, d1) is: ',num2str(GammaSGd241_min)]);
disp(['The lowest achivable peak for SGd341 (W_in, P_in, d1) is: ',num2str(GammaSGd341_min)]);
disp('##########################################################');
disp(['The lowest achivable peak for SGd122 (P_bh, P_wh, d2) is: ',num2str(GammaSGd122_min)]);
disp(['The lowest achivable peak for SGd132 (P_bh, W_in, d2) is: ',num2str(GammaSGd132_min)]);
disp(['The lowest achivable peak for SGd142 (P_bh, P_in, d2) is: ',num2str(GammaSGd142_min)]);
disp(['The lowest achivable peak for SGd232 (P_wh, W_in, d2) is: ',num2str(GammaSGd232_min)]);
disp(['The lowest achivable peak for SGd242 (P_wh, P_in, d2) is: ',num2str(GammaSGd242_min)]);
disp(['The lowest achivable peak for SGd342 (W_in, P_in, d2) is: ',num2str(GammaSGd342_min)]);

