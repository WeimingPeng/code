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
Gu = tf(sys({'P_bh','P_wh','W_in','P_in','P_rb','DP_r','P_t','Q_out','W_out','Rho_t','Alpha_L'},{'Z1','Z2'}));
Gd = tf(sys({'P_bh','P_wh','W_in','P_in','P_rb','DP_r','P_t','Q_out','W_out','Rho_t','Alpha_L'},{'d1','d2'}));
%% ****Tranfer function and steady state values ****
G1 =(Gu('P_bh',{'Z1','Z2'})); G2 =(Gu('P_wh',{'Z1','Z2'})); G3 =(Gu('W_in',{'Z1','Z2'})); G4 =(Gu('P_in',{'Z1','Z2'}));
G5 =(Gu('P_rb',{'Z1','Z2'})); G6 =(Gu('DP_r',{'Z1','Z2'})); G7 =(Gu('P_t',{'Z1','Z2'})); G8 =(Gu('Q_out',{'Z1','Z2'}));
G9 =(Gu('W_out',{'Z1','Z2'})); G10 =(Gu('Rho_t',{'Z1','Z2'})); G11 =(Gu('Alpha_L',{'Z1','Z2'}));
G1_0 = abs(freqresp(G1,0)); G2_0 = abs(freqresp(G2,0)); G3_0 = abs(freqresp(G3,0)); G4_0 = abs(freqresp(G4,0)); 
G5_0 = abs(freqresp(G5,0)); G6_0 = abs(freqresp(G6,0)); G7_0 = abs(freqresp(G7,0)); G8_0 = abs(freqresp(G8,0)); 
G9_0 = abs(freqresp(G9,0)); G10_0 = abs(freqresp(G10,0)); G11_0 = abs(freqresp(G11,0));
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
% Transfer functions from Inputs
[Gs11 Gm11 Gms11]=gmsi(G1(1));
[Gs12 Gm12 Gms12]=gmsi(G1(2));

[Gs21 Gm21 Gms21]=gmsi(G2(1));
[Gs22 Gm22 Gms22]=gmsi(G2(2));

[Gs31 Gm31 Gms31]=gmsi(G3(1));
[Gs32 Gm32 Gms32]=gmsi(G3(2));

[Gs41 Gm41 Gms41]=gmsi(G4(1));
[Gs42 Gm42 Gms42]=gmsi(G4(2));

[Gs51 Gm51 Gms51]=gmsi(G5(1));
[Gs52 Gm52 Gms52]=gmsi(G5(2));

[Gs61 Gm61 Gms61]=gmsi(G6(1));
[Gs62 Gm62 Gms62]=gmsi(G6(2));

[Gs71 Gm71 Gms71]=gmsi(G7(1));
[Gs72 Gm72 Gms72]=gmsi(G7(2));

[Gs81 Gm81 Gms81]=gmsi(G8(1));
[Gs82 Gm82 Gms82]=gmsi(G8(2));

[Gs91 Gm91 Gms91]=gmsi(G9(1));
[Gs92 Gm92 Gms92]=gmsi(G9(2));

[Gs101 Gm101 Gms101]=gmsi(G10(1));
[Gs102 Gm102 Gms102]=gmsi(G10(2));

[Gs111 Gm111 Gms111]=gmsi(G11(1));
[Gs112 Gm112 Gms112]=gmsi(G11(2));

% Transfer functions from Disturbances
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
%Sensitivity peak for combination of P_bh and P_in.
sys14 = sys({'P_bh','P_in'},{'Z1','Z2'});
G14 = tf(sys14);
p = pole(sys14);
z14 = zero(sys14);
RHPp14=p(find(p>0));
RHPz14=z14(find(z14>0));
l=2; % number of outputs
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

%% *** Bounds on KS ****
%To find the bound for KS the mirror image of the antistable part of the
%G`s are used. If a zero is very close to be negativ, this zero will be
%implemented in the antistable part.

%Lower bound on KS for the output P_bh & P_in
%Split G into stable and anti stable part, by using a matlab
%function. 
[GSu14, GNSu14]=stabsep(tf(sys14));

tf_mirror=(GNSu14)';
h14=hsvd(tf_mirror);
h14=h14(h14>0.01);
KSmin14=1/min(h14);

%% *** Bounds on KSGd **** form 6.26 page 230 (Tight))
%To find the bound for KSGd the mirror image of the antistable part of the
%1/(Gd,ms)*G`s are used. If a zero is very close to be negativ, this zero will be
%implemented in the antistable part.
s=tf('s');

%% Lower bound on KSGd1 for the outputs P_bh & P_in and disturbance d1
%Split G into stable and anti stable part, by using a matlab
%function. 
Gdms_inv = [1/Gdms11  1/Gdms41];
temp = minreal(Gdms_inv*G14);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu351 = gnsdij(temp);
[GSu141, GNSu141]=stabsep(temp/(1.e-3*s+1));

tf_mirror=(GNSu141)';
h141=hsvd(tf_mirror);
h141=h141(h141>0.01);
KSGdmin141=1/min(h141);
%% Lower bound on KSGd2 for the outputs P_bh & P_in and disturbance d2
%Split G into stable and anti stable part, by using a matlab
%function. 
Gdms_inv = [1/Gdms12  1/Gdms42];
temp = minreal(Gdms_inv*G14);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu352 = gnsdij(temp);
[GSu142, GNSu142]=stabsep(temp/(1.e-3*s+1));

tf_mirror=(GNSu142)';
h142=hsvd(tf_mirror);
h142=h142(h142>0.01);
KSGdmin142=1/min(h142);

%% *** Bounds on KSGd **** form 6.27 page 230 (Not tight)
% having only one input, input pole direction is always 1
%% Lower bound on KSGd141 (first disturbance) for the output P_bh & P_in and d1
[Q,Pi]=eig(A'); Q=Q(:,find(diag(Pi)>0));
UP = Bs1'*Q;
Gs14_inv = inv([Gs11 Gs12 ; Gs41 Gs42]);
Gdms141 = [Gdms11;Gdms41];
RHPp14=p(find(p>0));
np14 = length(RHPp14);
KSGd_min = zeros(1,np14);
for i = np14
    KSGd_min(i) = norm(ctranspose(UP)*evalfr(Gs14_inv,RHPp14(i))*evalfr(Gdms141,RHPp14(i)));
    %KSGd_min(i) = norm(evalfr(Gs14_inv,RHPp14(i))*evalfr(Gdms141,RHPp14(i)));
end
KSGd141_min = max(KSGd_min);
%% Lower bound on KSGd142 (second disturbance) for the output P_bh & P_in and d2
Gs14_inv = inv([Gs11 Gs12 ; Gs41 Gs42]);
Gdms142 = [Gdms12;Gdms42];
RHPp14=p(find(p>0));
np14 = length(RHPp14);
KSGd_min = zeros(1,np14);
for i = np14
    KSGd_min(i) = norm(evalfr(Gs14_inv,RHPp14(i))*evalfr(Gdms142,RHPp14(i)));
end
KSGd142_min = max(KSGd_min);
%% Bounds on SG14 (P_bh & P_in)
% Page 235, Bound on ||W1SW2|| with W1 = 1; and W2 = G14ms(s)
l=2; % Number of outputs

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
    
    W2 = [Gms11 Gms12 ; Gms41 Gms42];

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

%% Bounds on SGd141 (P_bh,P_in,d1)
% Page 235, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms21;Gdms21]
l=2; % Number of outputs
RHPp14=p(find(p>0));
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
            Qp2(i,j) = ctranspose(Yp14(:,i))*ctranspose(pinv(W2_pi))*pinv(W2_pj)*Yp14(:,j)/(conj(RHPp14(i))+RHPp14(j));
        end
    end
    for i=1:nz14
        for j=1:np14
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz14(:,i))*inv(W1_zi)*W1_pj*Yp14(:,j)/(RHPz14(i)-RHPp14(j));
            
            W2_zi = evalfr(W2,RHPz14(i));
            W2_pj = evalfr(W2,RHPp14(j));
            Qzp2(i,j) = ctranspose(Yz14(:,i))*W2_zi*pinv(W2_pj)*Yp14(:,j)/(RHPz14(i)-RHPp14(j));
        end
    end
    temp14 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd141_min = sqrt(max(eig(temp14)));
end

%% Bounds on SGd142 (P_bh,P_in,d2)
% Page 235, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms22;Gdms22]
RHPp14=p(find(p>0));
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
            Qp2(i,j) = ctranspose(Yp14(:,i))*ctranspose(pinv(W2_pi))*pinv(W2_pj)*Yp14(:,j)/(conj(RHPp14(i))+RHPp14(j));
        end
    end
    for i=1:nz14
        for j=1:np14
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz14(:,i))*pinv(W1_zi)*W1_pj*Yp14(:,j)/(RHPz14(i)-RHPp14(j));
            
            W2_zi = evalfr(W2,RHPz14(i));
            W2_pj = evalfr(W2,RHPp14(j));
            Qzp2(i,j) = ctranspose(Yz14(:,i))*W2_zi*pinv(W2_pj)*Yp14(:,j)/(RHPz14(i)-RHPp14(j));
        end
    end
    temp14 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd142_min = sqrt(max(eig(temp14)));
end

%% Pole Vectors
    np14=length(RHPp14); nz14=length(RHPz14);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp14 = zeros(l,np14);
    Yz14 = zeros(l,nz14);
    for i=1:np14
        Yp14(:,i) = Cs([1 4],:)*V(:,i); %Pole direction
    end
    Yp_max = max(max(abs(Yp14)),[],2);
%% ******** Controllability analysis*************
%*** Sensitivity and complementery sensitivity peak***
%Sensitivity peak for combination of P_bh and P_rb.
sys15 = sys({'P_bh','P_rb'},{'Z1','Z2'});
G15 = tf(sys15);
p = pole(sys15);
z15 = zero(sys15);
RHPp15=p(find(p>0));
RHPz15=z15(find(z15>0));
l=2; % number of outputs
if (isempty(RHPp15) || isempty(RHPz15))
    Msmin15=1;
else
    np15=length(RHPp15); nz15=length(RHPz15);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp15 = zeros(l,np15);
    Yz15 = zeros(l,nz15);
    for i=1:np15
        Yp15(:,i) = Cs([1 5],:)*V(:,i)/norm(Cs([1 5],:)*V(:,i)); %Pole direction
    end
    for j=1:nz15
        [U,S,V]=svd(evalfr(G15,RHPz15(j))); Yz15(:,j)=U(:,end); %zero direction
    end
    Qz15 = zeros(nz15,nz15); Qp15 = zeros(np15,np15); Qzp15 = zeros(nz15,np15); 
    for i=1:nz15
        for j=1:nz15
            Qz15(i,j) = ctranspose(Yz15(:,i))*Yz15(:,j)/(RHPz15(i)+conj(RHPz15(j)));
        end
    end
    for i=1:np15
        for j=1:np15
            Qp15(i,j) = ctranspose(Yp15(:,i))*Yp15(:,j)/(conj(RHPp15(i))+RHPp15(j));
        end
    end
    for i=1:nz15
        for j=1:np15
            Qzp15(i,j) = ctranspose(Yz15(:,i))*Yp15(:,j)/(RHPz15(i)-RHPp15(j));
        end
    end
%     Qp1=(Yp'*Yp).*(1./(diag(RHPp1')*ones(np1)+ones(np1)*diag(RHPp1)));
%     Qz1=(Yz'*Yz).*(1./(diag(RHPz1)*ones(nz1)+ones(nz1)*diag(RHPz1')));
%     Qzp1=(Yz'*Yp).*(1./(diag(RHPz1)*ones(nz1,np1)-ones(nz1,np1)*diag(RHPp1)));
    Msmin15=sqrt(1+norm(sqrtm(inv(Qz15))*Qzp15*sqrtm(inv(Qp15)))^2);
end

%% *** Bounds on KS ****
%To find the bound for KS the mirror image of the antistable part of the
%G`s are used. If a zero is very close to be negativ, this zero will be
%implemented in the antistable part.

%Lower bound on KS for the output P_bh & P_rb
%Split G into stable and anti stable part, by using a matlab
%function. 
[GSu15, GNSu15]=stabsep(tf(sys15));

tf_mirror=(GNSu15)';
h15=hsvd(tf_mirror);
h15=h15(h15>0.01);
KSmin15=1/min(h15);

%% *** Bounds on KSGd **** form 6.26 page 230 (Tight))
%To find the bound for KSGd the mirror image of the antistable part of the
%1/(Gd,ms)*G`s are used. If a zero is very close to be negativ, this zero will be
%implemented in the antistable part.
s=tf('s');

%% Lower bound on KSGd1 for the outputs P_bh & P_rb and disturbance d1
%Split G into stable and anti stable part, by using a matlab
%function. 
Gdms_inv = [1/Gdms11  1/Gdms51];
temp = minreal(Gdms_inv*G15);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu351 = gnsdij(temp);
[GSu151, GNSu151]=stabsep(temp/(1.e-3*s+1));

tf_mirror=(GNSu151)';
h151=hsvd(tf_mirror);
h151=h151(h151>0.01);
KSGdmin151=1/min(h151);
%% Lower bound on KSGd2 for the outputs P_bh & P_rb and disturbance d2
%Split G into stable and anti stable part, by using a matlab
%function. 
Gdms_inv = [1/Gdms12  1/Gdms52];
temp = minreal(Gdms_inv*G15);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu352 = gnsdij(temp);
[GSu152, GNSu152]=stabsep(temp/(1.e-3*s+1));

tf_mirror=(GNSu152)';
h152=hsvd(tf_mirror);
h152=h152(h152>0.01);
KSGdmin152=1/min(h152);

%% *** Bounds on KSGd **** form 6.27 page 230 (Not tight)
% having only one input, input pole direction is always 1
%% Lower bound on KSGd151 (first disturbance) for the output P_bh & P_rb and d1
[Q,Pi]=eig(A'); Q=Q(:,find(diag(Pi)>0));
UP = Bs1'*Q;
Gs15_inv = inv([Gs11 Gs12 ; Gs51 Gs52]);
Gdms151 = [Gdms11;Gdms51];
RHPp15=p(find(p>0));
np15 = length(RHPp15);
KSGd_min = zeros(1,np15);
for i = np15
    KSGd_min(i) = norm(ctranspose(UP)*evalfr(Gs15_inv,RHPp15(i))*evalfr(Gdms151,RHPp15(i)));
    %KSGd_min(i) = norm(evalfr(Gs15_inv,RHPp15(i))*evalfr(Gdms151,RHPp15(i)));
end
KSGd151_min = max(KSGd_min);
%% Lower bound on KSGd152 (second disturbance) for the output P_bh & P_rb and d2
Gs15_inv = inv([Gs11 Gs12 ; Gs51 Gs52]);
Gdms152 = [Gdms12;Gdms52];
RHPp15=p(find(p>0));
np15 = length(RHPp15);
KSGd_min = zeros(1,np15);
for i = np15
    KSGd_min(i) = norm(evalfr(Gs15_inv,RHPp15(i))*evalfr(Gdms152,RHPp15(i)));
end
KSGd152_min = max(KSGd_min);
%% Bounds on SG15 (P_bh & P_rb)
% Page 235, Bound on ||W1SW2|| with W1 = 1; and W2 = G15ms(s)
l=2; % Number of outputs

if (isempty(RHPp15) || isempty(RHPz15))
    GammaSG15_min=0;
else
    np15=length(RHPp15); nz15=length(RHPz15);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp15 = zeros(l,np15);
    Yz15 = zeros(l,nz15);
    for i=1:np15
        Yp15(:,i) = Cs([1 5],:)*V(:,i)/norm(Cs([1 5],:)*V(:,i)); %Pole direction
    end
    for j=1:nz15
        [U,S,V]=svd(evalfr(sys15,RHPz15(j))); 
        Yz15(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gms11 Gms12 ; Gms51 Gms52];

    Qz1 = zeros(nz15,nz15); Qp1 = zeros(np15,np15); Qzp1 = zeros(nz15,np15); 
    Qz2 = zeros(nz15,nz15); Qp2 = zeros(np15,np15); Qzp2 = zeros(nz15,np15); 
    for i=1:nz15
        for j=1:nz15
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz15(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz15(:,j)/(RHPz15(i)+conj(RHPz15(j)));
            
            W2_zi = evalfr(W2,RHPz15(i));
            W2_zj = evalfr(W2,RHPz15(j));
            Qz2(i,j) = ctranspose(Yz15(:,i))*W2_zi*ctranspose(W2_zj)*Yz15(:,j)/(RHPz15(i)+conj(RHPz15(j)));
        end
    end
    for i=1:np15
        for j=1:np15
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp15(:,i))*ctranspose(W1)*W1_pj*Yp15(:,j)/(conj(RHPp15(i))+RHPp15(j));
            
            W2_pi = evalfr(W2,RHPp15(i));
            W2_pj = evalfr(W2,RHPp15(j));
            Qp2(i,j) = ctranspose(Yp15(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp15(:,j)/(conj(RHPp15(i))+RHPp15(j));
        end
    end
    for i=1:nz15
        for j=1:np15
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz15(:,i))*inv(W1_zi)*W1_pj*Yp15(:,j)/(RHPz15(i)-RHPp15(j));
            
            W2_zi = evalfr(W2,RHPz15(i));
            W2_pj = evalfr(W2,RHPp15(j));
            Qzp2(i,j) = ctranspose(Yz15(:,i))*W2_zi*inv(W2_pj)*Yp15(:,j)/(RHPz15(i)-RHPp15(j));
        end
    end
    temp15 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSG15_min = real(sqrt(max(eig(temp15))));
end

%% Bounds on SGd151 (P_bh,P_rb,d1)
% Page 235, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms21;Gdms21]
l=2; % Number of outputs
RHPp15=p(find(p>0));
RHPz15=z15(find(z15>0));

if (isempty(RHPp15) || isempty(RHPz15))
    GammaSGd151_min=0;
else
    np15=length(RHPp15); nz15=length(RHPz15);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp15 = zeros(l,np15);
    Yz15 = zeros(l,nz15);
    for i=1:np15
        Yp15(:,i) = Cs([1 5],:)*V(:,i)/norm(Cs([1 5],:)*V(:,i)); %Pole direction
    end
    for j=1:nz15
        [U,S,V]=svd(evalfr(sys15,RHPz15(j))); 
        Yz15(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms11;Gdms51];

    Qz1 = zeros(nz15,nz15); Qp1 = zeros(np15,np15); Qzp1 = zeros(nz15,np15); 
    Qz2 = zeros(nz15,nz15); Qp2 = zeros(np15,np15); Qzp2 = zeros(nz15,np15); 
    for i=1:nz15
        for j=1:nz15
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz15(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz15(:,j)/(RHPz15(i)+conj(RHPz15(j)));
            
            W2_zi = evalfr(W2,RHPz15(i));
            W2_zj = evalfr(W2,RHPz15(j));
            Qz2(i,j) = ctranspose(Yz15(:,i))*W2_zi*ctranspose(W2_zj)*Yz15(:,j)/(RHPz15(i)+conj(RHPz15(j)));
        end
    end
    for i=1:np15
        for j=1:np15
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp15(:,i))*ctranspose(W1)*W1_pj*Yp15(:,j)/(conj(RHPp15(i))+RHPp15(j));
            
            W2_pi = evalfr(W2,RHPp15(i));
            W2_pj = evalfr(W2,RHPp15(j));
            Qp2(i,j) = ctranspose(Yp15(:,i))*ctranspose(pinv(W2_pi))*pinv(W2_pj)*Yp15(:,j)/(conj(RHPp15(i))+RHPp15(j));
        end
    end
    for i=1:nz15
        for j=1:np15
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz15(:,i))*inv(W1_zi)*W1_pj*Yp15(:,j)/(RHPz15(i)-RHPp15(j));
            
            W2_zi = evalfr(W2,RHPz15(i));
            W2_pj = evalfr(W2,RHPp15(j));
            Qzp2(i,j) = ctranspose(Yz15(:,i))*W2_zi*pinv(W2_pj)*Yp15(:,j)/(RHPz15(i)-RHPp15(j));
        end
    end
    temp15 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd151_min = real(sqrt(max(eig(temp15))));
end

%% Bounds on SGd152 (P_bh,P_rb,d2)
% Page 235, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms22;Gdms22]
RHPp15=p(find(p>0));
RHPz15=z15(find(z15>0));

if (isempty(RHPp15) || isempty(RHPz15))
    GammaSGd152_min=0;
else
    np15=length(RHPp15); nz15=length(RHPz15);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp15 = zeros(l,np15);
    Yz15 = zeros(l,nz15);
    for i=1:np15
        Yp15(:,i) = Cs([1 5],:)*V(:,i)/norm(Cs([1 5],:)*V(:,i)); %Pole direction
    end
    for j=1:nz15
        [U,S,V]=svd(evalfr(sys15,RHPz15(j))); 
        Yz15(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms12;Gdms52];

    Qz1 = zeros(nz15,nz15); Qp1 = zeros(np15,np15); Qzp1 = zeros(nz15,np15); 
    Qz2 = zeros(nz15,nz15); Qp2 = zeros(np15,np15); Qzp2 = zeros(nz15,np15); 
    for i=1:nz15
        for j=1:nz15
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz15(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz15(:,j)/(RHPz15(i)+conj(RHPz15(j)));
            
            W2_zi = evalfr(W2,RHPz15(i));
            W2_zj = evalfr(W2,RHPz15(j));
            Qz2(i,j) = ctranspose(Yz15(:,i))*W2_zi*ctranspose(W2_zj)*Yz15(:,j)/(RHPz15(i)+conj(RHPz15(j)));
        end
    end
    for i=1:np15
        for j=1:np15
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp15(:,i))*ctranspose(W1)*W1_pj*Yp15(:,j)/(conj(RHPp15(i))+RHPp15(j));
            
            W2_pi = evalfr(W2,RHPp15(i));
            W2_pj = evalfr(W2,RHPp15(j));
            Qp2(i,j) = ctranspose(Yp15(:,i))*ctranspose(pinv(W2_pi))*pinv(W2_pj)*Yp15(:,j)/(conj(RHPp15(i))+RHPp15(j));
        end
    end
    for i=1:nz15
        for j=1:np15
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz15(:,i))*pinv(W1_zi)*W1_pj*Yp15(:,j)/(RHPz15(i)-RHPp15(j));
            
            W2_zi = evalfr(W2,RHPz15(i));
            W2_pj = evalfr(W2,RHPp15(j));
            Qzp2(i,j) = ctranspose(Yz15(:,i))*W2_zi*pinv(W2_pj)*Yp15(:,j)/(RHPz15(i)-RHPp15(j));
        end
    end
    temp15 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd152_min = sqrt(max(eig(temp15)));
end

%% Pole Vectors
    np15=length(RHPp15); nz15=length(RHPz15);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp15 = zeros(l,np15);
    Yz15 = zeros(l,nz15);
    for i=1:np15
        Yp15(:,i) = Cs([1 5],:)*V(:,i); %Pole direction
    end
    Yp_max = max(max(abs(Yp15)),[],2);
%% ******** Controllability analysis*************
%*** Sensitivity and complementery sensitivity peak***
%Sensitivity peak for combination of P_bh and P_t.
sys17 = sys({'P_bh','P_t'},{'Z1','Z2'});
G17 = tf(sys17);
p = pole(sys17);
z17 = zero(sys17);
RHPp17=p(find(p>0));
RHPz17=z17(find(z17>0));
l=2; % number of outputs
if (isempty(RHPp17) || isempty(RHPz17))
    Msmin17=1;
else
    np17=length(RHPp17); nz17=length(RHPz17);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp17 = zeros(l,np17);
    Yz17 = zeros(l,nz17);
    for i=1:np17
        Yp17(:,i) = Cs([1 7],:)*V(:,i)/norm(Cs([1 7],:)*V(:,i)); %Pole direction
    end
    for j=1:nz17
        [U,S,V]=svd(evalfr(G17,RHPz17(j))); Yz17(:,j)=U(:,end); %zero direction
    end
    Qz17 = zeros(nz17,nz17); Qp17 = zeros(np17,np17); Qzp17 = zeros(nz17,np17); 
    for i=1:nz17
        for j=1:nz17
            Qz17(i,j) = ctranspose(Yz17(:,i))*Yz17(:,j)/(RHPz17(i)+conj(RHPz17(j)));
        end
    end
    for i=1:np17
        for j=1:np17
            Qp17(i,j) = ctranspose(Yp17(:,i))*Yp17(:,j)/(conj(RHPp17(i))+RHPp17(j));
        end
    end
    for i=1:nz17
        for j=1:np17
            Qzp17(i,j) = ctranspose(Yz17(:,i))*Yp17(:,j)/(RHPz17(i)-RHPp17(j));
        end
    end
%     Qp1=(Yp'*Yp).*(1./(diag(RHPp1')*ones(np1)+ones(np1)*diag(RHPp1)));
%     Qz1=(Yz'*Yz).*(1./(diag(RHPz1)*ones(nz1)+ones(nz1)*diag(RHPz1')));
%     Qzp1=(Yz'*Yp).*(1./(diag(RHPz1)*ones(nz1,np1)-ones(nz1,np1)*diag(RHPp1)));
    Msmin17=sqrt(1+norm(sqrtm(inv(Qz17))*Qzp17*sqrtm(inv(Qp17)))^2);
end

%% *** Bounds on KS ****
%To find the bound for KS the mirror image of the antistable part of the
%G`s are used. If a zero is very close to be negativ, this zero will be
%implemented in the antistable part.

%Lower bound on KS for the output P_bh & P_t
%Split G into stable and anti stable part, by using a matlab
%function. 
[GSu17, GNSu17]=stabsep(tf(sys17));

tf_mirror=(GNSu17)';
h17=hsvd(tf_mirror);
h17=h17(h17>0.01);
KSmin17=1/min(h17);

%% *** Bounds on KSGd **** form 6.26 page 230 (Tight))
%To find the bound for KSGd the mirror image of the antistable part of the
%1/(Gd,ms)*G`s are used. If a zero is very close to be negativ, this zero will be
%implemented in the antistable part.
s=tf('s');

%% Lower bound on KSGd1 for the outputs P_bh & P_t and disturbance d1
%Split G into stable and anti stable part, by using a matlab
%function. 
Gdms_inv = [1/Gdms11  1/Gdms71];
temp = minreal(Gdms_inv*G17);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu351 = gnsdij(temp);
[GSu171, GNSu171]=stabsep(temp/(1.e-3*s+1));

tf_mirror=(GNSu171)';
h171=hsvd(tf_mirror);
h171=h171(h171>0.01);
KSGdmin171=1/min(h171);
%% Lower bound on KSGd2 for the outputs P_bh & P_t and disturbance d2
%Split G into stable and anti stable part, by using a matlab
%function. 
Gdms_inv = [1/Gdms12  1/Gdms72];
temp = minreal(Gdms_inv*G17);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu352 = gnsdij(temp);
[GSu172, GNSu172]=stabsep(temp/(1.e-3*s+1));

tf_mirror=(GNSu172)';
h172=hsvd(tf_mirror);
h172=h172(h172>0.01);
KSGdmin172=1/min(h172);

%% *** Bounds on KSGd **** form 6.27 page 230 (Not tight)
% having only one input, input pole direction is always 1
%% Lower bound on KSGd171 (first disturbance) for the output P_bh & P_t and d1
[Q,Pi]=eig(A'); Q=Q(:,find(diag(Pi)>0));
UP = Bs1'*Q;
Gs17_inv = inv([Gs11 Gs12 ; Gs71 Gs72]);
Gdms171 = [Gdms11;Gdms71];
RHPp17=p(find(p>0));
np17 = length(RHPp17);
KSGd_min = zeros(1,np17);
for i = np17
    KSGd_min(i) = norm(ctranspose(UP)*evalfr(Gs17_inv,RHPp17(i))*evalfr(Gdms171,RHPp17(i)));
    %KSGd_min(i) = norm(evalfr(Gs17_inv,RHPp17(i))*evalfr(Gdms171,RHPp17(i)));
end
KSGd171_min = real(max(KSGd_min));
%% Lower bound on KSGd172 (second disturbance) for the output P_bh & P_t and d2
Gs17_inv = inv([Gs11 Gs12 ; Gs71 Gs72]);
Gdms172 = [Gdms12;Gdms72];
RHPp17=p(find(p>0));
np17 = length(RHPp17);
KSGd_min = zeros(1,np17);
for i = np17
    KSGd_min(i) = norm(evalfr(Gs17_inv,RHPp17(i))*evalfr(Gdms172,RHPp17(i)));
end
KSGd172_min = max(KSGd_min);
%% Bounds on SG17 (P_bh & P_t)
% Page 235, Bound on ||W1SW2|| with W1 = 1; and W2 = G17ms(s)
l=2; % Number of outputs

if (isempty(RHPp17) || isempty(RHPz17))
    GammaSG17_min=0;
else
    np17=length(RHPp17); nz17=length(RHPz17);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp17 = zeros(l,np17);
    Yz17 = zeros(l,nz17);
    for i=1:np17
        Yp17(:,i) = Cs([1 7],:)*V(:,i)/norm(Cs([1 7],:)*V(:,i)); %Pole direction
    end
    for j=1:nz17
        [U,S,V]=svd(evalfr(sys17,RHPz17(j))); 
        Yz17(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gms11 Gms12 ; Gms71 Gms72];

    Qz1 = zeros(nz17,nz17); Qp1 = zeros(np17,np17); Qzp1 = zeros(nz17,np17); 
    Qz2 = zeros(nz17,nz17); Qp2 = zeros(np17,np17); Qzp2 = zeros(nz17,np17); 
    for i=1:nz17
        for j=1:nz17
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz17(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz17(:,j)/(RHPz17(i)+conj(RHPz17(j)));
            
            W2_zi = evalfr(W2,RHPz17(i));
            W2_zj = evalfr(W2,RHPz17(j));
            Qz2(i,j) = ctranspose(Yz17(:,i))*W2_zi*ctranspose(W2_zj)*Yz17(:,j)/(RHPz17(i)+conj(RHPz17(j)));
        end
    end
    for i=1:np17
        for j=1:np17
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp17(:,i))*ctranspose(W1)*W1_pj*Yp17(:,j)/(conj(RHPp17(i))+RHPp17(j));
            
            W2_pi = evalfr(W2,RHPp17(i));
            W2_pj = evalfr(W2,RHPp17(j));
            Qp2(i,j) = ctranspose(Yp17(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp17(:,j)/(conj(RHPp17(i))+RHPp17(j));
        end
    end
    for i=1:nz17
        for j=1:np17
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz17(:,i))*inv(W1_zi)*W1_pj*Yp17(:,j)/(RHPz17(i)-RHPp17(j));
            
            W2_zi = evalfr(W2,RHPz17(i));
            W2_pj = evalfr(W2,RHPp17(j));
            Qzp2(i,j) = ctranspose(Yz17(:,i))*W2_zi*inv(W2_pj)*Yp17(:,j)/(RHPz17(i)-RHPp17(j));
        end
    end
    temp17 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSG17_min = real(sqrt(max(eig(temp17))));
end

%% Bounds on SGd171 (P_bh,P_t,d1)
% Page 235, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms21;Gdms21]
l=2; % Number of outputs
RHPp17=p(find(p>0));
RHPz17=z17(find(z17>0));

if (isempty(RHPp17) || isempty(RHPz17))
    GammaSGd171_min=0;
else
    np17=length(RHPp17); nz17=length(RHPz17);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp17 = zeros(l,np17);
    Yz17 = zeros(l,nz17);
    for i=1:np17
        Yp17(:,i) = Cs([1 7],:)*V(:,i)/norm(Cs([1 7],:)*V(:,i)); %Pole direction
    end
    for j=1:nz17
        [U,S,V]=svd(evalfr(sys17,RHPz17(j))); 
        Yz17(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms11;Gdms71];

    Qz1 = zeros(nz17,nz17); Qp1 = zeros(np17,np17); Qzp1 = zeros(nz17,np17); 
    Qz2 = zeros(nz17,nz17); Qp2 = zeros(np17,np17); Qzp2 = zeros(nz17,np17); 
    for i=1:nz17
        for j=1:nz17
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz17(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz17(:,j)/(RHPz17(i)+conj(RHPz17(j)));
            
            W2_zi = evalfr(W2,RHPz17(i));
            W2_zj = evalfr(W2,RHPz17(j));
            Qz2(i,j) = ctranspose(Yz17(:,i))*W2_zi*ctranspose(W2_zj)*Yz17(:,j)/(RHPz17(i)+conj(RHPz17(j)));
        end
    end
    for i=1:np17
        for j=1:np17
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp17(:,i))*ctranspose(W1)*W1_pj*Yp17(:,j)/(conj(RHPp17(i))+RHPp17(j));
            
            W2_pi = evalfr(W2,RHPp17(i));
            W2_pj = evalfr(W2,RHPp17(j));
            Qp2(i,j) = ctranspose(Yp17(:,i))*ctranspose(pinv(W2_pi))*pinv(W2_pj)*Yp17(:,j)/(conj(RHPp17(i))+RHPp17(j));
        end
    end
    for i=1:nz17
        for j=1:np17
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz17(:,i))*inv(W1_zi)*W1_pj*Yp17(:,j)/(RHPz17(i)-RHPp17(j));
            
            W2_zi = evalfr(W2,RHPz17(i));
            W2_pj = evalfr(W2,RHPp17(j));
            Qzp2(i,j) = ctranspose(Yz17(:,i))*W2_zi*pinv(W2_pj)*Yp17(:,j)/(RHPz17(i)-RHPp17(j));
        end
    end
    temp17 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd171_min = real(sqrt(max(eig(temp17))));
end

%% Bounds on SGd172 (P_bh,P_t,d2)
% Page 235, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms22;Gdms22]
RHPp17=p(find(p>0));
RHPz17=z17(find(z17>0));

if (isempty(RHPp17) || isempty(RHPz17))
    GammaSGd172_min=0;
else
    np17=length(RHPp17); nz17=length(RHPz17);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp17 = zeros(l,np17);
    Yz17 = zeros(l,nz17);
    for i=1:np17
        Yp17(:,i) = Cs([1 7],:)*V(:,i)/norm(Cs([1 7],:)*V(:,i)); %Pole direction
    end
    for j=1:nz17
        [U,S,V]=svd(evalfr(sys17,RHPz17(j))); 
        Yz17(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms12;Gdms72];

    Qz1 = zeros(nz17,nz17); Qp1 = zeros(np17,np17); Qzp1 = zeros(nz17,np17); 
    Qz2 = zeros(nz17,nz17); Qp2 = zeros(np17,np17); Qzp2 = zeros(nz17,np17); 
    for i=1:nz17
        for j=1:nz17
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz17(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz17(:,j)/(RHPz17(i)+conj(RHPz17(j)));
            
            W2_zi = evalfr(W2,RHPz17(i));
            W2_zj = evalfr(W2,RHPz17(j));
            Qz2(i,j) = ctranspose(Yz17(:,i))*W2_zi*ctranspose(W2_zj)*Yz17(:,j)/(RHPz17(i)+conj(RHPz17(j)));
        end
    end
    for i=1:np17
        for j=1:np17
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp17(:,i))*ctranspose(W1)*W1_pj*Yp17(:,j)/(conj(RHPp17(i))+RHPp17(j));
            
            W2_pi = evalfr(W2,RHPp17(i));
            W2_pj = evalfr(W2,RHPp17(j));
            Qp2(i,j) = ctranspose(Yp17(:,i))*ctranspose(pinv(W2_pi))*pinv(W2_pj)*Yp17(:,j)/(conj(RHPp17(i))+RHPp17(j));
        end
    end
    for i=1:nz17
        for j=1:np17
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz17(:,i))*pinv(W1_zi)*W1_pj*Yp17(:,j)/(RHPz17(i)-RHPp17(j));
            
            W2_zi = evalfr(W2,RHPz17(i));
            W2_pj = evalfr(W2,RHPp17(j));
            Qzp2(i,j) = ctranspose(Yz17(:,i))*W2_zi*pinv(W2_pj)*Yp17(:,j)/(RHPz17(i)-RHPp17(j));
        end
    end
    temp17 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd172_min = real(sqrt(max(eig(temp17))));
end

%% Pole Vectors
    np17=length(RHPp17); nz17=length(RHPz17);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp17 = zeros(l,np17);
    Yz17 = zeros(l,nz17);
    for i=1:np17
        Yp17(:,i) = Cs([1 7],:)*V(:,i); %Pole direction
    end
    Yp_max = max(max(abs(Yp17)),[],2);
    
    %% ******** Controllability analysis*************
%*** Sensitivity and complementery sensitivity peak***
%Sensitivity peak for combination of P_bh and W_out.
sys19 = sys({'P_bh','W_out'},{'Z1','Z2'});
G19 = tf(sys19);
p = pole(sys19);
z19 = zero(sys19);
RHPp19=p(find(p>0));
RHPz19=z19(find(z19>0));
l=2; % number of outputs
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

%% *** Bounds on KS ****
%To find the bound for KS the mirror image of the antistable part of the
%G`s are used. If a zero is very close to be negativ, this zero will be
%implemented in the antistable part.

%Lower bound on KS for the output P_bh & W_out
%Split G into stable and anti stable part, by using a matlab
%function. 
[GSu19, GNSu19]=stabsep(tf(sys19));

tf_mirror=(GNSu19)';
h19=hsvd(tf_mirror);
h19=h19(h19>0.01);
KSmin19=1/min(h19);

%% *** Bounds on KSGd **** form 6.26 page 230 (Tight))
%To find the bound for KSGd the mirror image of the antistable part of the
%1/(Gd,ms)*G`s are used. If a zero is very close to be negativ, this zero will be
%implemented in the antistable part.
s=tf('s');

%% Lower bound on KSGd1 for the outputs P_bh & W_out and disturbance d1
%Split G into stable and anti stable part, by using a matlab
%function. 
Gdms_inv = [1/Gdms11  1/Gdms91];
temp = minreal(Gdms_inv*G19);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu351 = gnsdij(temp);
[GSu191, GNSu191]=stabsep(temp/(1.e-3*s+1));

tf_mirror=(GNSu191)';
h191=hsvd(tf_mirror);
h191=h191(h191>0.01);
KSGdmin191=1/min(h191);
%% Lower bound on KSGd2 for the outputs P_bh & W_out and disturbance d2
%Split G into stable and anti stable part, by using a matlab
%function. 
Gdms_inv = [1/Gdms12  1/Gdms92];
temp = minreal(Gdms_inv*G19);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu352 = gnsdij(temp);
[GSu192, GNSu192]=stabsep(temp/(1.e-3*s+1));

tf_mirror=(GNSu192)';
h192=hsvd(tf_mirror);
h192=h192(h192>0.01);
KSGdmin192=1/min(h192);

%% *** Bounds on KSGd **** form 6.27 page 230 (Not tight)
% having only one input, input pole direction is always 1
%% Lower bound on KSGd191 (first disturbance) for the output P_bh & W_out and d1
[Q,Pi]=eig(A'); Q=Q(:,find(diag(Pi)>0));
UP = Bs1'*Q;
Gs19_inv = inv([Gs11 Gs12 ; Gs91 Gs92]);
Gdms191 = [Gdms11;Gdms91];
RHPp19=p(find(p>0));
np19 = length(RHPp19);
KSGd_min = zeros(1,np19);
for i = np19
    KSGd_min(i) = norm(ctranspose(UP)*evalfr(Gs19_inv,RHPp19(i))*evalfr(Gdms191,RHPp19(i)));
    %KSGd_min(i) = norm(evalfr(Gs19_inv,RHPp19(i))*evalfr(Gdms191,RHPp19(i)));
end
KSGd191_min = real(max(KSGd_min));
%% Lower bound on KSGd192 (second disturbance) for the output P_bh & W_out and d2
Gs19_inv = inv([Gs11 Gs12 ; Gs91 Gs92]);
Gdms192 = [Gdms12;Gdms92];
RHPp19=p(find(p>0));
np19 = length(RHPp19);
KSGd_min = zeros(1,np19);
for i = np19
    KSGd_min(i) = norm(evalfr(Gs19_inv,RHPp19(i))*evalfr(Gdms192,RHPp19(i)));
end
KSGd192_min = max(KSGd_min);
%% Bounds on SG19 (P_bh & W_out)
% Page 235, Bound on ||W1SW2|| with W1 = 1; and W2 = G19ms(s)
l=2; % Number of outputs

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
    
    W2 = [Gms11 Gms12 ; Gms91 Gms92];

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
    GammaSG19_min = real(sqrt(max(eig(temp19))));
end

%% Bounds on SGd191 (P_bh,W_out,d1)
% Page 235, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms21;Gdms21]
l=2; % Number of outputs
RHPp19=p(find(p>0));
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
            Qp2(i,j) = ctranspose(Yp19(:,i))*ctranspose(pinv(W2_pi))*pinv(W2_pj)*Yp19(:,j)/(conj(RHPp19(i))+RHPp19(j));
        end
    end
    for i=1:nz19
        for j=1:np19
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz19(:,i))*inv(W1_zi)*W1_pj*Yp19(:,j)/(RHPz19(i)-RHPp19(j));
            
            W2_zi = evalfr(W2,RHPz19(i));
            W2_pj = evalfr(W2,RHPp19(j));
            Qzp2(i,j) = ctranspose(Yz19(:,i))*W2_zi*pinv(W2_pj)*Yp19(:,j)/(RHPz19(i)-RHPp19(j));
        end
    end
    temp19 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd191_min = real(sqrt(max(eig(temp19))));
end

%% Bounds on SGd192 (P_bh,W_out,d2)
% Page 235, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms22;Gdms22]
RHPp19=p(find(p>0));
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
            Qp2(i,j) = ctranspose(Yp19(:,i))*ctranspose(pinv(W2_pi))*pinv(W2_pj)*Yp19(:,j)/(conj(RHPp19(i))+RHPp19(j));
        end
    end
    for i=1:nz19
        for j=1:np19
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz19(:,i))*pinv(W1_zi)*W1_pj*Yp19(:,j)/(RHPz19(i)-RHPp19(j));
            
            W2_zi = evalfr(W2,RHPz19(i));
            W2_pj = evalfr(W2,RHPp19(j));
            Qzp2(i,j) = ctranspose(Yz19(:,i))*W2_zi*pinv(W2_pj)*Yp19(:,j)/(RHPz19(i)-RHPp19(j));
        end
    end
    temp19 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd192_min = real(sqrt(max(eig(temp19))));
end

%% Pole Vectors
    np19=length(RHPp19); nz19=length(RHPz19);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp19 = zeros(l,np19);
    Yz19 = zeros(l,nz19);
    for i=1:np19
        Yp19(:,i) = Cs([1 9],:)*V(:,i); %Pole direction
    end
    Yp_max = max(max(abs(Yp19)),[],2);
    %% ******** Controllability analysis*************
%*** Sensitivity and complementery sensitivity peak***
%Sensitivity peak for combination of P_wh and P_rb.
sys25 = sys({'P_wh','P_rb'},{'Z1','Z2'});
G25 = tf(sys25);
p = pole(sys25);
z25 = zero(sys25);
RHPp25=p(find(p>0));
RHPz25=z25(find(z25>0));
l=2; % number of outputs
if (isempty(RHPp25) || isempty(RHPz25))
    Msmin25=1;
else
    np25=length(RHPp25); nz25=length(RHPz25);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp25 = zeros(l,np25);
    Yz25 = zeros(l,nz25);
    for i=1:np25
        Yp25(:,i) = Cs([2 5],:)*V(:,i)/norm(Cs([2 5],:)*V(:,i)); %Pole direction
    end
    for j=1:nz25
        [U,S,V]=svd(evalfr(G25,RHPz25(j))); Yz25(:,j)=U(:,end); %zero direction
    end
    Qz25 = zeros(nz25,nz25); Qp25 = zeros(np25,np25); Qzp25 = zeros(nz25,np25); 
    for i=1:nz25
        for j=1:nz25
            Qz25(i,j) = ctranspose(Yz25(:,i))*Yz25(:,j)/(RHPz25(i)+conj(RHPz25(j)));
        end
    end
    for i=1:np25
        for j=1:np25
            Qp25(i,j) = ctranspose(Yp25(:,i))*Yp25(:,j)/(conj(RHPp25(i))+RHPp25(j));
        end
    end
    for i=1:nz25
        for j=1:np25
            Qzp25(i,j) = ctranspose(Yz25(:,i))*Yp25(:,j)/(RHPz25(i)-RHPp25(j));
        end
    end
%     Qp1=(Yp'*Yp).*(1./(diag(RHPp1')*ones(np1)+ones(np1)*diag(RHPp1)));
%     Qz1=(Yz'*Yz).*(1./(diag(RHPz1)*ones(nz1)+ones(nz1)*diag(RHPz1')));
%     Qzp1=(Yz'*Yp).*(1./(diag(RHPz1)*ones(nz1,np1)-ones(nz1,np1)*diag(RHPp1)));
    Msmin25=sqrt(1+norm(sqrtm(inv(Qz25))*Qzp25*sqrtm(inv(Qp25)))^2);
end

%% *** Bounds on KS ****
%To find the bound for KS the mirror image of the antistable part of the
%G`s are used. If a zero is very close to be negativ, this zero will be
%implemented in the antistable part.

%Lower bound on KS for the output P_wh & P_rb
%Split G into stable and anti stable part, by using a matlab
%function. 
[GSu25, GNSu25]=stabsep(tf(sys25));

tf_mirror=(GNSu25)';
h25=hsvd(tf_mirror);
h25=h25(h25>0.01);
KSmin25=1/min(h25);

%% *** Bounds on KSGd **** form 6.26 page 230 (Tight))
%To find the bound for KSGd the mirror image of the antistable part of the
%1/(Gd,ms)*G`s are used. If a zero is very close to be negativ, this zero will be
%implemented in the antistable part.
s=tf('s');

%% Lower bound on KSGd1 for the outputs P_wh & P_rb and disturbance d1
%Split G into stable and anti stable part, by using a matlab
%function. 
Gdms_inv = [1/Gdms21  1/Gdms51];
temp = minreal(Gdms_inv*G25);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu351 = gnsdij(temp);
[GSu251, GNSu251]=stabsep(temp/(1.e-3*s+1));

tf_mirror=(GNSu251)';
h251=hsvd(tf_mirror);
h251=h251(h251>0.01);
KSGdmin251=1/min(h251);
%% Lower bound on KSGd2 for the outputs P_wh & P_rb and disturbance d2
%Split G into stable and anti stable part, by using a matlab
%function. 
Gdms_inv = [1/Gdms22  1/Gdms52];
temp = minreal(Gdms_inv*G25);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu352 = gnsdij(temp);
[GSu252, GNSu252]=stabsep(temp/(1.e-3*s+1));

tf_mirror=(GNSu252)';
h252=hsvd(tf_mirror);
h252=h252(h252>0.01);
KSGdmin252=1/min(h252);

%% *** Bounds on KSGd **** form 6.27 page 230 (Not tight)
% having only one input, input pole direction is always 1
%% Lower bound on KSGd251 (first disturbance) for the output P_wh & P_rb and d1
[Q,Pi]=eig(A'); Q=Q(:,find(diag(Pi)>0));
UP = Bs1'*Q;
Gs25_inv = inv([Gs21 Gs22 ; Gs51 Gs52]);
Gdms251 = [Gdms21;Gdms51];
RHPp25=p(find(p>0));
np25 = length(RHPp25);
KSGd_min = zeros(1,np25);
for i = np25
    KSGd_min(i) = norm(ctranspose(UP)*evalfr(Gs25_inv,RHPp25(i))*evalfr(Gdms251,RHPp25(i)));
    %KSGd_min(i) = norm(evalfr(Gs25_inv,RHPp25(i))*evalfr(Gdms251,RHPp25(i)));
end
KSGd251_min = max(KSGd_min);
%% Lower bound on KSGd252 (second disturbance) for the output P_wh & P_rb and d2
Gs25_inv = inv([Gs21 Gs22 ; Gs51 Gs52]);
Gdms252 = [Gdms22;Gdms52];
RHPp25=p(find(p>0));
np25 = length(RHPp25);
KSGd_min = zeros(1,np25);
for i = np25
    KSGd_min(i) = norm(evalfr(Gs25_inv,RHPp25(i))*evalfr(Gdms252,RHPp25(i)));
end
KSGd252_min = max(KSGd_min);
%% Bounds on SG25 (P_wh & P_rb)
% Page 235, Bound on ||W1SW2|| with W1 = 1; and W2 = G25ms(s)
l=2; % Number of outputs

if (isempty(RHPp25) || isempty(RHPz25))
    GammaSG25_min=0;
else
    np25=length(RHPp25); nz25=length(RHPz25);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp25 = zeros(l,np25);
    Yz25 = zeros(l,nz25);
    for i=1:np25
        Yp25(:,i) = Cs([2 5],:)*V(:,i)/norm(Cs([2 5],:)*V(:,i)); %Pole direction
    end
    for j=1:nz25
        [U,S,V]=svd(evalfr(sys25,RHPz25(j))); 
        Yz25(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gms21 Gms22 ; Gms51 Gms52];

    Qz1 = zeros(nz25,nz25); Qp1 = zeros(np25,np25); Qzp1 = zeros(nz25,np25); 
    Qz2 = zeros(nz25,nz25); Qp2 = zeros(np25,np25); Qzp2 = zeros(nz25,np25); 
    for i=1:nz25
        for j=1:nz25
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz25(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz25(:,j)/(RHPz25(i)+conj(RHPz25(j)));
            
            W2_zi = evalfr(W2,RHPz25(i));
            W2_zj = evalfr(W2,RHPz25(j));
            Qz2(i,j) = ctranspose(Yz25(:,i))*W2_zi*ctranspose(W2_zj)*Yz25(:,j)/(RHPz25(i)+conj(RHPz25(j)));
        end
    end
    for i=1:np25
        for j=1:np25
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp25(:,i))*ctranspose(W1)*W1_pj*Yp25(:,j)/(conj(RHPp25(i))+RHPp25(j));
            
            W2_pi = evalfr(W2,RHPp25(i));
            W2_pj = evalfr(W2,RHPp25(j));
            Qp2(i,j) = ctranspose(Yp25(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp25(:,j)/(conj(RHPp25(i))+RHPp25(j));
        end
    end
    for i=1:nz25
        for j=1:np25
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz25(:,i))*inv(W1_zi)*W1_pj*Yp25(:,j)/(RHPz25(i)-RHPp25(j));
            
            W2_zi = evalfr(W2,RHPz25(i));
            W2_pj = evalfr(W2,RHPp25(j));
            Qzp2(i,j) = ctranspose(Yz25(:,i))*W2_zi*inv(W2_pj)*Yp25(:,j)/(RHPz25(i)-RHPp25(j));
        end
    end
    temp25 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSG25_min = real(sqrt(max(eig(temp25))));
end

%% Bounds on SGd251 (P_wh,P_rb,d1)
% Page 235, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms21;Gdms21]
l=2; % Number of outputs
RHPp25=p(find(p>0));
RHPz25=z25(find(z25>0));

if (isempty(RHPp25) || isempty(RHPz25))
    GammaSGd251_min=0;
else
    np25=length(RHPp25); nz25=length(RHPz25);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp25 = zeros(l,np25);
    Yz25 = zeros(l,nz25);
    for i=1:np25
        Yp25(:,i) = Cs([2 5],:)*V(:,i)/norm(Cs([2 5],:)*V(:,i)); %Pole direction
    end
    for j=1:nz25
        [U,S,V]=svd(evalfr(sys25,RHPz25(j))); 
        Yz25(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms21;Gdms51];

    Qz1 = zeros(nz25,nz25); Qp1 = zeros(np25,np25); Qzp1 = zeros(nz25,np25); 
    Qz2 = zeros(nz25,nz25); Qp2 = zeros(np25,np25); Qzp2 = zeros(nz25,np25); 
    for i=1:nz25
        for j=1:nz25
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz25(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz25(:,j)/(RHPz25(i)+conj(RHPz25(j)));
            
            W2_zi = evalfr(W2,RHPz25(i));
            W2_zj = evalfr(W2,RHPz25(j));
            Qz2(i,j) = ctranspose(Yz25(:,i))*W2_zi*ctranspose(W2_zj)*Yz25(:,j)/(RHPz25(i)+conj(RHPz25(j)));
        end
    end
    for i=1:np25
        for j=1:np25
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp25(:,i))*ctranspose(W1)*W1_pj*Yp25(:,j)/(conj(RHPp25(i))+RHPp25(j));
            
            W2_pi = evalfr(W2,RHPp25(i));
            W2_pj = evalfr(W2,RHPp25(j));
            Qp2(i,j) = ctranspose(Yp25(:,i))*ctranspose(pinv(W2_pi))*pinv(W2_pj)*Yp25(:,j)/(conj(RHPp25(i))+RHPp25(j));
        end
    end
    for i=1:nz25
        for j=1:np25
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz25(:,i))*inv(W1_zi)*W1_pj*Yp25(:,j)/(RHPz25(i)-RHPp25(j));
            
            W2_zi = evalfr(W2,RHPz25(i));
            W2_pj = evalfr(W2,RHPp25(j));
            Qzp2(i,j) = ctranspose(Yz25(:,i))*W2_zi*pinv(W2_pj)*Yp25(:,j)/(RHPz25(i)-RHPp25(j));
        end
    end
    temp25 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd251_min = real(sqrt(max(eig(temp25))));
end

%% Bounds on SGd252 (P_wh,P_rb,d2)
% Page 235, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms22;Gdms22]
RHPp25=p(find(p>0));
RHPz25=z25(find(z25>0));

if (isempty(RHPp25) || isempty(RHPz25))
    GammaSGd252_min=0;
else
    np25=length(RHPp25); nz25=length(RHPz25);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp25 = zeros(l,np25);
    Yz25 = zeros(l,nz25);
    for i=1:np25
        Yp25(:,i) = Cs([2 5],:)*V(:,i)/norm(Cs([2 5],:)*V(:,i)); %Pole direction
    end
    for j=1:nz25
        [U,S,V]=svd(evalfr(sys25,RHPz25(j))); 
        Yz25(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms22;Gdms52];

    Qz1 = zeros(nz25,nz25); Qp1 = zeros(np25,np25); Qzp1 = zeros(nz25,np25); 
    Qz2 = zeros(nz25,nz25); Qp2 = zeros(np25,np25); Qzp2 = zeros(nz25,np25); 
    for i=1:nz25
        for j=1:nz25
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz25(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz25(:,j)/(RHPz25(i)+conj(RHPz25(j)));
            
            W2_zi = evalfr(W2,RHPz25(i));
            W2_zj = evalfr(W2,RHPz25(j));
            Qz2(i,j) = ctranspose(Yz25(:,i))*W2_zi*ctranspose(W2_zj)*Yz25(:,j)/(RHPz25(i)+conj(RHPz25(j)));
        end
    end
    for i=1:np25
        for j=1:np25
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp25(:,i))*ctranspose(W1)*W1_pj*Yp25(:,j)/(conj(RHPp25(i))+RHPp25(j));
            
            W2_pi = evalfr(W2,RHPp25(i));
            W2_pj = evalfr(W2,RHPp25(j));
            Qp2(i,j) = ctranspose(Yp25(:,i))*ctranspose(pinv(W2_pi))*pinv(W2_pj)*Yp25(:,j)/(conj(RHPp25(i))+RHPp25(j));
        end
    end
    for i=1:nz25
        for j=1:np25
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz25(:,i))*pinv(W1_zi)*W1_pj*Yp25(:,j)/(RHPz25(i)-RHPp25(j));
            
            W2_zi = evalfr(W2,RHPz25(i));
            W2_pj = evalfr(W2,RHPp25(j));
            Qzp2(i,j) = ctranspose(Yz25(:,i))*W2_zi*pinv(W2_pj)*Yp25(:,j)/(RHPz25(i)-RHPp25(j));
        end
    end
    temp25 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd252_min = sqrt(max(eig(temp25)));
end

%% Pole Vectors
    np25=length(RHPp25); nz25=length(RHPz25);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp25 = zeros(l,np25);
    Yz25 = zeros(l,nz25);
    for i=1:np25
        Yp25(:,i) = Cs([2 5],:)*V(:,i); %Pole direction
    end
    Yp_max = max(max(abs(Yp25)),[],2);
    %% ******** Controllability analysis*************
%*** Sensitivity and complementery sensitivity peak***
%Sensitivity peak for combination of P_wh and P_t.
sys27 = sys({'P_wh','P_t'},{'Z1','Z2'});
G27 = tf(sys27);
p = pole(sys27);
z27 = zero(sys27);
RHPp27=p(find(p>0));
RHPz27=z27(find(z27>0));
l=2; % number of outputs
if (isempty(RHPp27) || isempty(RHPz27))
    Msmin27=1;
else
    np27=length(RHPp27); nz27=length(RHPz27);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp27 = zeros(l,np27);
    Yz27 = zeros(l,nz27);
    for i=1:np27
        Yp27(:,i) = Cs([2 7],:)*V(:,i)/norm(Cs([2 7],:)*V(:,i)); %Pole direction
    end
    for j=1:nz27
        [U,S,V]=svd(evalfr(G27,RHPz27(j))); Yz27(:,j)=U(:,end); %zero direction
    end
    Qz27 = zeros(nz27,nz27); Qp27 = zeros(np27,np27); Qzp27 = zeros(nz27,np27); 
    for i=1:nz27
        for j=1:nz27
            Qz27(i,j) = ctranspose(Yz27(:,i))*Yz27(:,j)/(RHPz27(i)+conj(RHPz27(j)));
        end
    end
    for i=1:np27
        for j=1:np27
            Qp27(i,j) = ctranspose(Yp27(:,i))*Yp27(:,j)/(conj(RHPp27(i))+RHPp27(j));
        end
    end
    for i=1:nz27
        for j=1:np27
            Qzp27(i,j) = ctranspose(Yz27(:,i))*Yp27(:,j)/(RHPz27(i)-RHPp27(j));
        end
    end
%     Qp1=(Yp'*Yp).*(1./(diag(RHPp1')*ones(np1)+ones(np1)*diag(RHPp1)));
%     Qz1=(Yz'*Yz).*(1./(diag(RHPz1)*ones(nz1)+ones(nz1)*diag(RHPz1')));
%     Qzp1=(Yz'*Yp).*(1./(diag(RHPz1)*ones(nz1,np1)-ones(nz1,np1)*diag(RHPp1)));
    Msmin27=sqrt(1+norm(sqrtm(inv(Qz27))*Qzp27*sqrtm(inv(Qp27)))^2);
end

%% *** Bounds on KS ****
%To find the bound for KS the mirror image of the antistable part of the
%G`s are used. If a zero is very close to be negativ, this zero will be
%implemented in the antistable part.

%Lower bound on KS for the output P_wh & P_t
%Split G into stable and anti stable part, by using a matlab
%function. 
[GSu27, GNSu27]=stabsep(tf(sys27));

tf_mirror=(GNSu27)';
h27=hsvd(tf_mirror);
h27=h27(h27>0.01);
KSmin27=1/min(h27);

%% *** Bounds on KSGd **** form 6.26 page 230 (Tight))
%To find the bound for KSGd the mirror image of the antistable part of the
%1/(Gd,ms)*G`s are used. If a zero is very close to be negativ, this zero will be
%implemented in the antistable part.
s=tf('s');

%% Lower bound on KSGd1 for the outputs P_wh & P_t and disturbance d1
%Split G into stable and anti stable part, by using a matlab
%function. 
Gdms_inv = [1/Gdms21  1/Gdms71];
temp = minreal(Gdms_inv*G27);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu351 = gnsdij(temp);
[GSu271, GNSu271]=stabsep(temp/(1.e-3*s+1));

tf_mirror=(GNSu271)';
h271=hsvd(tf_mirror);
h271=h271(h271>0.01);
KSGdmin271=1/min(h271);
%% Lower bound on KSGd2 for the outputs P_wh & P_t and disturbance d2
%Split G into stable and anti stable part, by using a matlab
%function. 
Gdms_inv = [1/Gdms22  1/Gdms72];
temp = minreal(Gdms_inv*G27);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu352 = gnsdij(temp);
[GSu272, GNSu272]=stabsep(temp/(1.e-3*s+1));

tf_mirror=(GNSu272)';
h272=hsvd(tf_mirror);
h272=h272(h272>0.01);
KSGdmin272=1/min(h272);

%% *** Bounds on KSGd **** form 6.27 page 230 (Not tight)
% having only one input, input pole direction is always 1
%% Lower bound on KSGd271 (first disturbance) for the output P_wh & P_t and d1
[Q,Pi]=eig(A'); Q=Q(:,find(diag(Pi)>0));
UP = Bs1'*Q;
Gs27_inv = inv([Gs21 Gs22 ; Gs71 Gs72]);
Gdms271 = [Gdms21;Gdms71];
RHPp27=p(find(p>0));
np27 = length(RHPp27);
KSGd_min = zeros(1,np27);
for i = np27
    KSGd_min(i) = norm(ctranspose(UP)*evalfr(Gs27_inv,RHPp27(i))*evalfr(Gdms271,RHPp27(i)));
    %KSGd_min(i) = norm(evalfr(Gs27_inv,RHPp27(i))*evalfr(Gdms271,RHPp27(i)));
end
KSGd271_min = real(max(KSGd_min));
%% Lower bound on KSGd272 (second disturbance) for the output P_wh & P_t and d2
Gs27_inv = inv([Gs21 Gs22 ; Gs71 Gs72]);
Gdms272 = [Gdms22;Gdms72];
RHPp27=p(find(p>0));
np27 = length(RHPp27);
KSGd_min = zeros(1,np27);
for i = np27
    KSGd_min(i) = norm(evalfr(Gs27_inv,RHPp27(i))*evalfr(Gdms272,RHPp27(i)));
end
KSGd272_min = max(KSGd_min);
%% Bounds on SG27 (P_wh & P_t)
% Page 235, Bound on ||W1SW2|| with W1 = 1; and W2 = G27ms(s)
l=2; % Number of outputs

if (isempty(RHPp27) || isempty(RHPz27))
    GammaSG27_min=0;
else
    np27=length(RHPp27); nz27=length(RHPz27);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp27 = zeros(l,np27);
    Yz27 = zeros(l,nz27);
    for i=1:np27
        Yp27(:,i) = Cs([2 7],:)*V(:,i)/norm(Cs([2 7],:)*V(:,i)); %Pole direction
    end
    for j=1:nz27
        [U,S,V]=svd(evalfr(sys27,RHPz27(j))); 
        Yz27(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gms21 Gms22 ; Gms71 Gms72];

    Qz1 = zeros(nz27,nz27); Qp1 = zeros(np27,np27); Qzp1 = zeros(nz27,np27); 
    Qz2 = zeros(nz27,nz27); Qp2 = zeros(np27,np27); Qzp2 = zeros(nz27,np27); 
    for i=1:nz27
        for j=1:nz27
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz27(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz27(:,j)/(RHPz27(i)+conj(RHPz27(j)));
            
            W2_zi = evalfr(W2,RHPz27(i));
            W2_zj = evalfr(W2,RHPz27(j));
            Qz2(i,j) = ctranspose(Yz27(:,i))*W2_zi*ctranspose(W2_zj)*Yz27(:,j)/(RHPz27(i)+conj(RHPz27(j)));
        end
    end
    for i=1:np27
        for j=1:np27
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp27(:,i))*ctranspose(W1)*W1_pj*Yp27(:,j)/(conj(RHPp27(i))+RHPp27(j));
            
            W2_pi = evalfr(W2,RHPp27(i));
            W2_pj = evalfr(W2,RHPp27(j));
            Qp2(i,j) = ctranspose(Yp27(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp27(:,j)/(conj(RHPp27(i))+RHPp27(j));
        end
    end
    for i=1:nz27
        for j=1:np27
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz27(:,i))*inv(W1_zi)*W1_pj*Yp27(:,j)/(RHPz27(i)-RHPp27(j));
            
            W2_zi = evalfr(W2,RHPz27(i));
            W2_pj = evalfr(W2,RHPp27(j));
            Qzp2(i,j) = ctranspose(Yz27(:,i))*W2_zi*inv(W2_pj)*Yp27(:,j)/(RHPz27(i)-RHPp27(j));
        end
    end
    temp27 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSG27_min = real(sqrt(max(eig(temp27))));
end

%% Bounds on SGd271 (P_wh,P_t,d1)
% Page 235, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms21;Gdms21]
l=2; % Number of outputs
RHPp27=p(find(p>0));
RHPz27=z27(find(z27>0));

if (isempty(RHPp27) || isempty(RHPz27))
    GammaSGd271_min=0;
else
    np27=length(RHPp27); nz27=length(RHPz27);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp27 = zeros(l,np27);
    Yz27 = zeros(l,nz27);
    for i=1:np27
        Yp27(:,i) = Cs([2 7],:)*V(:,i)/norm(Cs([2 7],:)*V(:,i)); %Pole direction
    end
    for j=1:nz27
        [U,S,V]=svd(evalfr(sys27,RHPz27(j))); 
        Yz27(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms21;Gdms71];

    Qz1 = zeros(nz27,nz27); Qp1 = zeros(np27,np27); Qzp1 = zeros(nz27,np27); 
    Qz2 = zeros(nz27,nz27); Qp2 = zeros(np27,np27); Qzp2 = zeros(nz27,np27); 
    for i=1:nz27
        for j=1:nz27
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz27(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz27(:,j)/(RHPz27(i)+conj(RHPz27(j)));
            
            W2_zi = evalfr(W2,RHPz27(i));
            W2_zj = evalfr(W2,RHPz27(j));
            Qz2(i,j) = ctranspose(Yz27(:,i))*W2_zi*ctranspose(W2_zj)*Yz27(:,j)/(RHPz27(i)+conj(RHPz27(j)));
        end
    end
    for i=1:np27
        for j=1:np27
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp27(:,i))*ctranspose(W1)*W1_pj*Yp27(:,j)/(conj(RHPp27(i))+RHPp27(j));
            
            W2_pi = evalfr(W2,RHPp27(i));
            W2_pj = evalfr(W2,RHPp27(j));
            Qp2(i,j) = ctranspose(Yp27(:,i))*ctranspose(pinv(W2_pi))*pinv(W2_pj)*Yp27(:,j)/(conj(RHPp27(i))+RHPp27(j));
        end
    end
    for i=1:nz27
        for j=1:np27
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz27(:,i))*inv(W1_zi)*W1_pj*Yp27(:,j)/(RHPz27(i)-RHPp27(j));
            
            W2_zi = evalfr(W2,RHPz27(i));
            W2_pj = evalfr(W2,RHPp27(j));
            Qzp2(i,j) = ctranspose(Yz27(:,i))*W2_zi*pinv(W2_pj)*Yp27(:,j)/(RHPz27(i)-RHPp27(j));
        end
    end
    temp27 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd271_min = real(sqrt(max(eig(temp27))));
end

%% Bounds on SGd272 (P_wh,P_t,d2)
% Page 235, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms22;Gdms22]
RHPp27=p(find(p>0));
RHPz27=z27(find(z27>0));

if (isempty(RHPp27) || isempty(RHPz27))
    GammaSGd272_min=0;
else
    np27=length(RHPp27); nz27=length(RHPz27);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp27 = zeros(l,np27);
    Yz27 = zeros(l,nz27);
    for i=1:np27
        Yp27(:,i) = Cs([2 7],:)*V(:,i)/norm(Cs([2 7],:)*V(:,i)); %Pole direction
    end
    for j=1:nz27
        [U,S,V]=svd(evalfr(sys27,RHPz27(j))); 
        Yz27(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms22;Gdms72];

    Qz1 = zeros(nz27,nz27); Qp1 = zeros(np27,np27); Qzp1 = zeros(nz27,np27); 
    Qz2 = zeros(nz27,nz27); Qp2 = zeros(np27,np27); Qzp2 = zeros(nz27,np27); 
    for i=1:nz27
        for j=1:nz27
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz27(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz27(:,j)/(RHPz27(i)+conj(RHPz27(j)));
            
            W2_zi = evalfr(W2,RHPz27(i));
            W2_zj = evalfr(W2,RHPz27(j));
            Qz2(i,j) = ctranspose(Yz27(:,i))*W2_zi*ctranspose(W2_zj)*Yz27(:,j)/(RHPz27(i)+conj(RHPz27(j)));
        end
    end
    for i=1:np27
        for j=1:np27
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp27(:,i))*ctranspose(W1)*W1_pj*Yp27(:,j)/(conj(RHPp27(i))+RHPp27(j));
            
            W2_pi = evalfr(W2,RHPp27(i));
            W2_pj = evalfr(W2,RHPp27(j));
            Qp2(i,j) = ctranspose(Yp27(:,i))*ctranspose(pinv(W2_pi))*pinv(W2_pj)*Yp27(:,j)/(conj(RHPp27(i))+RHPp27(j));
        end
    end
    for i=1:nz27
        for j=1:np27
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz27(:,i))*pinv(W1_zi)*W1_pj*Yp27(:,j)/(RHPz27(i)-RHPp27(j));
            
            W2_zi = evalfr(W2,RHPz27(i));
            W2_pj = evalfr(W2,RHPp27(j));
            Qzp2(i,j) = ctranspose(Yz27(:,i))*W2_zi*pinv(W2_pj)*Yp27(:,j)/(RHPz27(i)-RHPp27(j));
        end
    end
    temp27 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd272_min = real(sqrt(max(eig(temp27))));
end

%% Pole Vectors
    np27=length(RHPp27); nz27=length(RHPz27);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp27 = zeros(l,np27);
    Yz27 = zeros(l,nz27);
    for i=1:np27
        Yp27(:,i) = Cs([2 7],:)*V(:,i); %Pole direction
    end
    Yp_max = max(max(abs(Yp27)),[],2);
    
    %% ******** Controllability analysis*************
%*** Sensitivity and complementery sensitivity peak***
%Sensitivity peak for combination of P_wh and W_out.
sys29 = sys({'P_wh','W_out'},{'Z1','Z2'});
G29 = tf(sys29);
p = pole(sys29);
z29 = zero(sys29);
RHPp29=p(find(p>0));
RHPz29=z29(find(z29>0));
l=2; % number of outputs
if (isempty(RHPp29) || isempty(RHPz29))
    Msmin29=1;
else
    np29=length(RHPp29); nz29=length(RHPz29);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp29 = zeros(l,np29);
    Yz29 = zeros(l,nz29);
    for i=1:np29
        Yp29(:,i) = Cs([2 9],:)*V(:,i)/norm(Cs([2 9],:)*V(:,i)); %Pole direction
    end
    for j=1:nz29
        [U,S,V]=svd(evalfr(G29,RHPz29(j))); Yz29(:,j)=U(:,end); %zero direction
    end
    Qz29 = zeros(nz29,nz29); Qp29 = zeros(np29,np29); Qzp29 = zeros(nz29,np29); 
    for i=1:nz29
        for j=1:nz29
            Qz29(i,j) = ctranspose(Yz29(:,i))*Yz29(:,j)/(RHPz29(i)+conj(RHPz29(j)));
        end
    end
    for i=1:np29
        for j=1:np29
            Qp29(i,j) = ctranspose(Yp29(:,i))*Yp29(:,j)/(conj(RHPp29(i))+RHPp29(j));
        end
    end
    for i=1:nz29
        for j=1:np29
            Qzp29(i,j) = ctranspose(Yz29(:,i))*Yp29(:,j)/(RHPz29(i)-RHPp29(j));
        end
    end
%     Qp1=(Yp'*Yp).*(1./(diag(RHPp1')*ones(np1)+ones(np1)*diag(RHPp1)));
%     Qz1=(Yz'*Yz).*(1./(diag(RHPz1)*ones(nz1)+ones(nz1)*diag(RHPz1')));
%     Qzp1=(Yz'*Yp).*(1./(diag(RHPz1)*ones(nz1,np1)-ones(nz1,np1)*diag(RHPp1)));
    Msmin29=sqrt(1+norm(sqrtm(inv(Qz29))*Qzp29*sqrtm(inv(Qp29)))^2);
end

%% *** Bounds on KS ****
%To find the bound for KS the mirror image of the antistable part of the
%G`s are used. If a zero is very close to be negativ, this zero will be
%implemented in the antistable part.

%Lower bound on KS for the output P_wh & W_out
%Split G into stable and anti stable part, by using a matlab
%function. 
[GSu29, GNSu29]=stabsep(tf(sys29));

tf_mirror=(GNSu29)';
h29=hsvd(tf_mirror);
h29=h29(h29>0.01);
KSmin29=1/min(h29);

%% *** Bounds on KSGd **** form 6.26 page 230 (Tight))
%To find the bound for KSGd the mirror image of the antistable part of the
%1/(Gd,ms)*G`s are used. If a zero is very close to be negativ, this zero will be
%implemented in the antistable part.
s=tf('s');

%% Lower bound on KSGd1 for the outputs P_wh & W_out and disturbance d1
%Split G into stable and anti stable part, by using a matlab
%function. 
Gdms_inv = [1/Gdms21  1/Gdms91];
temp = minreal(Gdms_inv*G29);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu351 = gnsdij(temp);
[GSu291, GNSu291]=stabsep(temp/(1.e-3*s+1));

tf_mirror=(GNSu291)';
h291=hsvd(tf_mirror);
h291=h291(h291>0.01);
KSGdmin291=1/min(h291);
%% Lower bound on KSGd2 for the outputs P_wh & W_out and disturbance d2
%Split G into stable and anti stable part, by using a matlab
%function. 
Gdms_inv = [1/Gdms22  1/Gdms92];
temp = minreal(Gdms_inv*G29);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu352 = gnsdij(temp);
[GSu292, GNSu292]=stabsep(temp/(1.e-3*s+1));

tf_mirror=(GNSu292)';
h292=hsvd(tf_mirror);
h292=h292(h292>0.01);
KSGdmin292=1/min(h292);

%% *** Bounds on KSGd **** form 6.27 page 230 (Not tight)
% having only one input, input pole direction is always 1
%% Lower bound on KSGd291 (first disturbance) for the output P_wh & W_out and d1
[Q,Pi]=eig(A'); Q=Q(:,find(diag(Pi)>0));
UP = Bs1'*Q;
Gs29_inv = inv([Gs21 Gs22 ; Gs91 Gs92]);
Gdms291 = [Gdms21;Gdms91];
RHPp29=p(find(p>0));
np29 = length(RHPp29);
KSGd_min = zeros(1,np29);
for i = np29
    KSGd_min(i) = norm(ctranspose(UP)*evalfr(Gs29_inv,RHPp29(i))*evalfr(Gdms291,RHPp29(i)));
    %KSGd_min(i) = norm(evalfr(Gs29_inv,RHPp29(i))*evalfr(Gdms291,RHPp29(i)));
end
KSGd291_min = real(max(KSGd_min));
%% Lower bound on KSGd292 (second disturbance) for the output P_wh & W_out and d2
Gs29_inv = inv([Gs21 Gs22 ; Gs91 Gs92]);
Gdms292 = [Gdms22;Gdms92];
RHPp29=p(find(p>0));
np29 = length(RHPp29);
KSGd_min = zeros(1,np29);
for i = np29
    KSGd_min(i) = norm(evalfr(Gs29_inv,RHPp29(i))*evalfr(Gdms292,RHPp29(i)));
end
KSGd292_min = max(KSGd_min);
%% Bounds on SG29 (P_wh & W_out)
% Page 235, Bound on ||W1SW2|| with W1 = 1; and W2 = G29ms(s)
l=2; % Number of outputs

if (isempty(RHPp29) || isempty(RHPz29))
    GammaSG29_min=0;
else
    np29=length(RHPp29); nz29=length(RHPz29);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp29 = zeros(l,np29);
    Yz29 = zeros(l,nz29);
    for i=1:np29
        Yp29(:,i) = Cs([2 9],:)*V(:,i)/norm(Cs([2 9],:)*V(:,i)); %Pole direction
    end
    for j=1:nz29
        [U,S,V]=svd(evalfr(sys29,RHPz29(j))); 
        Yz29(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gms21 Gms22 ; Gms91 Gms92];

    Qz1 = zeros(nz29,nz29); Qp1 = zeros(np29,np29); Qzp1 = zeros(nz29,np29); 
    Qz2 = zeros(nz29,nz29); Qp2 = zeros(np29,np29); Qzp2 = zeros(nz29,np29); 
    for i=1:nz29
        for j=1:nz29
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz29(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz29(:,j)/(RHPz29(i)+conj(RHPz29(j)));
            
            W2_zi = evalfr(W2,RHPz29(i));
            W2_zj = evalfr(W2,RHPz29(j));
            Qz2(i,j) = ctranspose(Yz29(:,i))*W2_zi*ctranspose(W2_zj)*Yz29(:,j)/(RHPz29(i)+conj(RHPz29(j)));
        end
    end
    for i=1:np29
        for j=1:np29
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp29(:,i))*ctranspose(W1)*W1_pj*Yp29(:,j)/(conj(RHPp29(i))+RHPp29(j));
            
            W2_pi = evalfr(W2,RHPp29(i));
            W2_pj = evalfr(W2,RHPp29(j));
            Qp2(i,j) = ctranspose(Yp29(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp29(:,j)/(conj(RHPp29(i))+RHPp29(j));
        end
    end
    for i=1:nz29
        for j=1:np29
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz29(:,i))*inv(W1_zi)*W1_pj*Yp29(:,j)/(RHPz29(i)-RHPp29(j));
            
            W2_zi = evalfr(W2,RHPz29(i));
            W2_pj = evalfr(W2,RHPp29(j));
            Qzp2(i,j) = ctranspose(Yz29(:,i))*W2_zi*inv(W2_pj)*Yp29(:,j)/(RHPz29(i)-RHPp29(j));
        end
    end
    temp29 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSG29_min = real(sqrt(max(eig(temp29))));
end

%% Bounds on SGd291 (P_wh,W_out,d1)
% Page 235, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms21;Gdms21]
l=2; % Number of outputs
RHPp29=p(find(p>0));
RHPz29=z29(find(z29>0));

if (isempty(RHPp29) || isempty(RHPz29))
    GammaSGd291_min=0;
else
    np29=length(RHPp29); nz29=length(RHPz29);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp29 = zeros(l,np29);
    Yz29 = zeros(l,nz29);
    for i=1:np29
        Yp29(:,i) = Cs([2 9],:)*V(:,i)/norm(Cs([2 9],:)*V(:,i)); %Pole direction
    end
    for j=1:nz29
        [U,S,V]=svd(evalfr(sys29,RHPz29(j))); 
        Yz29(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms21;Gdms91];

    Qz1 = zeros(nz29,nz29); Qp1 = zeros(np29,np29); Qzp1 = zeros(nz29,np29); 
    Qz2 = zeros(nz29,nz29); Qp2 = zeros(np29,np29); Qzp2 = zeros(nz29,np29); 
    for i=1:nz29
        for j=1:nz29
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz29(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz29(:,j)/(RHPz29(i)+conj(RHPz29(j)));
            
            W2_zi = evalfr(W2,RHPz29(i));
            W2_zj = evalfr(W2,RHPz29(j));
            Qz2(i,j) = ctranspose(Yz29(:,i))*W2_zi*ctranspose(W2_zj)*Yz29(:,j)/(RHPz29(i)+conj(RHPz29(j)));
        end
    end
    for i=1:np29
        for j=1:np29
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp29(:,i))*ctranspose(W1)*W1_pj*Yp29(:,j)/(conj(RHPp29(i))+RHPp29(j));
            
            W2_pi = evalfr(W2,RHPp29(i));
            W2_pj = evalfr(W2,RHPp29(j));
            Qp2(i,j) = ctranspose(Yp29(:,i))*ctranspose(pinv(W2_pi))*pinv(W2_pj)*Yp29(:,j)/(conj(RHPp29(i))+RHPp29(j));
        end
    end
    for i=1:nz29
        for j=1:np29
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz29(:,i))*inv(W1_zi)*W1_pj*Yp29(:,j)/(RHPz29(i)-RHPp29(j));
            
            W2_zi = evalfr(W2,RHPz29(i));
            W2_pj = evalfr(W2,RHPp29(j));
            Qzp2(i,j) = ctranspose(Yz29(:,i))*W2_zi*pinv(W2_pj)*Yp29(:,j)/(RHPz29(i)-RHPp29(j));
        end
    end
    temp29 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd291_min = real(sqrt(max(eig(temp29))));
end

%% Bounds on SGd292 (P_wh,W_out,d2)
% Page 235, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms22;Gdms22]
RHPp29=p(find(p>0));
RHPz29=z29(find(z29>0));

if (isempty(RHPp29) || isempty(RHPz29))
    GammaSGd292_min=0;
else
    np29=length(RHPp29); nz29=length(RHPz29);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp29 = zeros(l,np29);
    Yz29 = zeros(l,nz29);
    for i=1:np29
        Yp29(:,i) = Cs([2 9],:)*V(:,i)/norm(Cs([2 9],:)*V(:,i)); %Pole direction
    end
    for j=1:nz29
        [U,S,V]=svd(evalfr(sys29,RHPz29(j))); 
        Yz29(:,j)=U(:,end); %zero direction
    end
    W1 = 1;
    
    W2 = [Gdms22;Gdms92];

    Qz1 = zeros(nz29,nz29); Qp1 = zeros(np29,np29); Qzp1 = zeros(nz29,np29); 
    Qz2 = zeros(nz29,nz29); Qp2 = zeros(np29,np29); Qzp2 = zeros(nz29,np29); 
    for i=1:nz29
        for j=1:nz29
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz29(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz29(:,j)/(RHPz29(i)+conj(RHPz29(j)));
            
            W2_zi = evalfr(W2,RHPz29(i));
            W2_zj = evalfr(W2,RHPz29(j));
            Qz2(i,j) = ctranspose(Yz29(:,i))*W2_zi*ctranspose(W2_zj)*Yz29(:,j)/(RHPz29(i)+conj(RHPz29(j)));
        end
    end
    for i=1:np29
        for j=1:np29
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp29(:,i))*ctranspose(W1)*W1_pj*Yp29(:,j)/(conj(RHPp29(i))+RHPp29(j));
            
            W2_pi = evalfr(W2,RHPp29(i));
            W2_pj = evalfr(W2,RHPp29(j));
            Qp2(i,j) = ctranspose(Yp29(:,i))*ctranspose(pinv(W2_pi))*pinv(W2_pj)*Yp29(:,j)/(conj(RHPp29(i))+RHPp29(j));
        end
    end
    for i=1:nz29
        for j=1:np29
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz29(:,i))*pinv(W1_zi)*W1_pj*Yp29(:,j)/(RHPz29(i)-RHPp29(j));
            
            W2_zi = evalfr(W2,RHPz29(i));
            W2_pj = evalfr(W2,RHPp29(j));
            Qzp2(i,j) = ctranspose(Yz29(:,i))*W2_zi*pinv(W2_pj)*Yp29(:,j)/(RHPz29(i)-RHPp29(j));
        end
    end
    temp29 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd292_min = real(sqrt(max(eig(temp29))));
end

%% Pole Vectors
    np29=length(RHPp29); nz29=length(RHPz29);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp29 = zeros(l,np29);
    Yz29 = zeros(l,nz29);
    for i=1:np29
        Yp29(:,i) = Cs([2 9],:)*V(:,i); %Pole direction
    end
    Yp_max = max(max(abs(Yp29)),[],2);
%% ******** Controllability analysis*************
%*** Sensitivity and complementery sensitivity peak***
%Sensitivity peak for combination of P_in and P_t.
sys47 = sys({'P_in','P_t'},{'Z1','Z2'});
G47 = tf(sys47);
p = pole(sys47);
z47 = zero(sys47);
RHPp47=p(find(p>0));
RHPz47=z47(find(z47>0));
l=2; % number of outputs
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

%% *** Bounds on KS ****
%To find the bound for KS the mirror image of the antistable part of the
%G`s are used. If a zero is very close to be negativ, this zero will be
%implemented in the antistable part.

%Lower bound on KS for the output P_in & P_t
%Split G into stable and anti stable part, by using a matlab
%function. 
[GSu47, GNSu47]=stabsep(tf(sys47));

tf_mirror=(GNSu47)';
h47=hsvd(tf_mirror);
h47=h47(h47>0.01);
KSmin47=1/min(h47);

%% *** Bounds on KSGd **** form 6.26 page 230 (Tight))
%To find the bound for KSGd the mirror image of the antistable part of the
%1/(Gd,ms)*G`s are used. If a zero is very close to be negativ, this zero will be
%implemented in the antistable part.
s=tf('s');

%% Lower bound on KSGd1 for the outputs P_in & P_t and disturbance d1
%Split G into stable and anti stable part, by using a matlab
%function. 
Gdms_inv = [1/Gdms41  1/Gdms71];
temp = minreal(Gdms_inv*G47);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu351 = gnsdij(temp);
[GSu471, GNSu471]=stabsep(temp/(1.e-3*s+1));

tf_mirror=(GNSu471)';
h471=hsvd(tf_mirror);
h471=h471(h471>0.01);
KSGdmin471=1/min(h471);
%% Lower bound on KSGd2 for the outputs P_in & P_t and disturbance d2
%Split G into stable and anti stable part, by using a matlab
%function. 
Gdms_inv = [1/Gdms42  1/Gdms72];
temp = minreal(Gdms_inv*G47);
%Split G into stable and anti stable part, by using a matlab
%function. 
%GNSu352 = gnsdij(temp);
[GSu472, GNSu472]=stabsep(temp/(1.e-3*s+1));

tf_mirror=(GNSu472)';
h472=hsvd(tf_mirror);
h472=h472(h472>0.01);
KSGdmin472=1/min(h472);

%% *** Bounds on KSGd **** form 6.27 page 230 (Not tight)
% having only one input, input pole direction is always 1
%% Lower bound on KSGd471 (first disturbance) for the output P_in & P_t and d1
[Q,Pi]=eig(A'); Q=Q(:,find(diag(Pi)>0));
UP = Bs1'*Q;
Gs47_inv = inv([Gs41 Gs42 ; Gs71 Gs72]);
Gdms471 = [Gdms41;Gdms71];
RHPp47=p(find(p>0));
np47 = length(RHPp47);
KSGd_min = zeros(1,np47);
for i = np47
    KSGd_min(i) = norm(ctranspose(UP)*evalfr(Gs47_inv,RHPp47(i))*evalfr(Gdms471,RHPp47(i)));
    %KSGd_min(i) = norm(evalfr(Gs47_inv,RHPp47(i))*evalfr(Gdms471,RHPp47(i)));
end
KSGd471_min = max(KSGd_min);
%% Lower bound on KSGd472 (second disturbance) for the output P_in & P_t and d2
Gs47_inv = inv([Gs41 Gs42 ; Gs71 Gs72]);
Gdms472 = [Gdms42;Gdms72];
RHPp47=p(find(p>0));
np47 = length(RHPp47);
KSGd_min = zeros(1,np47);
for i = np47
    KSGd_min(i) = norm(evalfr(Gs47_inv,RHPp47(i))*evalfr(Gdms472,RHPp47(i)));
end
KSGd472_min = max(KSGd_min);
%% Bounds on SG47 (P_in & P_t)
% Page 235, Bound on ||W1SW2|| with W1 = 1; and W2 = G47ms(s)
l=2; % Number of outputs

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
    
    W2 = [Gms41 Gms42 ; Gms71 Gms72];

    Qz1 = zeros(nz47,nz47); Qp1 = zeros(np47,np47); Qzp1 = zeros(nz47,np47); 
    Qz2 = zeros(nz47,nz47); Qp2 = zeros(np47,np47); Qzp2 = zeros(nz47,np47); 
    for i=1:nz47
        for j=1:nz47
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz47(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz47(:,j)/(RHPz47(i)+conj(RHPz47(j)));
            
            W2_zi = evalfr(W2,RHPz47(i));
            W2_zj = evalfr(W2,RHPz47(j));
            Qz2(i,j) = ctranspose(Yz47(:,i))*W2_zi*ctranspose(W2_zj)*Yz47(:,j)/(RHPz47(i)+conj(RHPz47(j)));
        end
    end
    for i=1:np47
        for j=1:np47
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp47(:,i))*ctranspose(W1)*W1_pj*Yp47(:,j)/(conj(RHPp47(i))+RHPp47(j));
            
            W2_pi = evalfr(W2,RHPp47(i));
            W2_pj = evalfr(W2,RHPp47(j));
            Qp2(i,j) = ctranspose(Yp47(:,i))*ctranspose(inv(W2_pi))*inv(W2_pj)*Yp47(:,j)/(conj(RHPp47(i))+RHPp47(j));
        end
    end
    for i=1:nz47
        for j=1:np47
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz47(:,i))*inv(W1_zi)*W1_pj*Yp47(:,j)/(RHPz47(i)-RHPp47(j));
            
            W2_zi = evalfr(W2,RHPz47(i));
            W2_pj = evalfr(W2,RHPp47(j));
            Qzp2(i,j) = ctranspose(Yz47(:,i))*W2_zi*inv(W2_pj)*Yp47(:,j)/(RHPz47(i)-RHPp47(j));
        end
    end
    temp47 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSG47_min = sqrt(max(eig(temp47)));
end

%% Bounds on SGd471 (P_in,P_t,d1)
% Page 235, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms21;Gdms21]
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

    Qz1 = zeros(nz47,nz47); Qp1 = zeros(np47,np47); Qzp1 = zeros(nz47,np47); 
    Qz2 = zeros(nz47,nz47); Qp2 = zeros(np47,np47); Qzp2 = zeros(nz47,np47); 
    for i=1:nz47
        for j=1:nz47
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz47(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz47(:,j)/(RHPz47(i)+conj(RHPz47(j)));
            
            W2_zi = evalfr(W2,RHPz47(i));
            W2_zj = evalfr(W2,RHPz47(j));
            Qz2(i,j) = ctranspose(Yz47(:,i))*W2_zi*ctranspose(W2_zj)*Yz47(:,j)/(RHPz47(i)+conj(RHPz47(j)));
        end
    end
    for i=1:np47
        for j=1:np47
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp47(:,i))*ctranspose(W1)*W1_pj*Yp47(:,j)/(conj(RHPp47(i))+RHPp47(j));
            
            W2_pi = evalfr(W2,RHPp47(i));
            W2_pj = evalfr(W2,RHPp47(j));
            Qp2(i,j) = ctranspose(Yp47(:,i))*ctranspose(pinv(W2_pi))*pinv(W2_pj)*Yp47(:,j)/(conj(RHPp47(i))+RHPp47(j));
        end
    end
    for i=1:nz47
        for j=1:np47
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz47(:,i))*inv(W1_zi)*W1_pj*Yp47(:,j)/(RHPz47(i)-RHPp47(j));
            
            W2_zi = evalfr(W2,RHPz47(i));
            W2_pj = evalfr(W2,RHPp47(j));
            Qzp2(i,j) = ctranspose(Yz47(:,i))*W2_zi*pinv(W2_pj)*Yp47(:,j)/(RHPz47(i)-RHPp47(j));
        end
    end
    temp47 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd471_min = sqrt(max(eig(temp47)));
end

%% Bounds on SGd472 (P_in,P_t,d2)
% Page 235, Bound on ||W1SW2|| with W1 = 1; and W2 = [Gdms22;Gdms22]
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
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_zj = W1; %evalfr(W1,RHPz2(j));
            Qz1(i,j) = ctranspose(Yz47(:,i))*inv(W1_zi)*ctranspose(inv(W1_zj))*Yz47(:,j)/(RHPz47(i)+conj(RHPz47(j)));
            
            W2_zi = evalfr(W2,RHPz47(i));
            W2_zj = evalfr(W2,RHPz47(j));
            Qz2(i,j) = ctranspose(Yz47(:,i))*W2_zi*ctranspose(W2_zj)*Yz47(:,j)/(RHPz47(i)+conj(RHPz47(j)));
        end
    end
    for i=1:np47
        for j=1:np47
            W1_pi = W1; %evalfr(W1,RHPp2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qp1(i,j) = ctranspose(Yp47(:,i))*ctranspose(W1)*W1_pj*Yp47(:,j)/(conj(RHPp47(i))+RHPp47(j));
            
            W2_pi = evalfr(W2,RHPp47(i));
            W2_pj = evalfr(W2,RHPp47(j));
            Qp2(i,j) = ctranspose(Yp47(:,i))*ctranspose(pinv(W2_pi))*pinv(W2_pj)*Yp47(:,j)/(conj(RHPp47(i))+RHPp47(j));
        end
    end
    for i=1:nz47
        for j=1:np47
            W1_zi = W1; %evalfr(W1,RHPz2(i));
            W1_pj = W1; %evalfr(W1,RHPp2(j));
            Qzp1(i,j) = ctranspose(Yz47(:,i))*pinv(W1_zi)*W1_pj*Yp47(:,j)/(RHPz47(i)-RHPp47(j));
            
            W2_zi = evalfr(W2,RHPz47(i));
            W2_pj = evalfr(W2,RHPp47(j));
            Qzp2(i,j) = ctranspose(Yz47(:,i))*W2_zi*pinv(W2_pj)*Yp47(:,j)/(RHPz47(i)-RHPp47(j));
        end
    end
    temp47 = sqrtm(inv(Qz1))*(Qz2 + Qzp2*inv(Qp2)*ctranspose(Qzp2))*sqrtm(inv(Qz1));
    GammaSGd472_min = real(sqrt(max(eig(temp47))));
end

%% Pole Vectors
    np47=length(RHPp47); nz47=length(RHPz47);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp47 = zeros(l,np47);
    Yz47 = zeros(l,nz47);
    for i=1:np47
        Yp47(:,i) = Cs([4 7],:)*V(:,i); %Pole direction
    end
    Yp_max = max(max(abs(Yp47)),[],2);
%%
disp('Pole vectors');
disp(['Pole vectors for P_bh & P_rb: ', num2str(Yp_max)]);

disp('##########################################################');
disp('Bounds on S and T');
disp(['The lowest achivable peak for S14 and T14 is: ', num2str(Msmin14)]);
disp(['The lowest achivable peak for S15 and T15 is: ', num2str(Msmin15)]);
disp(['The lowest achivable peak for S17 and T17 is: ', num2str(Msmin17)]);
disp(['The lowest achivable peak for S19 and T19 is: ', num2str(Msmin19)]);
disp(['The lowest achivable peak for S25 and T25 is: ', num2str(Msmin25)]);
disp(['The lowest achivable peak for S27 and T27 is: ', num2str(Msmin27)]);
disp(['The lowest achivable peak for S29 and T29 is: ', num2str(Msmin29)]);
disp(['The lowest achivable peak for S47 and T47 is: ', num2str(Msmin47)]);
disp('##########################################################');
disp('Bounds on KS');
disp(['The lowest achivable peak for KS14 is: ', num2str(KSmin14)]);
disp(['The lowest achivable peak for KS15 is: ', num2str(KSmin15)]);
disp(['The lowest achivable peak for KS17 is: ', num2str(KSmin17)]);
disp(['The lowest achivable peak for KS19 is: ', num2str(KSmin19)]);
disp(['The lowest achivable peak for KS25 is: ', num2str(KSmin25)]);
disp(['The lowest achivable peak for KS27 is: ', num2str(KSmin27)]);
disp(['The lowest achivable peak for KS29 is: ', num2str(KSmin29)]);
disp(['The lowest achivable peak for KS47 is: ', num2str(KSmin47)]);
disp('##########################################################');
disp('Bounds on KSGd form equation 6.26 (tight)');
disp(['The lowest achivable peak for KSGd141 (P_bh,P_rb,d1) is: ', num2str(KSGdmin141)]);
disp(['The lowest achivable peak for KSGd151 (P_bh,P_rb,d1) is: ', num2str(KSGdmin151)]);
disp(['The lowest achivable peak for KSGd171 (P_bh,P_rb,d1) is: ', num2str(KSGdmin171)]);
disp(['The lowest achivable peak for KSGd191 (P_bh,P_rb,d1) is: ', num2str(KSGdmin191)]);
disp(['The lowest achivable peak for KSGd251 (P_bh,P_rb,d1) is: ', num2str(KSGdmin251)]);
disp(['The lowest achivable peak for KSGd271 (P_bh,P_rb,d1) is: ', num2str(KSGdmin271)]);
disp(['The lowest achivable peak for KSGd291 (P_bh,P_rb,d1) is: ', num2str(KSGdmin291)]);
disp(['The lowest achivable peak for KSGd471 (P_bh,P_rb,d1) is: ', num2str(KSGdmin471)]);
disp(['The lowest achivable peak for KSGd142 (P_bh,P_rb,d2) is: ', num2str(KSGdmin142)]);
disp(['The lowest achivable peak for KSGd152 (P_bh,P_rb,d2) is: ', num2str(KSGdmin152)]);
disp(['The lowest achivable peak for KSGd172 (P_bh,P_rb,d2) is: ', num2str(KSGdmin172)]);
disp(['The lowest achivable peak for KSGd192 (P_bh,P_rb,d2) is: ', num2str(KSGdmin192)]);
disp(['The lowest achivable peak for KSGd252 (P_bh,P_rb,d2) is: ', num2str(KSGdmin252)]);
disp(['The lowest achivable peak for KSGd272 (P_bh,P_rb,d2) is: ', num2str(KSGdmin272)]);
disp(['The lowest achivable peak for KSGd292 (P_bh,P_rb,d2) is: ', num2str(KSGdmin292)]);
disp(['The lowest achivable peak for KSGd472 (P_bh,P_rb,d2) is: ', num2str(KSGdmin472)]);
disp('##########################################################');
disp('Bounds on KSGd form equation 6.27 (not tight)');
disp(['The lowest achivable peak for KSGd141 (P_bh,P_rb,d1) is: ', num2str(KSGd141_min)]);
disp(['The lowest achivable peak for KSGd151 (P_bh,P_rb,d1) is: ', num2str(KSGd151_min)]);
disp(['The lowest achivable peak for KSGd171 (P_bh,P_rb,d1) is: ', num2str(KSGd171_min)]);
disp(['The lowest achivable peak for KSGd191 (P_bh,P_rb,d1) is: ', num2str(KSGd191_min)]);
disp(['The lowest achivable peak for KSGd251 (P_bh,P_rb,d1) is: ', num2str(KSGd251_min)]);
disp(['The lowest achivable peak for KSGd271 (P_bh,P_rb,d1) is: ', num2str(KSGd271_min)]);
disp(['The lowest achivable peak for KSGd291 (P_bh,P_rb,d1) is: ', num2str(KSGd291_min)]);
disp(['The lowest achivable peak for KSGd471 (P_bh,P_rb,d1) is: ', num2str(KSGd471_min)]);
disp(['The lowest achivable peak for KSGd142 (P_bh,P_rb,d2) is: ', num2str(KSGd142_min)]);
disp(['The lowest achivable peak for KSGd152 (P_bh,P_rb,d2) is: ', num2str(KSGd152_min)]);
disp(['The lowest achivable peak for KSGd172 (P_bh,P_rb,d2) is: ', num2str(KSGd172_min)]);
disp(['The lowest achivable peak for KSGd192 (P_bh,P_rb,d2) is: ', num2str(KSGd192_min)]);
disp(['The lowest achivable peak for KSGd252 (P_bh,P_rb,d2) is: ', num2str(KSGd252_min)]);
disp(['The lowest achivable peak for KSGd272 (P_bh,P_rb,d2) is: ', num2str(KSGd272_min)]);
disp(['The lowest achivable peak for KSGd292 (P_bh,P_rb,d2) is: ', num2str(KSGd292_min)]);
disp(['The lowest achivable peak for KSGd472 (P_bh,P_rb,d2) is: ', num2str(KSGd472_min)]);
disp('##########################################################');
disp('Bounds on SG');
disp(['The lowest achivable peak for SG14 (P_bh,P_rb) is: ', num2str(GammaSG14_min)]);
disp(['The lowest achivable peak for SG15 (P_bh,P_rb) is: ', num2str(GammaSG15_min)]);
disp(['The lowest achivable peak for SG17 (P_bh,P_rb) is: ', num2str(GammaSG17_min)]);
disp(['The lowest achivable peak for SG19 (P_bh,P_rb) is: ', num2str(GammaSG19_min)]);
disp(['The lowest achivable peak for SG25 (P_bh,P_rb) is: ', num2str(GammaSG25_min)]);
disp(['The lowest achivable peak for SG29 (P_bh,P_rb) is: ', num2str(GammaSG29_min)]);
disp(['The lowest achivable peak for SG27 (P_bh,P_rb) is: ', num2str(GammaSG27_min)]);
disp(['The lowest achivable peak for SG47 (P_bh,P_rb) is: ', num2str(GammaSG47_min)]);
disp('##########################################################');
disp('Bounds on SGd');
disp(['The lowest achivable peak for SGd141 (P_bh,P_rb,d1) is: ',num2str(GammaSGd141_min )]);
disp(['The lowest achivable peak for SGd151 (P_bh,P_rb,d1) is: ',num2str(GammaSGd151_min )]);
disp(['The lowest achivable peak for SGd171 (P_bh,P_rb,d1) is: ',num2str(GammaSGd171_min )]);
disp(['The lowest achivable peak for SGd191 (P_bh,P_rb,d1) is: ',num2str(GammaSGd191_min )]);
disp(['The lowest achivable peak for SGd251 (P_bh,P_rb,d1) is: ',num2str(GammaSGd251_min )]);
disp(['The lowest achivable peak for SGd271 (P_bh,P_rb,d1) is: ',num2str(GammaSGd271_min )]);
disp(['The lowest achivable peak for SGd291 (P_bh,P_rb,d1) is: ',num2str(GammaSGd291_min )]);
disp(['The lowest achivable peak for SGd471 (P_bh,P_rb,d1) is: ',num2str(GammaSGd471_min )]);
disp(['The lowest achivable peak for SGd142 (P_bh,P_rb,d2) is: ',num2str(GammaSGd142_min )]);
disp(['The lowest achivable peak for SGd152 (P_bh,P_rb,d2) is: ',num2str(GammaSGd152_min )]);
disp(['The lowest achivable peak for SGd172 (P_bh,P_rb,d2) is: ',num2str(GammaSGd172_min )]);
disp(['The lowest achivable peak for SGd192 (P_bh,P_rb,d2) is: ',num2str(GammaSGd192_min )]);
disp(['The lowest achivable peak for SGd252 (P_bh,P_rb,d2) is: ',num2str(GammaSGd252_min )]);
disp(['The lowest achivable peak for SGd272 (P_bh,P_rb,d2) is: ',num2str(GammaSGd272_min )]);
disp(['The lowest achivable peak for SGd292 (P_bh,P_rb,d2) is: ',num2str(GammaSGd292_min )]);
disp(['The lowest achivable peak for SGd472 (P_bh,P_rb,d2) is: ',num2str(GammaSGd472_min )]);