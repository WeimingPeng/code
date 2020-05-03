clc
clear all
warning ('off')
Current_Folder = cd('../WPR_model');
%% **** System Inputs ****
par = v1_new_6d_parametersSS(); %All parameters defined in one file
z1 = 0.3;
z2 = 0.2;
global u0;
u0 = [z1;z2;0;0];
%% ******** Initializing *************
disp('Initializing...')
tic
[x0,y0,par] =  v1_new_6d_initialize(u0,par);
toc
disp('Initialization completed.');
disp('#########################');
%% ******** Linearization *************
[A,B,C,D]=v2_new_6d_linmod(x0,y0,u0,par); %analytical
%[A,B,C,D] = v1_new_6d_linmodl_num(x0,u0,par); %numerical
sys =ss(A,B,C,D);
set(sys,'inputname',{'Z1','Z2','d1','d2'},'outputname',{'P_bh','P_wh','W_in','P_in','P_rb','DP_r','P_t','Q_out','W_out','Rho_t','Alpha_L'},'statename',{'m_Gw','m_Lw','m_Gp','m_Lp','m_Gr','m_Lr'});
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
set(sys,'inputname',{'Z1','Z2','d1','d2'},'outputname',{'P_bh','P_wh','W_in','P_in','P_rb','DP_r','P_t','Q_out','W_out','Rho_t','Alpha_L'},'statename',{'m_Gw','m_Lw','m_Gp','m_Lp','m_Gr','m_Lr'});
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
    GammaSG15_min = sqrt(max(eig(temp15)));
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
    GammaSGd151_min = sqrt(max(eig(temp15)));
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

%%
disp('Pole vectors');
disp(['Pole vectors for P_bh & P_rb: ', num2str(Yp_max)]);

disp('##########################################################');
disp('Bounds on S and T');
disp(['The lowest achivable peak for S15 and T15 is: ', num2str(Msmin15)]);
disp(['The lowest achivable peak for S17 and T17 is: ', num2str(Msmin17)]);
disp(['The lowest achivable peak for S25 and T25 is: ', num2str(Msmin25)]);
disp(['The lowest achivable peak for S27 and T27 is: ', num2str(Msmin27)]);
disp(['The lowest achivable peak for S47 and T47 is: ', num2str(Msmin47)]);
disp('##########################################################');
disp('Bounds on KS');
disp(['The lowest achivable peak for KS15 is: ', num2str(KSmin15)]);
disp(['The lowest achivable peak for KS17 is: ', num2str(KSmin17)]);
disp(['The lowest achivable peak for KS25 is: ', num2str(KSmin25)]);
disp(['The lowest achivable peak for KS27 is: ', num2str(KSmin27)]);
disp(['The lowest achivable peak for KS47 is: ', num2str(KSmin47)]);
disp('##########################################################');
disp('Bounds on KSGd form equation 6.26 (tight)');
disp(['The lowest achivable peak for KSGd151 (P_bh,P_rb,d1) is: ', num2str(KSGdmin151)]);
disp(['The lowest achivable peak for KSGd171 (P_bh,P_rb,d1) is: ', num2str(KSGdmin171)]);
disp(['The lowest achivable peak for KSGd251 (P_bh,P_rb,d1) is: ', num2str(KSGdmin251)]);
disp(['The lowest achivable peak for KSGd271 (P_bh,P_rb,d1) is: ', num2str(KSGdmin271)]);
disp(['The lowest achivable peak for KSGd471 (P_bh,P_rb,d1) is: ', num2str(KSGdmin471)]);
disp(['The lowest achivable peak for KSGd152 (P_bh,P_rb,d2) is: ', num2str(KSGdmin152)]);
disp(['The lowest achivable peak for KSGd172 (P_bh,P_rb,d2) is: ', num2str(KSGdmin172)]);
disp(['The lowest achivable peak for KSGd252 (P_bh,P_rb,d2) is: ', num2str(KSGdmin252)]);
disp(['The lowest achivable peak for KSGd272 (P_bh,P_rb,d2) is: ', num2str(KSGdmin272)]);
disp(['The lowest achivable peak for KSGd472 (P_bh,P_rb,d2) is: ', num2str(KSGdmin472)]);
disp('##########################################################');
disp('Bounds on KSGd form equation 6.27 (not tight)');
disp(['The lowest achivable peak for KSGd151 (P_bh,P_rb,d1) is: ', num2str(KSGd151_min)]);
disp(['The lowest achivable peak for KSGd171 (P_bh,P_rb,d1) is: ', num2str(KSGd171_min)]);
disp(['The lowest achivable peak for KSGd251 (P_bh,P_rb,d1) is: ', num2str(KSGd251_min)]);
disp(['The lowest achivable peak for KSGd271 (P_bh,P_rb,d1) is: ', num2str(KSGd271_min)]);
disp(['The lowest achivable peak for KSGd471 (P_bh,P_rb,d1) is: ', num2str(KSGd471_min)]);
disp(['The lowest achivable peak for KSGd152 (P_bh,P_rb,d2) is: ', num2str(KSGd152_min)]);
disp(['The lowest achivable peak for KSGd172 (P_bh,P_rb,d2) is: ', num2str(KSGd172_min)]);
disp(['The lowest achivable peak for KSGd252 (P_bh,P_rb,d2) is: ', num2str(KSGd252_min)]);
disp(['The lowest achivable peak for KSGd272 (P_bh,P_rb,d2) is: ', num2str(KSGd272_min)]);
disp(['The lowest achivable peak for KSGd472 (P_bh,P_rb,d2) is: ', num2str(KSGd472_min)]);
disp('##########################################################');
disp('Bounds on SG');
disp(['The lowest achivable peak for SG15 (P_bh,P_rb) is: ', num2str(GammaSG15_min)]);
disp(['The lowest achivable peak for SG17 (P_bh,P_rb) is: ', num2str(GammaSG17_min)]);
disp(['The lowest achivable peak for SG25 (P_bh,P_rb) is: ', num2str(GammaSG25_min)]);
disp(['The lowest achivable peak for SG27 (P_bh,P_rb) is: ', num2str(GammaSG27_min)]);
disp(['The lowest achivable peak for SG47 (P_bh,P_rb) is: ', num2str(GammaSG47_min)]);
disp('##########################################################');
disp('Bounds on SGd');
disp(['The lowest achivable peak for SGd151 (P_bh,P_rb,d1) is: ',num2str(GammaSGd151_min )]);
disp(['The lowest achivable peak for SGd171 (P_bh,P_rb,d1) is: ',num2str(GammaSGd171_min )]);
disp(['The lowest achivable peak for SGd251 (P_bh,P_rb,d1) is: ',num2str(GammaSGd251_min )]);
disp(['The lowest achivable peak for SGd271 (P_bh,P_rb,d1) is: ',num2str(GammaSGd271_min )]);
disp(['The lowest achivable peak for SGd471 (P_bh,P_rb,d1) is: ',num2str(GammaSGd471_min )]);
disp(['The lowest achivable peak for SGd152 (P_bh,P_rb,d2) is: ',num2str(GammaSGd152_min )]);
disp(['The lowest achivable peak for SGd172 (P_bh,P_rb,d2) is: ',num2str(GammaSGd172_min )]);
disp(['The lowest achivable peak for SGd252 (P_bh,P_rb,d2) is: ',num2str(GammaSGd252_min )]);
disp(['The lowest achivable peak for SGd272 (P_bh,P_rb,d2) is: ',num2str(GammaSGd272_min )]);
disp(['The lowest achivable peak for SGd472 (P_bh,P_rb,d2) is: ',num2str(GammaSGd472_min )]);