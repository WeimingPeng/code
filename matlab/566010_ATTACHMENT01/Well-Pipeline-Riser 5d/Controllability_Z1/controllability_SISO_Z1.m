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
% ******** Controllability analysis*************
%*** Sensitivity and complementery sensitivity peak***
%% Sensitivity peak for P_bh.
sys1 = sys('P_bh','Z1');
[p, z1]=pzmap(sys1);
RHPp1=p(find(p>0));
RHPz1=z1(find(z1>0));
l=1;
if (isempty(RHPp1) || isempty(RHPz1))
    Msmin1=1;
else
    np1=length(RHPp1); nz1=length(RHPz1);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp1 = zeros(l,np1);
    Yz1 = zeros(l,nz1);
    for i=1:np1
        Yp1(:,i) = Cs(1,:)*V(:,i)/norm(Cs(1,:)*V(:,i)); %Pole direction
    end
    for j=1:nz1
        [U,S,V]=svd(evalfr(G1,RHPz1(j))); Yz1(:,j)=U(:,end); %zero direction
    end
    Qz1 = zeros(nz1,nz1); Qp1 = zeros(np1,np1); Qzp1 = zeros(nz1,np1);
    for i=1:nz1
        for j=1:nz1
            Qz1(i,j) = ctranspose(Yz1(:,i))*Yz1(:,j)/(RHPz1(i)+conj(RHPz1(j)));
        end
    end
    for i=1:np1
        for j=1:np1
            Qp1(i,j) = ctranspose(Yp1(:,i))*Yp1(:,j)/(conj(RHPp1(i))+RHPp1(j));
        end
    end
    for i=1:nz1
        for j=1:np1
            Qzp1(i,j) = ctranspose(Yz1(:,i))*Yp1(:,j)/(RHPz1(i)-RHPp1(j));
        end
    end
    %     Qp1=(Yp'*Yp).*(1./(diag(RHPp1')*ones(np1)+ones(np1)*diag(RHPp1)));
    %     Qz1=(Yz'*Yz).*(1./(diag(RHPz1)*ones(nz1)+ones(nz1)*diag(RHPz1')));
    %     Qzp1=(Yz'*Yp).*(1./(diag(RHPz1)*ones(nz1,np1)-ones(nz1,np1)*diag(RHPp1)));
    Msmin1=sqrt(1+norm(sqrtm(inv(Qz1))*Qzp1*sqrtm(inv(Qp1)))^2);
end

%% Sensitivity peak for P_wh
sys2 = sys('P_wh','Z1');
l=1;
[p, z2]=pzmap(sys2);
RHPp2=p(find(p>0));
RHPz2=z2(find(z2>0));

if (isempty(RHPp2) || isempty(RHPz2))
    Msmin2=1;
else
    np2=length(RHPp2); nz2=length(RHPz2);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp2 = zeros(l,np2);
    Yz2 = zeros(l,nz2);
    for i=1:np2
        Yp2(:,i) = Cs(2,:)*V(:,i)/norm(Cs(2,:)*V(:,i)); %Pole direction
    end
    for j=1:nz2
        [U,S,V]=svd(evalfr(G2,RHPz2(j))); Yz2(:,j)=U(:,end); %zero direction
    end
    Qz2 = zeros(nz2,nz2); Qp2 = zeros(np2,np2); Qzp2 = zeros(nz2,np2);
    for i=1:nz2
        for j=1:nz2
            Qz2(i,j) = ctranspose(Yz2(:,i))*Yz2(:,j)/(RHPz2(i)+conj(RHPz2(j)));
        end
    end
    for i=1:np2
        for j=1:np2
            Qp2(i,j) = ctranspose(Yp2(:,i))*Yp2(:,j)/(conj(RHPp2(i))+RHPp2(j));
        end
    end
    for i=1:nz2
        for j=1:np2
            Qzp2(i,j) = ctranspose(Yz2(:,i))*Yp2(:,j)/(RHPz2(i)-RHPp2(j));
        end
    end
    
    %     Qp2=(Yp2'*Yp2).*(1./(diag(RHPp2')*ones(np2)+ones(np2)*diag(RHPp2)));
    %     Qz2=(Yz2'*Yz2).*(1./(diag(RHPz2)*ones(nz2)+ones(nz2)*diag(RHPz2')));
    %     Qzp2=(Yz2'*Yp2).*(1./(diag(RHPz2)*ones(nz2,np2)-ones(nz2,np2)*diag(RHPp2)));
    Msmin2=sqrt(1+norm((Qz2^(-0.5))*Qzp2*(Qp2^(-0.5)))^2);
end


%% Sensitivity peak for W_in
sys3 = sys('W_in','Z1');
l=1;
[p, z3]=pzmap(sys3);
RHPp3=p(find(p>0));
RHPz3=z3(find(z3>0));

if (isempty(RHPp3) || isempty(RHPz3))
    Msmin3=1;
else
    np3=length(RHPp3); nz3=length(RHPz3);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp3 = zeros(l,np3);
    Yz3 = zeros(l,nz3);
    for i=1:np3
        Yp3(:,i) = Cs(3,:)*V(:,i)/norm(Cs(3,:)*V(:,i)); %Pole direction
    end
    for j=1:nz3
        [U,S,V]=svd(evalfr(G3,RHPz3(j))); Yz3(:,j)=U(:,end); %zero direction
    end
    
    Qz3 = zeros(nz3,nz3); Qp3 = zeros(np3,np3); Qzp3 = zeros(nz3,np3);
    for i=1:nz3
        for j=1:nz3
            Qz3(i,j) = ctranspose(Yz3(:,i))*Yz3(:,j)/(RHPz3(i)+conj(RHPz3(j)));
        end
    end
    for i=1:np3
        for j=1:np3
            Qp3(i,j) = ctranspose(Yp3(:,i))*Yp3(:,j)/(conj(RHPp3(i))+RHPp3(j));
        end
    end
    for i=1:nz3
        for j=1:np3
            Qzp3(i,j) = ctranspose(Yz3(:,i))*Yp3(:,j)/(RHPz3(i)-RHPp3(j));
        end
    end
    
    %     Qp3=(Yp3'*Yp3).*(1./(diag(RHPp3')*ones(np3)+ones(np3)*diag(RHPp3)));
    %     Qz3=(Yz3'*Yz3).*(1./(diag(RHPz3)*ones(nz3)+ones(nz3)*diag(RHPz3')));
    %     Qzp3=(Yz3'*Yp3).*(1./(diag(RHPz3)*ones(nz3,np3)-ones(nz3,np3)*diag(RHPp3)));
    Msmin3=sqrt(1+norm((Qz3^(-0.5))*Qzp3*(Qp3^(-0.5)))^2);
end



%% Sensitivity peak for P_in
sys4 = sys('P_in','Z1');
l=1;
[p, z4]=pzmap(sys4);
RHPp4=p(find(p>0));
RHPz4=z4(find(z4>0));

if (isempty(RHPp4) || isempty(RHPz4))
    Msmin4=1;
else
    np4=length(RHPp4); nz4=length(RHPz4);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp4 = zeros(l,np4);
    Yz4 = zeros(l,nz4);
    for i=1:np4
        Yp4(:,i) = Cs(4,:)*V(:,i)/norm(Cs(4,:)*V(:,i)); %Pole direction
    end
    for j=1:nz4
        [U,S,V]=svd(evalfr(G4,RHPz4(j))); Yz4(:,j)=U(:,end); %zero direction
    end
    
    Qz4 = zeros(nz4,nz4); Qp4 = zeros(np4,np4); Qzp4 = zeros(nz4,np4);
    for i=1:nz4
        for j=1:nz4
            Qz4(i,j) = ctranspose(Yz4(:,i))*Yz4(:,j)/(RHPz4(i)+conj(RHPz4(j)));
        end
    end
    for i=1:np4
        for j=1:np4
            Qp4(i,j) = ctranspose(Yp4(:,i))*Yp4(:,j)/(conj(RHPp4(i))+RHPp4(j));
        end
    end
    for i=1:nz4
        for j=1:np4
            Qzp4(i,j) = ctranspose(Yz4(:,i))*Yp4(:,j)/(RHPz4(i)-RHPp4(j));
        end
    end
    
    %     Qp4=(Yp4'*Yp4).*(1./(diag(RHPp4')*ones(np4)+ones(np4)*diag(RHPp4)));
    %     Qz4=(Yz4'*Yz4).*(1./(diag(RHPz4)*ones(nz4)+ones(nz4)*diag(RHPz4')));
    %     Qzp4=(Yz4'*Yp4).*(1./(diag(RHPz4)*ones(nz4,np4)-ones(nz4,np4)*diag(RHPp4)));
    Msmin4=sqrt(1+norm((Qz4^(-0.5))*Qzp4*(Qp4^(-0.5)))^2);
end

%% Sensitivity peak for P_rb.
sys5 = sys('P_rb','Z1');
[p, z5]=pzmap(sys5);
RHPp5=p(find(p>0));
RHPz5=z5(find(z5>0));

if (isempty(RHPp5) || isempty(RHPz5))
    Msmin5=1;
else
    np5=length(RHPp5); nz5=length(RHPz5);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp5 = zeros(l,np5);
    Yz5 = zeros(l,nz5);
    for i=1:np5
        Yp5(:,i) = Cs(5,:)*V(:,i)/norm(Cs(5,:)*V(:,i)); %Pole direction
    end
    for j=1:nz5
        [U,S,V]=svd(evalfr(G5,RHPz5(j))); Yz5(:,j)=U(:,end); %zero direction
    end
    Qz5 = zeros(nz5,nz5); Qp5 = zeros(np5,np5); Qzp5 = zeros(nz5,np5);
    for i=1:nz5
        for j=1:nz5
            Qz5(i,j) = ctranspose(Yz5(:,i))*Yz5(:,j)/(RHPz5(i)+conj(RHPz5(j)));
        end
    end
    for i=1:np5
        for j=1:np5
            Qp5(i,j) = ctranspose(Yp5(:,i))*Yp5(:,j)/(conj(RHPp5(i))+RHPp5(j));
        end
    end
    for i=1:nz5
        for j=1:np5
            Qzp5(i,j) = ctranspose(Yz5(:,i))*Yp5(:,j)/(RHPz5(i)-RHPp5(j));
        end
    end
    %     Qp5=(Yp'*Yp).*(1./(diag(RHPp5')*ones(np5)+ones(np5)*diag(RHPp5)));
    %     Qz5=(Yz'*Yz).*(1./(diag(RHPz5)*ones(nz5)+ones(nz5)*diag(RHPz5')));
    %     Qzp5=(Yz'*Yp).*(1./(diag(RHPz5)*ones(nz5,np5)-ones(nz5,np5)*diag(RHPp5)));
    Msmin5=sqrt(1+norm(sqrtm(inv(Qz5))*Qzp5*sqrtm(inv(Qp5)))^2);
end


%% Sensitivity peak for DPr.
sys6 = sys('DP_r','Z1');
[p, z6]=pzmap(sys6);
RHPp6=p(find(p>0));
RHPz6=z6(find(z6>0));

if (isempty(RHPp6) || isempty(RHPz6))
    Msmin6=1;
else
    np6=length(RHPp6); nz6=length(RHPz6);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp6 = zeros(l,np6);
    Yz6 = zeros(l,nz6);
    for i=1:np6
        Yp6(:,i) = Cs(6,:)*V(:,i)/norm(Cs(6,:)*V(:,i)); %Pole direction
    end
    for j=1:nz6
        [U,S,V]=svd(evalfr(sys6,RHPz6(j))); Yz6(:,j)=U(:,end); %zero direction
    end
    Qz6 = zeros(nz6,nz6); Qp6 = zeros(np6,np6); Qzp6 = zeros(nz6,np6);
    for i=1:nz6
        for j=1:nz6
            Qz6(i,j) = ctranspose(Yz6(:,i))*Yz6(:,j)/(RHPz6(i)+conj(RHPz6(j)));
        end
    end
    for i=1:np6
        for j=1:np6
            Qp6(i,j) = ctranspose(Yp6(:,i))*Yp6(:,j)/(conj(RHPp6(i))+RHPp6(j));
        end
    end
    for i=1:nz6
        for j=1:np6
            Qzp6(i,j) = ctranspose(Yz6(:,i))*Yp6(:,j)/(RHPz6(i)-RHPp6(j));
        end
    end
    %     Qp6=(Yp'*Yp).*(1./(diag(RHPp6')*ones(np6)+ones(np6)*diag(RHPp6)));
    %     Qz6=(Yz'*Yz).*(1./(diag(RHPz6)*ones(nz6)+ones(nz6)*diag(RHPz6')));
    %     Qzp6=(Yz'*Yp).*(1./(diag(RHPz6)*ones(nz6,np6)-ones(nz6,np6)*diag(RHPp6)));
    Msmin6=sqrt(1+norm(sqrtm(inv(Qz6))*Qzp6*sqrtm(inv(Qp6)))^2);
end


%% Sensitivity peak for P_t.
sys7 = sys('P_t','Z1');
[p, z7]=pzmap(sys7);
RHPp7=p(find(p>0));
RHPz7=z7(find(z7>0));

if (isempty(RHPp7) || isempty(RHPz7))
    Msmin7=1;
else
    np7=length(RHPp7); nz7=length(RHPz7);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp7 = zeros(l,np7);
    Yz7 = zeros(l,nz7);
    for i=1:np7
        Yp7(:,i) = Cs(7,:)*V(:,i)/norm(Cs(7,:)*V(:,i)); %Pole direction
    end
    for j=1:nz7
        [U,S,V]=svd(evalfr(sys7,RHPz7(j))); Yz7(:,j)=U(:,end); %zero direction
    end
    Qz7 = zeros(nz7,nz7); Qp7 = zeros(np7,np7); Qzp7 = zeros(nz7,np7);
    for i=1:nz7
        for j=1:nz7
            Qz7(i,j) = ctranspose(Yz7(:,i))*Yz7(:,j)/(RHPz7(i)+conj(RHPz7(j)));
        end
    end
    for i=1:np7
        for j=1:np7
            Qp7(i,j) = ctranspose(Yp7(:,i))*Yp7(:,j)/(conj(RHPp7(i))+RHPp7(j));
        end
    end
    for i=1:nz7
        for j=1:np7
            Qzp7(i,j) = ctranspose(Yz7(:,i))*Yp7(:,j)/(RHPz7(i)-RHPp7(j));
        end
    end
    %     Qp7=(Yp'*Yp).*(1./(diag(RHPp7')*ones(np7)+ones(np7)*diag(RHPp7)));
    %     Qz7=(Yz'*Yz).*(1./(diag(RHPz7)*ones(nz7)+ones(nz7)*diag(RHPz7')));
    %     Qzp7=(Yz'*Yp).*(1./(diag(RHPz7)*ones(nz7,np7)-ones(nz7,np7)*diag(RHPp7)));
    Msmin7=sqrt(1+norm(sqrtm(inv(Qz7))*Qzp7*sqrtm(inv(Qp7)))^2);
end

%% Sensitivity peak for Q_out.
sys8 = sys('Q_out','Z1');
[p, z8]=pzmap(sys8);
RHPp8=p(find(p>0));
RHPz8=z8(find(z8>0));

if (isempty(RHPp8) || isempty(RHPz8))
    Msmin8=1;
else
    np8=length(RHPp8); nz8=length(RHPz8);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp8 = zeros(l,np8);
    Yz8 = zeros(l,nz8);
    for i=1:np8
        Yp8(:,i) = Cs(8,:)*V(:,i)/norm(Cs(8,:)*V(:,i)); %Pole direction
    end
    for j=1:nz8
        [U,S,V]=svd(evalfr(sys8,RHPz8(j))); Yz8(:,j)=U(:,end); %zero direction
    end
    Qz8 = zeros(nz8,nz8); Qp8 = zeros(np8,np8); Qzp8 = zeros(nz8,np8);
    for i=1:nz8
        for j=1:nz8
            Qz8(i,j) = ctranspose(Yz8(:,i))*Yz8(:,j)/(RHPz8(i)+conj(RHPz8(j)));
        end
    end
    for i=1:np8
        for j=1:np8
            Qp8(i,j) = ctranspose(Yp8(:,i))*Yp8(:,j)/(conj(RHPp8(i))+RHPp8(j));
        end
    end
    for i=1:nz8
        for j=1:np8
            Qzp8(i,j) = ctranspose(Yz8(:,i))*Yp8(:,j)/(RHPz8(i)-RHPp8(j));
        end
    end
    %     Qp8=(Yp'*Yp).*(1./(diag(RHPp8')*ones(np8)+ones(np8)*diag(RHPp8)));
    %     Qz8=(Yz'*Yz).*(1./(diag(RHPz8)*ones(nz8)+ones(nz8)*diag(RHPz8')));
    %     Qzp8=(Yz'*Yp).*(1./(diag(RHPz8)*ones(nz8,np8)-ones(nz8,np8)*diag(RHPp8)));
    Msmin8=sqrt(1+norm(sqrtm(inv(Qz8))*Qzp8*sqrtm(inv(Qp8)))^2);
end
%% Sensitivity peak for W_out.
sys9 = sys('W_out','Z1');
[p, z9]=pzmap(sys9);
RHPp9=p(find(p>0));
RHPz9=z9(find(z9>0));

if (isempty(RHPp9) || isempty(RHPz9))
    Msmin9=1;
else
    np9=length(RHPp9); nz9=length(RHPz9);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp9 = zeros(l,np9);
    Yz9 = zeros(l,nz9);
    for i=1:np9
        Yp9(:,i) = Cs(9,:)*V(:,i)/norm(Cs(9,:)*V(:,i)); %Pole direction
    end
    for j=1:nz9
        [U,S,V]=svd(evalfr(sys9,RHPz9(j))); Yz9(:,j)=U(:,end); %zero direction
    end
    Qz9 = zeros(nz9,nz9); Qp9 = zeros(np9,np9); Qzp9 = zeros(nz9,np9);
    for i=1:nz9
        for j=1:nz9
            Qz9(i,j) = ctranspose(Yz9(:,i))*Yz9(:,j)/(RHPz9(i)+conj(RHPz9(j)));
        end
    end
    for i=1:np9
        for j=1:np9
            Qp9(i,j) = ctranspose(Yp9(:,i))*Yp9(:,j)/(conj(RHPp9(i))+RHPp9(j));
        end
    end
    for i=1:nz9
        for j=1:np9
            Qzp9(i,j) = ctranspose(Yz9(:,i))*Yp9(:,j)/(RHPz9(i)-RHPp9(j));
        end
    end
    %     Qp9=(Yp'*Yp).*(1./(diag(RHPp9')*ones(np9)+ones(np9)*diag(RHPp9)));
    %     Qz9=(Yz'*Yz).*(1./(diag(RHPz9)*ones(nz9)+ones(nz9)*diag(RHPz9')));
    %     Qzp9=(Yz'*Yp).*(1./(diag(RHPz9)*ones(nz9,np9)-ones(nz9,np9)*diag(RHPp9)));
    Msmin9=sqrt(1+norm(sqrtm(inv(Qz9))*Qzp9*sqrtm(inv(Qp9)))^2);
end
%% Sensitivity peak for Rho_t.
sys10 = sys('Rho_t','Z1');
[p, z10]=pzmap(sys10);
RHPp10=p(find(p>0));
RHPz10=z10(find(z10>0));

if (isempty(RHPp10) || isempty(RHPz10))
    Msmin10=1;
else
    np10=length(RHPp10); nz10=length(RHPz10);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp10 = zeros(l,np10);
    Yz10 = zeros(l,nz10);
    for i=1:np10
        Yp10(:,i) = Cs(10,:)*V(:,i)/norm(Cs(10,:)*V(:,i)); %Pole direction
    end
    for j=1:nz10
        [U,S,V]=svd(evalfr(sys10,RHPz10(j))); Yz10(:,j)=U(:,end); %zero direction
    end
    Qz10 = zeros(nz10,nz10); Qp10 = zeros(np10,np10); Qzp10 = zeros(nz10,np10);
    for i=1:nz10
        for j=1:nz10
            Qz10(i,j) = ctranspose(Yz10(:,i))*Yz10(:,j)/(RHPz10(i)+conj(RHPz10(j)));
        end
    end
    for i=1:np10
        for j=1:np10
            Qp10(i,j) = ctranspose(Yp10(:,i))*Yp10(:,j)/(conj(RHPp10(i))+RHPp10(j));
        end
    end
    for i=1:nz10
        for j=1:np10
            Qzp10(i,j) = ctranspose(Yz10(:,i))*Yp10(:,j)/(RHPz10(i)-RHPp10(j));
        end
    end
    %     Qp10=(Yp'*Yp).*(1./(diag(RHPp10')*ones(np10)+ones(np10)*diag(RHPp10)));
    %     Qz10=(Yz'*Yz).*(1./(diag(RHPz10)*ones(nz10)+ones(nz10)*diag(RHPz10')));
    %     Qzp10=(Yz'*Yp).*(1./(diag(RHPz10)*ones(nz10,np10)-ones(nz10,np10)*diag(RHPp10)));
    Msmin10=sqrt(1+norm(sqrtm(inv(Qz10))*Qzp10*sqrtm(inv(Qp10)))^2);
end

%% Sensitivity peak for Alpha_L.
sys11 = sys('Alpha_L','Z1');
[p, z11]=pzmap(sys11);
RHPp11=p(find(p>0));
RHPz11=z11(find(z11>0));

if (isempty(RHPp11) || isempty(RHPz11))
    Msmin11=1;
else
    np11=length(RHPp11); nz11=length(RHPz11);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp11 = zeros(l,np11);
    Yz11 = zeros(l,nz11);
    for i=1:np11
        Yp11(:,i) = Cs(11,:)*V(:,i)/norm(Cs(11,:)*V(:,i)); %Pole direction
    end
    for j=1:nz11
        [U,S,V]=svd(evalfr(sys11,RHPz11(j))); Yz11(:,j)=U(:,end); %zero direction
    end
    Qz11 = zeros(nz11,nz11); Qp11 = zeros(np11,np11); Qzp11 = zeros(nz11,np11);
    for i=1:nz11
        for j=1:nz11
            Qz11(i,j) = ctranspose(Yz11(:,i))*Yz11(:,j)/(RHPz11(i)+conj(RHPz11(j)));
        end
    end
    for i=1:np11
        for j=1:np11
            Qp11(i,j) = ctranspose(Yp11(:,i))*Yp11(:,j)/(conj(RHPp11(i))+RHPp11(j));
        end
    end
    for i=1:nz11
        for j=1:np11
            Qzp11(i,j) = ctranspose(Yz11(:,i))*Yp11(:,j)/(RHPz11(i)-RHPp11(j));
        end
    end
    %     Qp11=(Yp'*Yp).*(1./(diag(RHPp11')*ones(np11)+ones(np11)*diag(RHPp11)));
    %     Qz11=(Yz'*Yz).*(1./(diag(RHPz11)*ones(nz11)+ones(nz11)*diag(RHPz11')));
    %     Qzp11=(Yz'*Yp).*(1./(diag(RHPz11)*ones(nz11,np11)-ones(nz11,np11)*diag(RHPp11)));
    Msmin11=sqrt(1+norm(sqrtm(inv(Qz11))*Qzp11*sqrtm(inv(Qp11)))^2);
end
%% *** Bounds on KS ****
%To find the bound for KS the mirror image of the antistable part of the
%G`s are used. If a zero is very close to be negativ, this zero will be
%implemented in the antistable part.

%Split G into stable and anti stable part, by using a matlab
%function.
[GSu1, GNSu1]=gnsi(G1);
[GSu2, GNSu2]=gnsi(G2);
[GSu3, GNSu3]=gnsi(G3);
[GSu4, GNSu4]=gnsi(G4);
[GSu5, GNSu5]=gnsi(G5);
[GSu6, GNSu6]=gnsi(G6);
[GSu7, GNSu7]=gnsi(G7);
[GSu8, GNSu8]=gnsi(G8);
[GSu9, GNSu9]=gnsi(G9);
[GSu10, GNSu10]=gnsi(G10);
[GSu11, GNSu11]=gnsi(G11);


%Lower bound on KS for the output P_bh
tf1=(GNSu1)';
h1=hsvd(tf1);
KSmin1=1/min(h1);

%Lower bound on KS for the output P_wh
tf2=(GNSu2)';
h2=hsvd(tf2);
KSmin2=1/min(h2);

%Lower bound on KS for the output W_in
tf3=(GNSu3)';
h3=hsvd(tf3);
KSmin3=1/min(h3);

%Lower bound on KS for the output P_in
tf4=(GNSu4)';
h4=hsvd(tf4);
KSmin4=1/min(h4);

%Lower bound on KS for the output P_rb
tf5=(GNSu5)';
h5=hsvd(tf5);
KSmin5=1/min(h5);

%Lower bound on KS for the output DPr
tf6=(GNSu6)';
h6=hsvd(tf6);
KSmin6=1/min(h6);

%Lower bound on KS for the output P_t
tf7=(GNSu7)';
h7=hsvd(tf7);
KSmin7=1/min(h7);

%Lower bound on KS for the output Q
tf8=(GNSu8)';
h8=hsvd(tf8);
KSmin8=1/min(h8);

%Lower bound on KS for the output W
tf9=(GNSu9)';
h9=hsvd(tf9);
KSmin9=1/min(h9);

%Lower bound on KS for the output Rho_t
tf10=(GNSu10)';
h10=hsvd(tf10);
KSmin10=1/min(h10);

%Lower bound on KS for the output Alpha_L
tf11=(GNSu11)';
h11=hsvd(tf11);
KSmin11=1/min(h11);
%% *** bounds on KSGd ****
s=tf('s') ;
%There are two disturbances for each output, and they will be considered
%separately. d1=d1 and d2=d2

%Lower bound on KSGd1 for the output P_bh
solve1=((inv(Gdms11))*G1);
[t11s t11ns]=gnsi(solve1);
tfKSd11=(t11ns)';
hd11=hsvd(tfKSd11);
KSGd11=1/min(hd11);

%Lower bound on KSGd2 for the output P_bh
solve2=((inv(Gdms12))*G1);
[t12s t12ns]=gnsi(solve2);
tfKSd12=(t12ns)';
hd12=hsvd(tfKSd12);
KSGd12=1/min(hd12);

%Lower bound on KSGd1 for the output P_wh
solve3=((inv(Gdms21))*G2);
[t21s t21ns]=gnsi(solve3);
tfKSd21=(t21ns)';
hd21=hsvd(tfKSd21);
KSGd21=1/min(hd21);

%Lower bound on KSGd2 for the output P_wh
solve4=((inv(Gdms22))*G2);
[t22s t22ns]=gnsi(solve4);
tfKSd22=(t22ns)';
hd22=hsvd(tfKSd22);
KSGd22=1/min(hd22);

%Lower bound on KSGd1 for the output W_in
solve5=((inv(Gdms31))*G3);
[t31s t31ns]=gnsi(solve5);
tfKSd31=(t31ns)';
hd31=hsvd(tfKSd31);
KSGd31=1/min(hd31);

%Lower bound on KSGd2 for the output W_in
solve6=((inv(Gdms32))*G3);
[t32s t32ns]=gnsi(solve6);
tfKSd32=(t32ns)';
hd32=hsvd(tfKSd32);
KSGd32=1/min(hd32);

%Lower bound on KSGd1 for the output P_in
solve7=((inv(Gdms41))*G4);
[t41s t41ns]=gnsi(solve7);
tfKSd41=(t41ns)';
hd41=hsvd(tfKSd41);
KSGd41=1/min(hd41);

%Lower bound on KSGd1 for the output P_in
solve8=((inv(Gdms42))*G4);
[t42s t42ns]=gnsi(solve8);
tfKSd42=(t42ns)';
hd42=hsvd(tfKSd42);
KSGd42=1/min(hd42);

%Lower bound on KSGd1 for the output P_rb
solve9=((inv(Gdms51))*G5);
[t51s t51ns]=gnsi(solve9);
tfKSd51=(t51ns)';
hd51=hsvd(tfKSd51);
KSGd51=1/min(hd51);

%Lower bound on KSGd1 for the output P_rb
solve10=((inv(Gdms52))*G5);
[t52s t52ns]=gnsi(solve10);
tfKSd52=(t52ns)';
hd52=hsvd(tfKSd52);
KSGd52=1/min(hd52);

%Lower bound on KSGd1 for the output DP_r
solve11=((inv(Gdms61))*G6);
[t61s t61ns]=gnsi(solve11);
tfKSd61=(t61ns)';
hd61=hsvd(tfKSd61);
KSGd61=1/min(hd61);

%Lower bound on KSGd2 for the output DP_r
solve12=((inv(Gdms62))*G6);
[t62s t62ns]=gnsi(solve12);
tfKSd62=(t62ns)';
hd62=hsvd(tfKSd62);
KSGd62=1/min(hd62);

%Lower bound on KSGd1 for the output P_t
solve13=((inv(Gdms71))*G7);
[t71s t71ns]=gnsi(solve13);
tfKSd71=(t71ns)';
hd71=hsvd(tfKSd71);
KSGd71=1/min(hd71);

%Lower bound on KSGd2 for the output P_t
solve14=((inv(Gdms72))*G7);
[t72s t72ns]=gnsi(solve14);
tfKSd72=(t72ns)';
hd72=hsvd(tfKSd72);
KSGd72=1/min(hd72);

%Lower bound on KSGd1 for the output Q_out
%Put in one stable zero, so gnsi works
solve15=((inv(Gdms81))*G8);
[t81s t81ns]=gnsi(solve5);
tfKSd81=(t81ns)';
hd81=hsvd(tfKSd81);
KSGd81=1/min(hd81);

%Lower bound on KSGd2 for the output Q_out
solve16=((inv(Gdms82))*G8);
[t82s t82ns]=gnsi(solve16);
tfKSd82=(t82ns)';
hd82=hsvd(tfKSd82);
KSGd82=1/min(hd82);

%Lower bound on KSGd1 for the output W_out
solve17=(inv(Gdms91))*G9;
[t91s t91ns]=gnsi(solve17);
tfKSd91=(t91ns)';
hd91=hsvd(tfKSd91);
KSGd91=1/min(hd91);

%Lower bound on KSGd2 for the output W_out
solve18=(inv(Gdms92))*G9;
[t92s t92ns]=gnsi(solve18);
tfKSd92=(t92ns)';
hd92=hsvd(tfKSd92);
KSGd92=1/min(hd92);

%Lower bound on KSGd1 for the output Rho_t
solve19=((inv(Gdms101))*G10);
[t101s t101ns]=gnsi(solve19);
tfKSd101=(t101ns)';
hd101=hsvd(tfKSd101);
KSGd101=1/min(hd101);

%Lower bound on KSGd2 for the output Rho_t
solve20=((inv(Gdms102))*G10);
[t102s t102ns]=gnsi(solve20);
tfKSd102=(t102ns)';
hd102=hsvd(tfKSd102);
KSGd102=1/min(hd102);

%Lower bound on KSGd1 for the output Alpa_L
solve21=((inv(Gdms111))*G11);
[t111s t111ns]=gnsi(solve21);
tfKSd111=(t111ns)';
hd111=hsvd(tfKSd111);
KSGd111=1/min(hd111);

%Lower bound on KSGd2 for the output Alpa_L
solve22=((inv(Gdms112))*G11);
[t112s t112ns]=gnsi(solve22);
tfKSd112=(t112ns)';
hd112=hsvd(tfKSd112);
KSGd112=1/min(hd112);

%% SG *****
%Used the poles and zeros that made the highest bound, since this bound is
%only tight for one RHP-pole and one RHP-zero.
%% SG for P_bh
[p, z1]=pzmap(sys1);
RHPp1=p(find(p>0));
RHPz1=z1(find(z1>0));
l = 1;
if (length(RHPp1)>0 && length(RHPz1)>0)
    np1=length(RHPp1); nz1=length(RHPz1);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp1 =  zeros(l,np1);
    Yz1 = zeros(l,nz1);
    for i=1:np1
        Yp1(:,i) = Cs(1,:)*V(:,i)/norm(Cs(1,:)*V(:,i)); %Pole direction
    end
    for j=1:nz1
        [U,S,V]=svd(evalfr(sys1,RHPz1(j))); Yz1(:,j)=U(:,end); %zero direction
    end
    
    SGmin = zeros(np1,nz1);
    for i=1:np1
        for j=1:nz1
            num1=abs((ctranspose(Yz1(:,j))*(evalfr(Gms1,RHPz1(j))))*((evalfr(Gms1,RHPp1(i)))^-1)*Yp1(:,i));
            den1=norm(ctranspose(Yz1(:,j))*(evalfr(Gms1,RHPz1(j))),2)*norm(((evalfr(Gms1,RHPp1(i)))^-1)*Yp1(:,i),2);
            cos2=(num1/den1)^2;
            sin2=1-cos2;
            
            term1 = (abs(RHPz1(j)+conj(RHPp1(i))))^2/(abs(RHPz1(j)-RHPp1(i)))^2;
            SGmin(i,j)=norm((ctranspose(Yz1(:,j))*(evalfr(Gms1,RHPz1(j)))),2)*sqrt(sin2+(term1*cos2));
        end
    end
    SGmin1 = max(max((SGmin)));
else
    SGmin1=0;
end

%% SG for P_wh
[p, z2]=pzmap(sys2);
RHPp2=p(find(p>0));
RHPz2=z2(find(z2>0));

if length(RHPp2)>0 && length(RHPz2)>0
    np2=length(RHPp2); nz2=length(RHPz2);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp2 =  zeros(l,np2);
    Yz2 = zeros(l,nz2);
    for i=1:np2
        Yp2(:,i) = Cs(2,:)*V(:,i)/norm(Cs(2,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz2
        [U,S,V]=svd(evalfr(sys2,RHPz2(j))); Yz2(:,j)=U(:,end); %zero direction
    end
    SGmin = zeros(np2,nz2);
    for i=1:np2
        for j=1:nz2
            num2=abs((ctranspose(Yz2(:,j))*(evalfr(Gms2,RHPz2(j))))*((evalfr(Gms2,RHPp2(i)))^-1)*Yp2(:,i));
            den2=norm(ctranspose(Yz2(:,j))*(evalfr(Gms2,RHPz2(j))),2)*norm(((evalfr(Gms2,RHPp2(i)))^-1)*Yp2(:,i),2);
            cos2=(num2/den2)^2;
            sin2=1-cos2;
            
            term2 = (abs(RHPz2(j)+conj(RHPp2(i))))^2/(abs(RHPz2(j)-RHPp2(i)))^2;
            SGmin(i,j)=norm((ctranspose(Yz2(:,j))*(evalfr(Gms2,RHPz2(j)))),2)*sqrt(sin2+(term2*cos2));
        end
    end
    SGmin2 = max(max((SGmin)));
else
    SGmin2=0;
end

%% SG for W_in
[p, z3]=pzmap(sys3);
RHPp3=p(find(p>0));
RHPz3=z3(find(z3>0));

if length(RHPp3)>0 && length(RHPz3)>0
    np3=length(RHPp3); nz3=length(RHPz3);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp3 =  zeros(l,np3);
    Yz3 = zeros(l,nz3);
    for i=1:np3
        Yp3(:,i) = Cs(3,:)*V(:,i)/norm(Cs(3,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz3
        [U,S,V]=svd(evalfr(sys3,RHPz3(j))); Yz3(:,j)=U(:,end); %zero direction
    end
    SGmin = zeros(np3,nz3);
    for i=1:np3
        for j=1:nz3
            num3=abs((ctranspose(Yz3(:,j))*(evalfr(Gms3,RHPz3(j))))*((evalfr(Gms3,RHPp3(i)))^-1)*Yp3(:,i));
            den3=norm(ctranspose(Yz3(:,j))*(evalfr(Gms3,RHPz3(j))),2)*norm(((evalfr(Gms3,RHPp3(i)))^-1)*Yp3(:,i),2);
            cos2=(num3/den3)^2;
            sin2=1-cos2;
            
            term3 = (abs(RHPz3(j)+conj(RHPp3(i))))^2/(abs(RHPz3(j)-RHPp3(i)))^2;
            SGmin(i,j)=norm((ctranspose(Yz3(:,j))*(evalfr(Gms3,RHPz3(j)))),2)*sqrt(sin2+(term3*cos2));
        end
    end
    SGmin3 = max(max((SGmin)));
else
    SGmin3=0;
end

%% SG for P_in
[p, z4]=pzmap(sys4);
RHPp4=p(find(p>0));
RHPz4=z4(find(z4>0));

if length(RHPp4)>0 && length(RHPz4)>0
    np4=length(RHPp4); nz4=length(RHPz4);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp4 =  zeros(l,np4);
    Yz4 = zeros(l,nz4);
    for i=1:np4
        Yp4(:,i) = Cs(4,:)*V(:,i)/norm(Cs(4,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz4
        [U,S,V]=svd(evalfr(sys4,RHPz4(j))); Yz4(:,j)=U(:,end); %zero direction
    end
    SGmin = zeros(np4,nz4);
    for i=1:np4
        for j=1:nz4
            num4=abs((ctranspose(Yz4(:,j))*(evalfr(Gms4,RHPz4(j))))*((evalfr(Gms4,RHPp4(i)))^-1)*Yp4(:,i));
            den4=norm(ctranspose(Yz4(:,j))*(evalfr(Gms4,RHPz4(j))),2)*norm(((evalfr(Gms4,RHPp4(i)))^-1)*Yp4(:,i),2);
            cos2=(num4/den4)^2;
            sin2=1-cos2;
            
            term4 = (abs(RHPz4(j)+conj(RHPp4(i))))^2/(abs(RHPz4(j)-RHPp4(i)))^2;
            SGmin(i,j)=norm((ctranspose(Yz4(:,j))*(evalfr(Gms4,RHPz4(j)))),2)*sqrt(sin2+(term4*cos2));
        end
    end
    SGmin4 = max(max((SGmin)));
else
    SGmin4=0;
end

%% SG for P_rb
[p, z5]=pzmap(sys5);
RHPp5=p(find(p>0));
RHPz5=z5(find(z5>0));

if length(RHPp5)>0 && length(RHPz5)>0
    np5=length(RHPp5); nz5=length(RHPz5);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp5 =  zeros(l,np5);
    Yz5 = zeros(l,nz5);
    for i=1:np5
        Yp5(:,i) = Cs(5,:)*V(:,i)/norm(Cs(5,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz5
        [U,S,V]=svd(evalfr(sys5,RHPz5(j))); Yz5(:,j)=U(:,end); %zero direction
    end
    SGmin = zeros(np5,nz5);
    for i=1:np5
        for j=1:nz5
            num5=abs((ctranspose(Yz5(:,j))*(evalfr(Gms5,RHPz5(j))))*((evalfr(Gms5,RHPp5(i)))^-1)*Yp5(:,i));
            den5=norm(ctranspose(Yz5(:,j))*(evalfr(Gms5,RHPz5(j))),2)*norm(((evalfr(Gms5,RHPp5(i)))^-1)*Yp5(:,i),2);
            cos2=(num5/den5)^2;
            sin2=1-cos2;
            
            term5 = (abs(RHPz5(j)+conj(RHPp5(i))))^2/(abs(RHPz5(j)-RHPp5(i)))^2;
            SGmin(i,j)=norm((ctranspose(Yz5(:,j))*(evalfr(Gms5,RHPz5(j)))),2)*sqrt(sin2+(term5*cos2));
        end
    end
    SGmin5 = max(max((SGmin)));
else
    SGmin5=0;
end

%% SG for DP_r
[p, z6]=pzmap(sys6);
RHPp6=p(find(p>0));
RHPz6=z6(find(z6>0));

if length(RHPp6)>0 && length(RHPz6)>0
    np6=length(RHPp6); nz6=length(RHPz6);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp6 =  zeros(l,np6);
    Yz6 = zeros(l,nz6);
    for i=1:np6
        Yp6(:,i) = Cs(6,:)*V(:,i)/norm(Cs(6,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz6
        [U,S,V]=svd(evalfr(sys6,RHPz6(j))); Yz6(:,j)=U(:,end); %zero direction
    end
    SGmin = zeros(np6,nz6);
    for i=1:np6
        for j=1:nz6
            num6=abs((ctranspose(Yz6(:,j))*(evalfr(Gms6,RHPz6(j))))*((evalfr(Gms6,RHPp6(i)))^-1)*Yp6(:,i));
            den6=norm(ctranspose(Yz6(:,j))*(evalfr(Gms6,RHPz6(j))),2)*norm(((evalfr(Gms6,RHPp6(i)))^-1)*Yp6(:,i),2);
            cos2=(num6/den6)^2;
            sin2=1-cos2;
            
            term6 = (abs(RHPz6(j)+conj(RHPp6(i))))^2/(abs(RHPz6(j)-RHPp6(i)))^2;
            SGmin(i,j)=norm((ctranspose(Yz6(:,j))*(evalfr(Gms6,RHPz6(j)))),2)*sqrt(sin2+(term6*cos2));
        end
    end
    SGmin6 = max(max((SGmin)));
else
    SGmin6=0;
end

%% SG for P_t
[p, z7]=pzmap(sys7);
RHPp7=p(find(p>0));
RHPz7=z7(find(z7>0));

if length(RHPp7)>0 && length(RHPz7)>0
    np7=length(RHPp7); nz7=length(RHPz7);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp7 =  zeros(l,np7);
    Yz7 = zeros(l,nz7);
    for i=1:np7
        Yp7(:,i) = Cs(7,:)*V(:,i)/norm(Cs(7,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz7
        [U,S,V]=svd(evalfr(sys7,RHPz7(j))); Yz7(:,j)=U(:,end); %zero direction
    end
    SGmin = zeros(np7,nz7);
    for i=1:np7
        for j=1:nz7
            num7=abs((ctranspose(Yz7(:,j))*(evalfr(Gms7,RHPz7(j))))*((evalfr(Gms7,RHPp7(i)))^-1)*Yp7(:,i));
            den7=norm(ctranspose(Yz7(:,j))*(evalfr(Gms7,RHPz7(j))),2)*norm(((evalfr(Gms7,RHPp7(i)))^-1)*Yp7(:,i),2);
            cos2=(num7/den7)^2;
            sin2=1-cos2;
            
            term7 = (abs(RHPz7(j)+conj(RHPp7(i))))^2/(abs(RHPz7(j)-RHPp7(i)))^2;
            SGmin(i,j)=norm((ctranspose(Yz7(:,j))*(evalfr(Gms7,RHPz7(j)))),2)*sqrt(sin2+(term7*cos2));
        end
    end
    SGmin7 = max(max((SGmin)));
else
    SGmin7=0;
end

%% SG for Q_out
[p, z8]=pzmap(sys8);
RHPp8=p(find(p>0));
RHPz8=z8(find(z8>0));

if length(RHPp8)>0 && length(RHPz8)>0
    np8=length(RHPp8); nz8=length(RHPz8);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp8 =  zeros(l,np8);
    Yz8 = zeros(l,nz8);
    for i=1:np8
        Yp8(:,i) = Cs(8,:)*V(:,i)/norm(Cs(8,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz8
        [U,S,V]=svd(evalfr(sys8,RHPz8(j))); Yz8(:,j)=U(:,end); %zero direction
    end
    SGmin = zeros(np8,nz8);
    for i=1:np8
        for j=1:nz8
            num8=abs((ctranspose(Yz8(:,j))*(evalfr(Gms8,RHPz8(j))))*((evalfr(Gms8,RHPp8(i)))^-1)*Yp8(:,i));
            den8=norm(ctranspose(Yz8(:,j))*(evalfr(Gms8,RHPz8(j))),2)*norm(((evalfr(Gms8,RHPp8(i)))^-1)*Yp8(:,i),2);
            cos2=(num8/den8)^2;
            sin2=1-cos2;
            
            term8 = (abs(RHPz8(j)+conj(RHPp8(i))))^2/(abs(RHPz8(j)-RHPp8(i)))^2;
            SGmin(i,j)=norm((ctranspose(Yz8(:,j))*(evalfr(Gms8,RHPz8(j)))),2)*sqrt(sin2+(term8*cos2));
        end
    end
    SGmin8 = max(max((SGmin)));
else
    SGmin8=0;
end

%% SG for W_out
[p, z9]=pzmap(sys9);
RHPp9=p(find(p>0));
RHPz9=z9(find(z9>0));

if length(RHPp9)>0 && length(RHPz9)>0
    np9=length(RHPp9); nz9=length(RHPz9);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp9 =  zeros(l,np9);
    Yz9 = zeros(l,nz9);
    for i=1:np9
        Yp9(:,i) = Cs(9,:)*V(:,i)/norm(Cs(9,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz9
        [U,S,V]=svd(evalfr(sys9,RHPz9(j))); Yz9(:,j)=U(:,end); %zero direction
    end
    SGmin = zeros(np9,nz9);
    for i=1:np9
        for j=1:nz9
            num9=abs((ctranspose(Yz9(:,j))*(evalfr(Gms9,RHPz9(j))))*((evalfr(Gms9,RHPp9(i)))^-1)*Yp9(:,i));
            den9=norm(ctranspose(Yz9(:,j))*(evalfr(Gms9,RHPz9(j))),2)*norm(((evalfr(Gms9,RHPp9(i)))^-1)*Yp9(:,i),2);
            cos2=(num9/den9)^2;
            sin2=1-cos2;
            
            term9 = (abs(RHPz9(j)+conj(RHPp9(i))))^2/(abs(RHPz9(j)-RHPp9(i)))^2;
            SGmin(i,j)=norm((ctranspose(Yz9(:,j))*(evalfr(Gms9,RHPz9(j)))),2)*sqrt(sin2+(term9*cos2));
        end
    end
    SGmin9 = max(max((SGmin)));
else
    SGmin9=0;
end
%% SG for Rho_t
[p, z10]=pzmap(sys10);
RHPp10=p(find(p>0));
RHPz10=z10(find(z10>0));

if length(RHPp10)>0 && length(RHPz10)>0
    np10=length(RHPp10); nz10=length(RHPz10);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp10 =  zeros(l,np10);
    Yz10 = zeros(l,nz10);
    for i=1:np10
        Yp10(:,i) = Cs(10,:)*V(:,i)/norm(Cs(10,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz10
        [U,S,V]=svd(evalfr(sys10,RHPz10(j))); Yz10(:,j)=U(:,end); %zero direction
    end
    SGmin = zeros(np10,nz10);
    for i=1:np10
        for j=1:nz10
            num10=abs((ctranspose(Yz10(:,j))*(evalfr(Gms10,RHPz10(j))))*((evalfr(Gms10,RHPp10(i)))^-1)*Yp10(:,i));
            den10=norm(ctranspose(Yz10(:,j))*(evalfr(Gms10,RHPz10(j))),2)*norm(((evalfr(Gms10,RHPp10(i)))^-1)*Yp10(:,i),2);
            cos2=(num10/den10)^2;
            sin2=1-cos2;
            
            term10 = (abs(RHPz10(j)+conj(RHPp10(i))))^2/(abs(RHPz10(j)-RHPp10(i)))^2;
            SGmin(i,j)=norm((ctranspose(Yz10(:,j))*(evalfr(Gms10,RHPz10(j)))),2)*sqrt(sin2+(term10*cos2));
        end
    end
    SGmin10 = max(max((SGmin)));
else
    SGmin10=0;
end
%% SG for Alpha_L
[p, z11]=pzmap(sys11);
RHPp11=p(find(p>0));
RHPz11=z11(find(z11>0));

if length(RHPp11)>0 && length(RHPz11)>0
    np11=length(RHPp11); nz11=length(RHPz11);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp11 =  zeros(l,np11);
    Yz11 = zeros(l,nz11);
    for i=1:np11
        Yp11(:,i) = Cs(11,:)*V(:,i)/norm(Cs(11,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz11
        [U,S,V]=svd(evalfr(sys11,RHPz11(j))); Yz11(:,j)=U(:,end); %zero direction
    end
    SGmin = zeros(np11,nz11);
    for i=1:np11
        for j=1:nz11
            num11=abs((ctranspose(Yz11(:,j))*(evalfr(Gms11,RHPz11(j))))*((evalfr(Gms11,RHPp11(i)))^-1)*Yp11(:,i));
            den11=norm(ctranspose(Yz11(:,j))*(evalfr(Gms11,RHPz11(j))),2)*norm(((evalfr(Gms11,RHPp11(i)))^-1)*Yp11(:,i),2);
            cos2=(num11/den11)^2;
            sin2=1-cos2;
            
            term11 = (abs(RHPz11(j)+conj(RHPp11(i))))^2/(abs(RHPz11(j)-RHPp11(i)))^2;
            SGmin(i,j)=norm((ctranspose(Yz11(:,j))*(evalfr(Gms11,RHPz11(j)))),2)*sqrt(sin2+(term11*cos2));
        end
    end
    SGmin11 = max(max((SGmin)));
else
    SGmin11=0;
end
%% ***** SGd ********
%% SGd1 for P_bh
[p, z1]=pzmap(sys1);
RHPp1=p(find(p>0));
RHPz1=z1(find(z1>0));
l = 1;
if (length(RHPp1)>0 && length(RHPz1)>0)
    np1=length(RHPp1); nz1=length(RHPz1);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp1 =  zeros(l,np1);
    Yz1 = zeros(l,nz1);
    for i=1:np1
        Yp1(:,i) = Cs(1,:)*V(:,i)/norm(Cs(1,:)*V(:,i)); %Pole direction
    end
    for j=1:nz1
        [U,S,V]=svd(evalfr(sys1,RHPz1(j))); Yz1(:,j)=U(:,end); %zero direction
    end
    
    SGdmin = zeros(np1,nz1);
    for i=1:np1
        for j=1:nz1
            num1=abs((ctranspose(Yz1(:,j))*(evalfr(Gdms11,RHPz1(j))))*((evalfr(Gdms11,RHPp1(i)))^-1)*Yp1(:,i));
            den1=norm(ctranspose(Yz1(:,j))*(evalfr(Gdms11,RHPz1(j))),2)*norm(((evalfr(Gdms11,RHPp1(i)))^-1)*Yp1(:,i),2);
            cos2=(num1/den1)^2;
            sin2=1-cos2;
            
            term1 = (abs(RHPz1(j)+conj(RHPp1(i))))^2/(abs(RHPz1(j)-RHPp1(i)))^2;
            SGdmin(i,j)=norm((ctranspose(Yz1(:,j))*(evalfr(Gdms11,RHPz1(j)))),2)*sqrt(sin2+(term1*cos2));
        end
    end
    SGdmin11 = max(max((SGdmin)));
else
    SGdmin11=0;
end

%% SGd2 for P_bh
[p, z1]=pzmap(sys1);
RHPp1=p(find(p>0));
RHPz1=z1(find(z1>0));
l = 1;
if (length(RHPp1)>0 && length(RHPz1)>0)
    np1=length(RHPp1); nz1=length(RHPz1);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp1 =  zeros(l,np1);
    Yz1 = zeros(l,nz1);
    for i=1:np1
        Yp1(:,i) = Cs(1,:)*V(:,i)/norm(Cs(1,:)*V(:,i)); %Pole direction
    end
    for j=1:nz1
        [U,S,V]=svd(evalfr(sys1,RHPz1(j))); Yz1(:,j)=U(:,end); %zero direction
    end
    
    SGdmin = zeros(np1,nz1);
    for i=1:np1
        for j=1:nz1
            num1=abs((ctranspose(Yz1(:,j))*(evalfr(Gdms12,RHPz1(j))))*((evalfr(Gdms12,RHPp1(i)))^-1)*Yp1(:,i));
            den1=norm(ctranspose(Yz1(:,j))*(evalfr(Gdms12,RHPz1(j))),2)*norm(((evalfr(Gdms12,RHPp1(i)))^-1)*Yp1(:,i),2);
            cos2=(num1/den1)^2;
            sin2=1-cos2;
            
            term1 = (abs(RHPz1(j)+conj(RHPp1(i))))^2/(abs(RHPz1(j)-RHPp1(i)))^2;
            SGdmin(i,j)=norm((ctranspose(Yz1(:,j))*(evalfr(Gdms12,RHPz1(j)))),2)*sqrt(sin2+(term1*cos2));
        end
    end
    SGdmin12 = max(max((SGdmin)));
else
    SGdmin12=0;
end
%% SGd1 for P_wh
[p, z2]=pzmap(sys2);
RHPp2=p(find(p>0));
RHPz2=z2(find(z2>0));

if length(RHPp2)>0 && length(RHPz2)>0
    np2=length(RHPp2); nz2=length(RHPz2);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp2 =  zeros(l,np2);
    Yz2 = zeros(l,nz2);
    for i=1:np2
        Yp2(:,i) = Cs(2,:)*V(:,i)/norm(Cs(2,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz2
        [U,S,V]=svd(evalfr(sys2,RHPz2(j))); Yz2(:,j)=U(:,end); %zero direction
    end
    SGdmin = zeros(np2,nz2);
    for i=1:np2
        for j=1:nz2
            num2=abs((ctranspose(Yz2(:,j))*(evalfr(Gdms21,RHPz2(j))))*((evalfr(Gdms21,RHPp2(i)))^-1)*Yp2(:,i));
            den2=norm(ctranspose(Yz2(:,j))*(evalfr(Gdms21,RHPz2(j))),2)*norm(((evalfr(Gdms21,RHPp2(i)))^-1)*Yp2(:,i),2);
            cos2=(num2/den2)^2;
            sin2=1-cos2;
            
            term2 = (abs(RHPz2(j)+conj(RHPp2(i))))^2/(abs(RHPz2(j)-RHPp2(i)))^2;
            SGdmin(i,j)=norm((ctranspose(Yz2(:,j))*(evalfr(Gdms21,RHPz2(j)))),2)*sqrt(sin2+(term2*cos2));
        end
    end
    SGdmin21 = max(max((SGdmin)));
else
    SGdmin21=0;
end
%% SGd2 for P_wh
[p, z2]=pzmap(sys2);
RHPp2=p(find(p>0));
RHPz2=z2(find(z2>0));

if length(RHPp2)>0 && length(RHPz2)>0
    np2=length(RHPp2); nz2=length(RHPz2);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp2 =  zeros(l,np2);
    Yz2 = zeros(l,nz2);
    for i=1:np2
        Yp2(:,i) = Cs(2,:)*V(:,i)/norm(Cs(2,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz2
        [U,S,V]=svd(evalfr(sys2,RHPz2(j))); Yz2(:,j)=U(:,end); %zero direction
    end
    SGdmin = zeros(np2,nz2);
    for i=1:np2
        for j=1:nz2
            num2=abs((ctranspose(Yz2(:,j))*(evalfr(Gdms22,RHPz2(j))))*((evalfr(Gdms22,RHPp2(i)))^-1)*Yp2(:,i));
            den2=norm(ctranspose(Yz2(:,j))*(evalfr(Gdms22,RHPz2(j))),2)*norm(((evalfr(Gdms22,RHPp2(i)))^-1)*Yp2(:,i),2);
            cos2=(num2/den2)^2;
            sin2=1-cos2;
            
            term2 = (abs(RHPz2(j)+conj(RHPp2(i))))^2/(abs(RHPz2(j)-RHPp2(i)))^2;
            SGdmin(i,j)=norm((ctranspose(Yz2(:,j))*(evalfr(Gdms22,RHPz2(j)))),2)*sqrt(sin2+(term2*cos2));
        end
    end
    SGdmin22 = max(max((SGdmin)));
else
    SGdmin22=0;
end
%% SGd1 for W_in
[p, z3]=pzmap(sys3);
RHPp3=p(find(p>0));
RHPz3=z3(find(z3>0));

if length(RHPp3)>0 && length(RHPz3)>0
    np3=length(RHPp3); nz3=length(RHPz3);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp3 =  zeros(l,np3);
    Yz3 = zeros(l,nz3);
    for i=1:np3
        Yp3(:,i) = Cs(3,:)*V(:,i)/norm(Cs(3,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz3
        [U,S,V]=svd(evalfr(sys3,RHPz3(j))); Yz3(:,j)=U(:,end); %zero direction
    end
    SGdmin = zeros(np3,nz3);
    for i=1:np3
        for j=1:nz3
            num3=abs((ctranspose(Yz3(:,j))*(evalfr(Gdms31,RHPz3(j))))*((evalfr(Gdms31,RHPp3(i)))^-1)*Yp3(:,i));
            den3=norm(ctranspose(Yz3(:,j))*(evalfr(Gdms31,RHPz3(j))),2)*norm(((evalfr(Gdms31,RHPp3(i)))^-1)*Yp3(:,i),2);
            cos2=(num3/den3)^2;
            sin2=1-cos2;
            
            term3 = (abs(RHPz3(j)+conj(RHPp3(i))))^2/(abs(RHPz3(j)-RHPp3(i)))^2;
            SGdmin(i,j)=norm((ctranspose(Yz3(:,j))*(evalfr(Gdms31,RHPz3(j)))),2)*sqrt(sin2+(term3*cos2));
        end
    end
    SGdmin31 = max(max((SGdmin)));
else
    SGdmin31=0;
end
%% SGd2 for W_in
[p, z3]=pzmap(sys3);
RHPp3=p(find(p>0));
RHPz3=z3(find(z3>0));

if length(RHPp3)>0 && length(RHPz3)>0
    np3=length(RHPp3); nz3=length(RHPz3);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp3 =  zeros(l,np3);
    Yz3 = zeros(l,nz3);
    for i=1:np3
        Yp3(:,i) = Cs(3,:)*V(:,i)/norm(Cs(3,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz3
        [U,S,V]=svd(evalfr(sys3,RHPz3(j))); Yz3(:,j)=U(:,end); %zero direction
    end
    SGdmin = zeros(np3,nz3);
    for i=1:np3
        for j=1:nz3
            num3=abs((ctranspose(Yz3(:,j))*(evalfr(Gdms32,RHPz3(j))))*((evalfr(Gdms32,RHPp3(i)))^-1)*Yp3(:,i));
            den3=norm(ctranspose(Yz3(:,j))*(evalfr(Gdms32,RHPz3(j))),2)*norm(((evalfr(Gdms32,RHPp3(i)))^-1)*Yp3(:,i),2);
            cos2=(num3/den3)^2;
            sin2=1-cos2;
            
            term3 = (abs(RHPz3(j)+conj(RHPp3(i))))^2/(abs(RHPz3(j)-RHPp3(i)))^2;
            SGdmin(i,j)=norm((ctranspose(Yz3(:,j))*(evalfr(Gdms32,RHPz3(j)))),2)*sqrt(sin2+(term3*cos2));
        end
    end
    SGdmin32 = max(max((SGdmin)));
else
    SGdmin32=0;
end
%% SGd1 for P_in
[p, z4]=pzmap(sys4);
RHPp4=p(find(p>0));
RHPz4=z4(find(z4>0));

if length(RHPp4)>0 && length(RHPz4)>0
    np4=length(RHPp4); nz4=length(RHPz4);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp4 =  zeros(l,np4);
    Yz4 = zeros(l,nz4);
    for i=1:np4
        Yp4(:,i) = Cs(4,:)*V(:,i)/norm(Cs(4,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz4
        [U,S,V]=svd(evalfr(sys4,RHPz4(j))); Yz4(:,j)=U(:,end); %zero direction
    end
    SGdmin = zeros(np4,nz4);
    for i=1:np4
        for j=1:nz4
            num4=abs((ctranspose(Yz4(:,j))*(evalfr(Gdms41,RHPz4(j))))*((evalfr(Gdms41,RHPp4(i)))^-1)*Yp4(:,i));
            den4=norm(ctranspose(Yz4(:,j))*(evalfr(Gdms41,RHPz4(j))),2)*norm(((evalfr(Gdms41,RHPp4(i)))^-1)*Yp4(:,i),2);
            cos2=(num4/den4)^2;
            sin2=1-cos2;
            
            term4 = (abs(RHPz4(j)+conj(RHPp4(i))))^2/(abs(RHPz4(j)-RHPp4(i)))^2;
            SGdmin(i,j)=norm((ctranspose(Yz4(:,j))*(evalfr(Gdms41,RHPz4(j)))),2)*sqrt(sin2+(term4*cos2));
        end
    end
    SGdmin41 = max(max((SGdmin)));
else
    SGdmin41 = 0;
end
%% SGd2 for P_in
[p, z4]=pzmap(sys4);
RHPp4=p(find(p>0));
RHPz4=z4(find(z4>0));

if length(RHPp4)>0 && length(RHPz4)>0
    np4=length(RHPp4); nz4=length(RHPz4);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp4 =  zeros(l,np4);
    Yz4 = zeros(l,nz4);
    for i=1:np4
        Yp4(:,i) = Cs(4,:)*V(:,i)/norm(Cs(4,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz4
        [U,S,V]=svd(evalfr(sys4,RHPz4(j))); Yz4(:,j)=U(:,end); %zero direction
    end
    SGdmin = zeros(np4,nz4);
    for i=1:np4
        for j=1:nz4
            num4=abs((ctranspose(Yz4(:,j))*(evalfr(Gdms42,RHPz4(j))))*((evalfr(Gdms42,RHPp4(i)))^-1)*Yp4(:,i));
            den4=norm(ctranspose(Yz4(:,j))*(evalfr(Gdms42,RHPz4(j))),2)*norm(((evalfr(Gdms42,RHPp4(i)))^-1)*Yp4(:,i),2);
            cos2=(num4/den4)^2;
            sin2=1-cos2;
            
            term4 = (abs(RHPz4(j)+conj(RHPp4(i))))^2/(abs(RHPz4(j)-RHPp4(i)))^2;
            SGdmin(i,j)=norm((ctranspose(Yz4(:,j))*(evalfr(Gdms42,RHPz4(j)))),2)*sqrt(sin2+(term4*cos2));
        end
    end
    SGdmin42 = max(max((SGdmin)));
else
    SGdmin42 = 0;
end
%% SGd1 for P_rb
[p, z5]=pzmap(sys5);
RHPp5=p(find(p>0));
RHPz5=z5(find(z5>0));

if length(RHPp5)>0 && length(RHPz5)>0
    np5=length(RHPp5); nz5=length(RHPz5);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp5 =  zeros(l,np5);
    Yz5 = zeros(l,nz5);
    for i=1:np5
        Yp5(:,i) = Cs(5,:)*V(:,i)/norm(Cs(5,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz5
        [U,S,V]=svd(evalfr(sys5,RHPz5(j))); Yz5(:,j)=U(:,end); %zero direction
    end
    SGdmin = zeros(np5,nz5);
    for i=1:np5
        for j=1:nz5
            num5=abs((ctranspose(Yz5(:,j))*(evalfr(Gdms51,RHPz5(j))))*((evalfr(Gdms51,RHPp5(i)))^-1)*Yp5(:,i));
            den5=norm(ctranspose(Yz5(:,j))*(evalfr(Gdms51,RHPz5(j))),2)*norm(((evalfr(Gdms51,RHPp5(i)))^-1)*Yp5(:,i),2);
            cos2=(num5/den5)^2;
            sin2=1-cos2;
            
            term5 = (abs(RHPz5(j)+conj(RHPp5(i))))^2/(abs(RHPz5(j)-RHPp5(i)))^2;
            SGdmin(i,j)=norm((ctranspose(Yz5(:,j))*(evalfr(Gdms51,RHPz5(j)))),2)*sqrt(sin2+(term5*cos2));
        end
    end
    SGdmin51 = max(max((SGdmin)));
else
    SGdmin51=0;
end
%% SGd2 for P_rb
[p, z5]=pzmap(sys5);
RHPp5=p(find(p>0));
RHPz5=z5(find(z5>0));

if length(RHPp5)>0 && length(RHPz5)>0
    np5=length(RHPp5); nz5=length(RHPz5);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp5 =  zeros(l,np5);
    Yz5 = zeros(l,nz5);
    for i=1:np5
        Yp5(:,i) = Cs(5,:)*V(:,i)/norm(Cs(5,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz5
        [U,S,V]=svd(evalfr(sys5,RHPz5(j))); Yz5(:,j)=U(:,end); %zero direction
    end
    SGdmin = zeros(np5,nz5);
    for i=1:np5
        for j=1:nz5
            num5=abs((ctranspose(Yz5(:,j))*(evalfr(Gdms52,RHPz5(j))))*((evalfr(Gdms52,RHPp5(i)))^-1)*Yp5(:,i));
            den5=norm(ctranspose(Yz5(:,j))*(evalfr(Gdms52,RHPz5(j))),2)*norm(((evalfr(Gdms52,RHPp5(i)))^-1)*Yp5(:,i),2);
            cos2=(num5/den5)^2;
            sin2=1-cos2;
            
            term5 = (abs(RHPz5(j)+conj(RHPp5(i))))^2/(abs(RHPz5(j)-RHPp5(i)))^2;
            SGdmin(i,j)=norm((ctranspose(Yz5(:,j))*(evalfr(Gdms52,RHPz5(j)))),2)*sqrt(sin2+(term5*cos2));
        end
    end
    SGdmin52 = max(max((SGdmin)));
else
    SGdmin52=0;
end
%% SGd1 for DP_r
[p, z6]=pzmap(sys6);
RHPp6=p(find(p>0));
RHPz6=z6(find(z6>0));

if length(RHPp6)>0 && length(RHPz6)>0
    np6=length(RHPp6); nz6=length(RHPz6);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp6 =  zeros(l,np6);
    Yz6 = zeros(l,nz6);
    for i=1:np6
        Yp6(:,i) = Cs(6,:)*V(:,i)/norm(Cs(6,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz6
        [U,S,V]=svd(evalfr(sys6,RHPz6(j))); Yz6(:,j)=U(:,end); %zero direction
    end
    SGdmin = zeros(np6,nz6);
    for i=1:np6
        for j=1:nz6
            num6=abs((ctranspose(Yz6(:,j))*(evalfr(Gdms61,RHPz6(j))))*((evalfr(Gdms61,RHPp6(i)))^-1)*Yp6(:,i));
            den6=norm(ctranspose(Yz6(:,j))*(evalfr(Gdms61,RHPz6(j))),2)*norm(((evalfr(Gdms61,RHPp6(i)))^-1)*Yp6(:,i),2);
            cos2=(num6/den6)^2;
            sin2=1-cos2;
            
            term6 = (abs(RHPz6(j)+conj(RHPp6(i))))^2/(abs(RHPz6(j)-RHPp6(i)))^2;
            SGdmin(i,j)=norm((ctranspose(Yz6(:,j))*(evalfr(Gdms61,RHPz6(j)))),2)*sqrt(sin2+(term6*cos2));
        end
    end
    SGdmin61 = max(max((SGdmin)));
else
    SGdmin61=0;
end
%% SGd2 for DP_r
[p, z6]=pzmap(sys6);
RHPp6=p(find(p>0));
RHPz6=z6(find(z6>0));

if length(RHPp6)>0 && length(RHPz6)>0
    np6=length(RHPp6); nz6=length(RHPz6);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp6 =  zeros(l,np6);
    Yz6 = zeros(l,nz6);
    for i=1:np6
        Yp6(:,i) = Cs(6,:)*V(:,i)/norm(Cs(6,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz6
        [U,S,V]=svd(evalfr(sys6,RHPz6(j))); Yz6(:,j)=U(:,end); %zero direction
    end
    SGdmin = zeros(np6,nz6);
    for i=1:np6
        for j=1:nz6
            num6=abs((ctranspose(Yz6(:,j))*(evalfr(Gdms62,RHPz6(j))))*((evalfr(Gdms62,RHPp6(i)))^-1)*Yp6(:,i));
            den6=norm(ctranspose(Yz6(:,j))*(evalfr(Gdms62,RHPz6(j))),2)*norm(((evalfr(Gdms62,RHPp6(i)))^-1)*Yp6(:,i),2);
            cos2=(num6/den6)^2;
            sin2=1-cos2;
            
            term6 = (abs(RHPz6(j)+conj(RHPp6(i))))^2/(abs(RHPz6(j)-RHPp6(i)))^2;
            SGdmin(i,j)=norm((ctranspose(Yz6(:,j))*(evalfr(Gdms62,RHPz6(j)))),2)*sqrt(sin2+(term6*cos2));
        end
    end
    SGdmin62 = max(max((SGdmin)));
else
    SGdmin62=0;
end

%% SGd1 for P_t
[p, z7]=pzmap(sys7);
RHPp7=p(find(p>0));
RHPz7=z7(find(z7>0));

if length(RHPp7)>0 && length(RHPz7)>0
    np7=length(RHPp7); nz7=length(RHPz7);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp7 =  zeros(l,np7);
    Yz7 = zeros(l,nz7);
    for i=1:np7
        Yp7(:,i) = Cs(7,:)*V(:,i)/norm(Cs(7,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz7
        [U,S,V]=svd(evalfr(sys7,RHPz7(j))); Yz7(:,j)=U(:,end); %zero direction
    end
    SGdmin = zeros(np7,nz7);
    for i=1:np7
        for j=1:nz7
            num7=abs((ctranspose(Yz7(:,j))*(evalfr(Gdms71,RHPz7(j))))*((evalfr(Gdms71,RHPp7(i)))^-1)*Yp7(:,i));
            den7=norm(ctranspose(Yz7(:,j))*(evalfr(Gdms71,RHPz7(j))),2)*norm(((evalfr(Gdms71,RHPp7(i)))^-1)*Yp7(:,i),2);
            cos2=(num7/den7)^2;
            sin2=1-cos2;
            
            term7 = (abs(RHPz7(j)+conj(RHPp7(i))))^2/(abs(RHPz7(j)-RHPp7(i)))^2;
            SGdmin(i,j)=norm((ctranspose(Yz7(:,j))*(evalfr(Gdms71,RHPz7(j)))),2)*sqrt(sin2+(term7*cos2));
        end
    end
    SGdmin71 = max(max((SGdmin)));
else
    SGdmin71=0;
end

%% SGd2 for P_t
[p, z7]=pzmap(sys7);
RHPp7=p(find(p>0));
RHPz7=z7(find(z7>0));

if length(RHPp7)>0 && length(RHPz7)>0
    np7=length(RHPp7); nz7=length(RHPz7);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp7 =  zeros(l,np7);
    Yz7 = zeros(l,nz7);
    for i=1:np7
        Yp7(:,i) = Cs(7,:)*V(:,i)/norm(Cs(7,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz7
        [U,S,V]=svd(evalfr(sys7,RHPz7(j))); Yz7(:,j)=U(:,end); %zero direction
    end
    SGdmin = zeros(np7,nz7);
    for i=1:np7
        for j=1:nz7
            num7=abs((ctranspose(Yz7(:,j))*(evalfr(Gdms72,RHPz7(j))))*((evalfr(Gdms72,RHPp7(i)))^-1)*Yp7(:,i));
            den7=norm(ctranspose(Yz7(:,j))*(evalfr(Gdms72,RHPz7(j))),2)*norm(((evalfr(Gdms72,RHPp7(i)))^-1)*Yp7(:,i),2);
            cos2=(num7/den7)^2;
            sin2=1-cos2;
            
            term7 = (abs(RHPz7(j)+conj(RHPp7(i))))^2/(abs(RHPz7(j)-RHPp7(i)))^2;
            SGdmin(i,j)=norm((ctranspose(Yz7(:,j))*(evalfr(Gdms72,RHPz7(j)))),2)*sqrt(sin2+(term7*cos2));
        end
    end
    SGdmin72 = max(max((SGdmin)));
else
    SGdmin72=0;
end

%% SGd1 for Q_out
[p, z8]=pzmap(sys8);
RHPp8=p(find(p>0));
RHPz8=z8(find(z8>0));

if length(RHPp8)>0 && length(RHPz8)>0
    np8=length(RHPp8); nz8=length(RHPz8);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp8 =  zeros(l,np8);
    Yz8 = zeros(l,nz8);
    for i=1:np8
        Yp8(:,i) = Cs(8,:)*V(:,i)/norm(Cs(8,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz8
        [U,S,V]=svd(evalfr(sys8,RHPz8(j))); Yz8(:,j)=U(:,end); %zero direction
    end
    SGdmin = zeros(np8,nz8);
    for i=1:np8
        for j=1:nz8
            num8=abs((ctranspose(Yz8(:,j))*(evalfr(Gdms81,RHPz8(j))))*((evalfr(Gdms81,RHPp8(i)))^-1)*Yp8(:,i));
            den8=norm(ctranspose(Yz8(:,j))*(evalfr(Gdms81,RHPz8(j))),2)*norm(((evalfr(Gdms81,RHPp8(i)))^-1)*Yp8(:,i),2);
            cos2=(num8/den8)^2;
            sin2=1-cos2;
            
            term8 = (abs(RHPz8(j)+conj(RHPp8(i))))^2/(abs(RHPz8(j)-RHPp8(i)))^2;
            SGdmin(i,j)=norm((ctranspose(Yz8(:,j))*(evalfr(Gdms81,RHPz8(j)))),2)*sqrt(sin2+(term8*cos2));
        end
    end
    SGdmin81 = max(max((SGdmin)));
else
    SGdmin81=0;
end

%% SGd2 for Q_out
[p, z8]=pzmap(sys8);
RHPp8=p(find(p>0));
RHPz8=z8(find(z8>0));

if length(RHPp8)>0 && length(RHPz8)>0
    np8=length(RHPp8); nz8=length(RHPz8);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp8 =  zeros(l,np8);
    Yz8 = zeros(l,nz8);
    for i=1:np8
        Yp8(:,i) = Cs(8,:)*V(:,i)/norm(Cs(8,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz8
        [U,S,V]=svd(evalfr(sys8,RHPz8(j))); Yz8(:,j)=U(:,end); %zero direction
    end
    SGdmin = zeros(np8,nz8);
    for i=1:np8
        for j=1:nz8
            num8=abs((ctranspose(Yz8(:,j))*(evalfr(Gdms82,RHPz8(j))))*((evalfr(Gdms82,RHPp8(i)))^-1)*Yp8(:,i));
            den8=norm(ctranspose(Yz8(:,j))*(evalfr(Gdms82,RHPz8(j))),2)*norm(((evalfr(Gdms82,RHPp8(i)))^-1)*Yp8(:,i),2);
            cos2=(num8/den8)^2;
            sin2=1-cos2;
            
            term8 = (abs(RHPz8(j)+conj(RHPp8(i))))^2/(abs(RHPz8(j)-RHPp8(i)))^2;
            SGdmin(i,j)=norm((ctranspose(Yz8(:,j))*(evalfr(Gdms82,RHPz8(j)))),2)*sqrt(sin2+(term8*cos2));
        end
    end
    SGdmin82 = max(max((SGdmin)));
else
    SGdmin82=0;
end
%% SGd1 for W_out
[p, z9]=pzmap(sys9);
RHPp9=p(find(p>0));
RHPz9=z9(find(z9>0));

if length(RHPp9)>0 && length(RHPz9)>0
    np9=length(RHPp9); nz9=length(RHPz9);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp9 =  zeros(l,np9);
    Yz9 = zeros(l,nz9);
    for i=1:np9
        Yp9(:,i) = Cs(9,:)*V(:,i)/norm(Cs(9,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz9
        [U,S,V]=svd(evalfr(sys9,RHPz9(j))); Yz9(:,j)=U(:,end); %zero direction
    end
    SGdmin = zeros(np9,nz9);
    for i=1:np9
        for j=1:nz9
            num9=abs((ctranspose(Yz9(:,j))*(evalfr(Gdms91,RHPz9(j))))*((evalfr(Gdms91,RHPp9(i)))^-1)*Yp9(:,i));
            den9=norm(ctranspose(Yz9(:,j))*(evalfr(Gdms91,RHPz9(j))),2)*norm(((evalfr(Gdms91,RHPp9(i)))^-1)*Yp9(:,i),2);
            cos2=(num9/den9)^2;
            sin2=1-cos2;
            
            term9 = (abs(RHPz9(j)+conj(RHPp9(i))))^2/(abs(RHPz9(j)-RHPp9(i)))^2;
            SGdmin(i,j)=norm((ctranspose(Yz9(:,j))*(evalfr(Gdms91,RHPz9(j)))),2)*sqrt(sin2+(term9*cos2));
        end
    end
    SGdmin91 = max(max((SGdmin)));
else
    SGdmin91=0;
end

%% SGd2 for W_out
[p, z9]=pzmap(sys9);
RHPp9=p(find(p>0));
RHPz9=z9(find(z9>0));

if length(RHPp9)>0 && length(RHPz9)>0
    np9=length(RHPp9); nz9=length(RHPz9);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp9 =  zeros(l,np9);
    Yz9 = zeros(l,nz9);
    for i=1:np9
        Yp9(:,i) = Cs(9,:)*V(:,i)/norm(Cs(9,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz9
        [U,S,V]=svd(evalfr(sys9,RHPz9(j))); Yz9(:,j)=U(:,end); %zero direction
    end
    SGdmin = zeros(np9,nz9);
    for i=1:np9
        for j=1:nz9
            num9=abs((ctranspose(Yz9(:,j))*(evalfr(Gdms92,RHPz9(j))))*((evalfr(Gdms92,RHPp9(i)))^-1)*Yp9(:,i));
            den9=norm(ctranspose(Yz9(:,j))*(evalfr(Gdms92,RHPz9(j))),2)*norm(((evalfr(Gdms92,RHPp9(i)))^-1)*Yp9(:,i),2);
            cos2=(num9/den9)^2;
            sin2=1-cos2;
            
            term9 = (abs(RHPz9(j)+conj(RHPp9(i))))^2/(abs(RHPz9(j)-RHPp9(i)))^2;
            SGdmin(i,j)=norm((ctranspose(Yz9(:,j))*(evalfr(Gdms92,RHPz9(j)))),2)*sqrt(sin2+(term9*cos2));
        end
    end
    SGdmin92 = max(max((SGdmin)));
else
    SGdmin92=0;
end

%% SGd1 for Rho_t
[p, z10]=pzmap(sys10);
RHPp10=p(find(p>0));
RHPz10=z10(find(z10>0));

if length(RHPp10)>0 && length(RHPz10)>0
    np10=length(RHPp10); nz10=length(RHPz10);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp10 =  zeros(l,np10);
    Yz10 = zeros(l,nz10);
    for i=1:np10
        Yp10(:,i) = Cs(10,:)*V(:,i)/norm(Cs(10,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz10
        [U,S,V]=svd(evalfr(sys10,RHPz10(j))); Yz10(:,j)=U(:,end); %zero direction
    end
    SGdmin = zeros(np10,nz10);
    for i=1:np10
        for j=1:nz10
            num10=abs((ctranspose(Yz10(:,j))*(evalfr(Gdms101,RHPz10(j))))*((evalfr(Gdms101,RHPp10(i)))^-1)*Yp10(:,i));
            den10=norm(ctranspose(Yz10(:,j))*(evalfr(Gdms101,RHPz10(j))),2)*norm(((evalfr(Gdms101,RHPp10(i)))^-1)*Yp10(:,i),2);
            cos2=(num10/den10)^2;
            sin2=1-cos2;
            
            term10 = (abs(RHPz10(j)+conj(RHPp10(i))))^2/(abs(RHPz10(j)-RHPp10(i)))^2;
            SGdmin(i,j)=norm((ctranspose(Yz10(:,j))*(evalfr(Gdms101,RHPz10(j)))),2)*sqrt(sin2+(term10*cos2));
        end
    end
    SGdmin101 = max(max((SGdmin)));
else
    SGdmin101=0;
end

%% SGd2 for Rho_t
[p, z10]=pzmap(sys10);
RHPp10=p(find(p>0));
RHPz10=z10(find(z10>0));

if length(RHPp10)>0 && length(RHPz10)>0
    np10=length(RHPp10); nz10=length(RHPz10);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp10 =  zeros(l,np10);
    Yz10 = zeros(l,nz10);
    for i=1:np10
        Yp10(:,i) = Cs(10,:)*V(:,i)/norm(Cs(10,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz10
        [U,S,V]=svd(evalfr(sys10,RHPz10(j))); Yz10(:,j)=U(:,end); %zero direction
    end
    SGdmin = zeros(np10,nz10);
    for i=1:np10
        for j=1:nz10
            num10=abs((ctranspose(Yz10(:,j))*(evalfr(Gdms102,RHPz10(j))))*((evalfr(Gdms102,RHPp10(i)))^-1)*Yp10(:,i));
            den10=norm(ctranspose(Yz10(:,j))*(evalfr(Gdms102,RHPz10(j))),2)*norm(((evalfr(Gdms102,RHPp10(i)))^-1)*Yp10(:,i),2);
            cos2=(num10/den10)^2;
            sin2=1-cos2;
            
            term10 = (abs(RHPz10(j)+conj(RHPp10(i))))^2/(abs(RHPz10(j)-RHPp10(i)))^2;
            SGdmin(i,j)=norm((ctranspose(Yz10(:,j))*(evalfr(Gdms102,RHPz10(j)))),2)*sqrt(sin2+(term10*cos2));
        end
    end
    SGdmin102 = max(max((SGdmin)));
else
    SGdmin102=0;
end
%% SGd1 for Alpha_L
[p, z11]=pzmap(sys11);
RHPp11=p(find(p>0));
RHPz11=z11(find(z11>0));

if length(RHPp11)>0 && length(RHPz11)>0
    np11=length(RHPp11); nz11=length(RHPz11);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp11 =  zeros(l,np11);
    Yz11 = zeros(l,nz11);
    for i=1:np11
        Yp11(:,i) = Cs(11,:)*V(:,i)/norm(Cs(11,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz11
        [U,S,V]=svd(evalfr(sys11,RHPz11(j))); Yz11(:,j)=U(:,end); %zero direction
    end
    SGdmin = zeros(np11,nz11);
    for i=1:np11
        for j=1:nz11
            num11=abs((ctranspose(Yz11(:,j))*(evalfr(Gdms111,RHPz11(j))))*((evalfr(Gdms111,RHPp11(i)))^-1)*Yp11(:,i));
            den11=norm(ctranspose(Yz11(:,j))*(evalfr(Gdms111,RHPz11(j))),2)*norm(((evalfr(Gdms111,RHPp11(i)))^-1)*Yp11(:,i),2);
            cos2=(num11/den11)^2;
            sin2=1-cos2;
            
            term11 = (abs(RHPz11(j)+conj(RHPp11(i))))^2/(abs(RHPz11(j)-RHPp11(i)))^2;
            SGdmin(i,j)=norm((ctranspose(Yz11(:,j))*(evalfr(Gdms111,RHPz11(j)))),2)*sqrt(sin2+(term11*cos2));
        end
    end
    SGdmin111 = max(max((SGdmin)));
else
    SGdmin111=0;
end

%% SGd2 for Alpha_L
[p, z11]=pzmap(sys11);
RHPp11=p(find(p>0));
RHPz11=z11(find(z11>0));

if length(RHPp11)>0 && length(RHPz11)>0
    np11=length(RHPp11); nz11=length(RHPz11);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp11 =  zeros(l,np11);
    Yz11 = zeros(l,nz11);
    for i=1:np11
        Yp11(:,i) = Cs(11,:)*V(:,i)/norm(Cs(11,:)*V(:,i)); %Pole direction
    end
    
    for j=1:nz11
        [U,S,V]=svd(evalfr(sys11,RHPz11(j))); Yz11(:,j)=U(:,end); %zero direction
    end
    SGdmin = zeros(np11,nz11);
    for i=1:np11
        for j=1:nz11
            num11=abs((ctranspose(Yz11(:,j))*(evalfr(Gdms112,RHPz11(j))))*((evalfr(Gdms112,RHPp11(i)))^-1)*Yp11(:,i));
            den11=norm(ctranspose(Yz11(:,j))*(evalfr(Gdms112,RHPz11(j))),2)*norm(((evalfr(Gdms112,RHPp11(i)))^-1)*Yp11(:,i),2);
            cos2=(num11/den11)^2;
            sin2=1-cos2;
            
            term11 = (abs(RHPz11(j)+conj(RHPp11(i))))^2/(abs(RHPz11(j)-RHPp11(i)))^2;
            SGdmin(i,j)=norm((ctranspose(Yz11(:,j))*(evalfr(Gdms112,RHPz11(j)))),2)*sqrt(sin2+(term11*cos2));
        end
    end
    SGdmin112 = max(max((SGdmin)));
else
    SGdmin112=0;
end
%% Pole Vectors
    np1=length(RHPp1); nz1=length(RHPz1);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp1 = zeros(l,np1);
    Yz1 = zeros(l,nz1);
    for i=1:np1
        Yp1(:,i) = Cs(1,:)*V(:,i)/norm(V(:,i)); %Pole vector
    end

    np2=length(RHPp2); nz2=length(RHPz2);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp2 = zeros(l,np2);
    Yz2 = zeros(l,nz2);
    for i=1:np2
        Yp2(:,i) = Cs(2,:)*V(:,i)/norm(V(:,i)); %Pole vector
    end

    np3=length(RHPp3); nz3=length(RHPz3);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp3 = zeros(l,np3);
    Yz3 = zeros(l,nz3);
    for i=1:np3
        Yp3(:,i) = Cs(3,:)*V(:,i)/norm(V(:,i)); %Pole vector
    end

    np4=length(RHPp4); nz4=length(RHPz4);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp4 = zeros(l,np4);
    Yz4 = zeros(l,nz4);
    for i=1:np4
        Yp4(:,i) = Cs(4,:)*V(:,i)/norm(V(:,i)); %Pole vector
    end

    np5=length(RHPp5); nz5=length(RHPz5);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp5 = zeros(l,np5);
    Yz5 = zeros(l,nz5);
    for i=1:np5
        Yp5(:,i) = Cs(5,:)*V(:,i)/norm(V(:,i)); %Pole vector
    end
    
    np6=length(RHPp6); nz6=length(RHPz6);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp6 = zeros(l,np6);
    Yz6 = zeros(l,nz6);
    for i=1:np6
        Yp6(:,i) = Cs(6,:)*V(:,i)/norm(V(:,i)); %Pole vector
    end
    
    np7=length(RHPp7); nz7=length(RHPz7);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp7 = zeros(l,np7);
    Yz7 = zeros(l,nz7);
    for i=1:np7
        Yp7(:,i) = Cs(7,:)*V(:,i)/norm(V(:,i)); %Pole vector
    end
    
     np8=length(RHPp8); nz8=length(RHPz8);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp8 = zeros(l,np8);
    Yz8 = zeros(l,nz8);
    for i=1:np8
        Yp8(:,i) = Cs(8,:)*V(:,i)/norm(V(:,i)); %Pole vector
    end
    
    np9=length(RHPp9); nz9=length(RHPz9);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp9 = zeros(l,np9);
    Yz9 = zeros(l,nz9);
    for i=1:np9
        Yp9(:,i) = Cs(9,:)*V(:,i)/norm(V(:,i)); %Pole vector
    end
    
    np10=length(RHPp10); nz10=length(RHPz10);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp10 = zeros(l,np10);
    Yz10 = zeros(l,nz10);
    for i=1:np10
        Yp10(:,i) = Cs(10,:)*V(:,i)/norm(V(:,i)); %Pole vector
    end
    
    np11=length(RHPp11); nz11=length(RHPz11);
    [V,E]=eig(A); V=V(:,find(diag(E)>0));
    Yp11 = zeros(l,np11);
    Yz11 = zeros(l,nz11);
    for i=1:np11
        Yp11(:,i) = Cs(11,:)*V(:,i)/norm(V(:,i)); %Pole vector
    end


%%
disp('Steady-state gains');
disp(['The steady-state gain of G1 (P_bh) is: ', num2str(freqresp(G1,0))]);
disp(['The steady-state gain of G2 (P_wh) is: ', num2str(freqresp(G2,0))]);
disp(['The steady-state gain of G3 (W_in) is: ', num2str(freqresp(G3,0))]);
disp(['The steady-state gain of G4 (P_in) is: ', num2str(freqresp(G4,0))]);
disp(['The steady-state gain of G5 (P_rb) is: ', num2str(freqresp(G5,0))]);
disp(['The steady-state gain of G6 (DP_r) is: ', num2str(freqresp(G6,0))]);
disp(['The steady-state gain of G7 (P_t) is: ', num2str(freqresp(G7,0))]);
disp(['The steady-state gain of G8 (Q_out) is: ', num2str(freqresp(G8,0))]);
disp(['The steady-state gain of G9 (W_out) is: ', num2str(freqresp(G9,0))]);
disp(['The steady-state gain of G10 (Rho_t) is: ', num2str(freqresp(G10,0))]);
disp(['The steady-state gain of G11 (Alpha_L) is: ', num2str(freqresp(G11,0))]);
disp('##########################################################');
disp('Pole vectors');
disp(['Pole vectors of (P_bh) is: ', num2str(min(abs(Yp1)))]);
disp(['Pole vectors of (P_wh) is: ', num2str(min(abs(Yp2)))]);
disp(['Pole vectors of (W_in) is: ', num2str(min(abs(Yp3)))]);
disp(['Pole vectors of (P_in) is: ', num2str(min(abs(Yp4)))]);
disp(['Pole vectors of (P_rb) is: ', num2str(min(abs(Yp5)))]);
disp(['Pole vectors of (DP_r) is: ', num2str(min(abs(Yp6)))]);
disp(['Pole vectors of (P_t) is: ', num2str(min(abs(Yp7)))]);
disp(['Pole vectors of (Q_out) is: ', num2str(min(abs(Yp8)))]);
disp(['Pole vectors of (W_out) is: ', num2str(min(abs(Yp9)))]);
disp(['Pole vectors of (Rho_t) is: ', num2str(min(abs(Yp10)))]);
disp(['Pole vectors of (Alpha_L) is: ', num2str(min(abs(Yp11)))]);
disp('##########################################################');
disp('Bounds on S and T');
disp(['The lowest achivable peak for Msmin1 (P_bh) is: ', num2str(Msmin1)]);
disp(['The lowest achivable peak for Msmin2 (P_wh) is: ', num2str(Msmin2)]);
disp(['The lowest achivable peak for Msmin3 (W_in) is: ', num2str(Msmin3)]);
disp(['The lowest achivable peak for Msmin4 (P_in) is: ', num2str(Msmin4)]);
disp(['The lowest achivable peak for Msmin5 (P_rb) is: ', num2str(Msmin5)]);
disp(['The lowest achivable peak for Msmin6 (DP_r) is: ', num2str(Msmin6)]);
disp(['The lowest achivable peak for Msmin7 (P_t) is: ', num2str(Msmin7)]);
disp(['The lowest achivable peak for Msmin8 (Q_out) is: ', num2str(Msmin8)]);
disp(['The lowest achivable peak for Msmin9 (W_out) is: ', num2str(Msmin9)]);
disp(['The lowest achivable peak for Msmin8 (Rho_t) is: ', num2str(Msmin10)]);
disp(['The lowest achivable peak for Msmin9 (Alpha_L) is: ', num2str(Msmin11)]);
disp('##########################################################');
disp('Bounds on KS');
disp(['The bound for KSmin1 (P_bh) is: ', num2str(KSmin1)]);
disp(['The bound for KSmin2 (P_wh) is: ', num2str(KSmin2)]);
disp(['The bound for KSmin3 (W_in) is: ', num2str(KSmin3)]);
disp(['The bound for KSmin4 (P_in) is: ', num2str(KSmin4)]);
disp(['The bound for KSmin5 (P_rb) is: ', num2str(KSmin5)]);
disp(['The bound for KSmin6 (DP_r) is: ', num2str(KSmin6)]);
disp(['The bound for KSmin7 (P_t) is: ', num2str(KSmin7)]);
disp(['The bound for KSmin8 (Q_out) is: ', num2str(KSmin8)]);
disp(['The bound for KSmin9 (W_out) is: ', num2str(KSmin9)]);
disp(['The bound for KSmin10 (Rho_t) is: ', num2str(KSmin10)]);
disp(['The bound for KSmin11 (Alpha_L) is: ', num2str(KSmin11)]);
disp('##########################################################');
disp('Bounds on KSGd1');
disp(['The bound for KSGd11 (P_bh) is: ', num2str(KSGd11)]);
disp(['The bound for KSGd21 (P_wh) is: ', num2str(KSGd21)]);
disp(['The bound for KSGd31 (W_in) is: ', num2str(KSGd31)]);
disp(['The bound for KSGd41 (P_in) is: ', num2str(KSGd41)]);
disp(['The bound for KSGd51 (P_rb) is: ', num2str(KSGd51)]);
disp(['The bound for KSGd61 (DP_r) is: ', num2str(KSGd61)]);
disp(['The bound for KSGd71 (P_t) is: ', num2str(KSGd71)]);
disp(['The bound for KSGd81 (Q_out) is: ', num2str(KSGd81)]);
disp(['The bound for KSGd91 (W_out) is: ', num2str(KSGd91)]);
disp(['The bound for KSGd101 (Rho_t) is: ', num2str(KSGd101)]);
disp(['The bound for KSGd111 (Alpha_L) is: ', num2str(KSGd111)]);
disp('Bounds on KSGd2');
disp(['The bound for KSGd12 (P_bh) is: ', num2str(KSGd12)]);
disp(['The bound for KSGd22 (P_wh) is: ', num2str(KSGd22)]);
disp(['The bound for KSGd32 (W_in) is: ', num2str(KSGd32)]);
disp(['The bound for KSGd42 (P_in) is: ', num2str(KSGd42)]);
disp(['The bound for KSGd52 (P_rb) is: ', num2str(KSGd52)]);
disp(['The bound for KSGd62 (DP_r) is: ', num2str(KSGd62)]);
disp(['The bound for KSGd72 (P_t) is: ', num2str(KSGd72)]);
disp(['The bound for KSGd82 (Q_out) is: ', num2str(KSGd82)]);
disp(['The bound for KSGd92 (W_out) is: ', num2str(KSGd92)]);
disp(['The bound for KSGd102 (Rho_t) is: ', num2str(KSGd102)]);
disp(['The bound for KSGd112 (Alpha_L) is: ', num2str(KSGd112)]);
disp('##########################################################');
disp('Bounds on SG');
disp(['The bound for SGmin1 (P_bh) is: ', num2str(SGmin1)]);
disp(['The bound for SGmin2 (P_wh) is: ', num2str(SGmin2)]);
disp(['The bound for SGmin3 (W_in) is: ', num2str(SGmin3)]);
disp(['The bound for SGmin4 (P_in) is: ', num2str(SGmin4)]);
disp(['The bound for SGmin5 (P_rb) is: ', num2str(SGmin5)]);
disp(['The bound for SGmin6 (DP_r) is: ', num2str(SGmin6)]);
disp(['The bound for SGmin7 (P_t) is: ', num2str(SGmin7)]);
disp(['The bound for SGmin8 (Q_out) is: ', num2str(SGmin8)]);
disp(['The bound for SGmin9 (W_out) is: ', num2str(SGmin9)]);
disp(['The bound for SGmin10 (Rho_t) is: ', num2str(SGmin10)]);
disp(['The bound for SGmin11 (Alpha_L) is: ', num2str(SGmin11)]);
disp('##########################################################');
disp('Bounds on SGd1');
disp(['The bound for SGdmin11 (P_bh) is: ', num2str(SGdmin11)]);
disp(['The bound for SGdmin21 (P_wh) is: ', num2str(SGdmin21)]);
disp(['The bound for SGdmin31 (W_in) is: ', num2str(SGdmin31)]);
disp(['The bound for SGdmin41 (P_in) is: ', num2str(SGdmin41)]);
disp(['The bound for SGdmin51 (P_rb) is: ', num2str(SGdmin51)]);
disp(['The bound for SGdmin61 (DP_r) is: ', num2str(SGdmin61)]);
disp(['The bound for SGdmin71 (P_t) is: ', num2str(SGdmin71)]);
disp(['The bound for SGdmin81 (Q_out) is: ', num2str(SGdmin81)]);
disp(['The bound for SGdmin91 (W_out) is: ', num2str(SGdmin91)]);
disp(['The bound for SdGmin101 (Rho_t) is: ', num2str(SGdmin101)]);
disp(['The bound for SdGmin111 (Alpha_L) is: ', num2str(SGdmin111)]);
disp('Bounds on SGd2');
disp(['The bound for SGdmin12 (P_bh) is: ', num2str(SGdmin12)]);
disp(['The bound for SGdmin22 (P_wh) is: ', num2str(SGdmin22)]);
disp(['The bound for SGdmin32 (W_in) is: ', num2str(SGdmin32)]);
disp(['The bound for SGdmin42 (P_in) is: ', num2str(SGdmin42)]);
disp(['The bound for SGdmin52 (P_rb) is: ', num2str(SGdmin52)]);
disp(['The bound for SGdmin62 (DP_r) is: ', num2str(SGdmin62)]);
disp(['The bound for SGdmin72 (P_t) is: ', num2str(SGdmin72)]);
disp(['The bound for SGdmin82 (Q_out) is: ', num2str(SGdmin82)]);
disp(['The bound for SGdmin92 (W_out) is: ', num2str(SGdmin92)]);
disp(['The bound for SdGmin102 (Rho_t) is: ', num2str(SGdmin102)]);
disp(['The bound for SdGmin112 (Alpha_L) is: ', num2str(SGdmin112)]);