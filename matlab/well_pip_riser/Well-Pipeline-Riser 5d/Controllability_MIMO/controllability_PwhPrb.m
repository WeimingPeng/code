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
    GammaSG25_min = sqrt(max(eig(temp25)));
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
    GammaSGd251_min = sqrt(max(eig(temp25)));
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