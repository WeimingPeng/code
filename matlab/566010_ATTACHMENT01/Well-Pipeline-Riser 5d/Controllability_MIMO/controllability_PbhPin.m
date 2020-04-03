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