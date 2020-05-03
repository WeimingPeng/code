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
KSGd171_min = max(KSGd_min);
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
    GammaSG17_min = sqrt(max(eig(temp17)));
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
    GammaSGd171_min = sqrt(max(eig(temp17)));
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
    GammaSGd172_min = sqrt(max(eig(temp17)));
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