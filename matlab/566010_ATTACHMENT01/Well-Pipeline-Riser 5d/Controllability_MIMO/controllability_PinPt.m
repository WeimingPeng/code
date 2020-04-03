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
    GammaSGd472_min = sqrt(max(eig(temp47)));
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