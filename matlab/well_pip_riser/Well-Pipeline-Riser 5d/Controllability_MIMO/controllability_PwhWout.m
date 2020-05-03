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
KSGd291_min = max(KSGd_min);
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
    GammaSG29_min = sqrt(max(eig(temp29)));
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
    GammaSGd291_min = sqrt(max(eig(temp29)));
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
    GammaSGd292_min = sqrt(max(eig(temp29)));
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