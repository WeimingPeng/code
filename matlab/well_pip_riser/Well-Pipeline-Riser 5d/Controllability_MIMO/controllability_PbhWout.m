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
KSGd191_min = max(KSGd_min);
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
    GammaSG19_min = sqrt(max(eig(temp19)));
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
    GammaSGd191_min = sqrt(max(eig(temp19)));
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
    GammaSGd192_min = sqrt(max(eig(temp19)));
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