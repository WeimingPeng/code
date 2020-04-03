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
KSGd271_min = max(KSGd_min);
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
    GammaSG27_min = sqrt(max(eig(temp27)));
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
    GammaSGd271_min = sqrt(max(eig(temp27)));
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
    GammaSGd272_min = sqrt(max(eig(temp27)));
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