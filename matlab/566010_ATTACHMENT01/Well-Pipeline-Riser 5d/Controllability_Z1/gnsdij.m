function Gns = gnsdij(Gdij)

ztot=zero(Gdij); ptot=pole(Gdij);
z=ztot(find(ztot>0)); p=ptot(find(ptot>0));
nz=length(z); np=length(p);

s=tf('s');
b=1;a=1;     
for i = 1:np
    b=b/(s-p(i));
end;
for i = 1:nz
    a=a*(s-z(i));
end    

Gns = a*b;