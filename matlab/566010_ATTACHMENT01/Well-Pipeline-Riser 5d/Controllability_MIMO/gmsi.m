function [Gsi Gmi Gmsi]=gmsi(Gi)

ztot=zero(Gi); ptot=pole(Gi);
z=ztot(find(real(ztot>0))); p=ptot(find(real(ptot>0)));
nz=length(z); np=length(p);

s=tf('s');
b=1;a=1;     
for i = 1:np
    b=b*(s-p(i))/(s+p(i));
end;
for i = 1:nz
    a=a*(s+z(i))/(s-z(i));
end    

Gsi=minreal(tf(b*Gi));Gmi=minreal(tf(Gi*a)); Gmsi=minreal(tf(b*Gi*a));