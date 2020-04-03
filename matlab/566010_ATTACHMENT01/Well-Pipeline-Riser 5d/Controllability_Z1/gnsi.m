function [Gsi,Gnsi]=gnsi(Gi)

ptot=pole(Gi);
p=ptot(find(real(ptot>0)));
np=length(p);

s=tf('s');
b=1;
c=1;
for i = 1:np
    b=b/(s-p(i));
    c=c*(s-p(i));
end
Gsi = minreal(c*Gi);
Gnsi=abs(freqresp(Gsi,p(1)))*b;