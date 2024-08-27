function [df,ct]=LinearC(T,c,precision) 

for i=1:length(c)
 ct=0;%mise à zero du compteur
 w=2*pi/T;
 k=w/c(i);
 g=9.81;%gravite
 do=1000; %valeur quelconque grande
 d=c(i)^2/g;
 while(abs(do-d)>precision)
 ct=ct+1;
do=d;
dispe=w^2-g*k*tanh(k*d);
fdispe=-g*(k^2)./(cosh(k*d)^2);
d=d-dispe./fdispe;
 end
df(i)=d;
 
end
end