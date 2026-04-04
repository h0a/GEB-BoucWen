function [gg] = Residuum_1 (DATA, q2, q1, z1, z2, fe, g)

x2 = q2(1); t2=q2(2); x1=q1(1); t1=q1(2);

m  = DATA.m;
K  = DATA.K;
L0 = DATA.L0;
C1 = DATA.C1;
alpha = DATA.Alpha;


gg = zeros(2,1);

gg(1,1) = - (alpha*K*(0.5*(x2+x1)) + (1-alpha)*C1*0.5*(z2+z1)) + fe*sin(0.5*(t2+t1)) + m*g*cos(0.5*(t2+t1));
gg(2,1) = 0.5*fe*((L0+x2)*cos(t2) + (L0+x1)*cos(t1)) - 0.5*m*g*((L0+x2)*sin(t2) + (L0+x1)*sin(t1));




end

