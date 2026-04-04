function [JJ] = Jacobian_1 (DATA, q2, q1, z1, z2, fe, g)

x2 = q2(1); t2=q2(2); x1=q1(1); t1=q1(2);

m  = DATA.m;
K  = DATA.K;
L0 = DATA.L0;
C1 = DATA.C1;
B1 = DATA.B1;
B2 = DATA.B2;
B3 = DATA.B3;
B4 = DATA.B4;
alpha = DATA.Alpha;
rho = DATA.Rho;
n  =  DATA.n;
h  = DATA.h;

JJ = zeros(2,2);

DGx = B1 * B4 / h - B2 * B4 * tanh(rho * (z2 + z1) / 0.2e1) ^ n * (z2 / 0.2e1 + z1 / 0.2e1) ^ n / h + B3 * B4 * tanh(rho * (z2 + z1) / 0.2e1) ^ n * tanh(rho * (z2 + z1) * (-x2 + x1) / h / 0.2e1) * (z2 / 0.2e1 + z1 / 0.2e1) ^ n / h - B3 * B4 * (x2 - x1) * tanh(rho * (z2 + z1) / 0.2e1) ^ n * (0.1e1 - tanh(rho * (z2 + z1) * (-x2 + x1) / h / 0.2e1) ^ 2) * rho * (z2 + z1) * (z2 / 0.2e1 + z1 / 0.2e1) ^ n / h ^ 2 / 0.2e1;
DGz = -B2 * B4 * (x2 - x1) * tanh(rho * (z2 + z1) / 0.2e1) ^ (n - 1) * n * (0.1e1 - tanh(rho * (z2 + z1) / 0.2e1) ^ 2) * rho * (z2 / 0.2e1 + z1 / 0.2e1) ^ n / h / 0.2e1 - B2 * B4 * (x2 - x1) * tanh(rho * (z2 + z1) / 0.2e1) ^ n * (z2 / 0.2e1 + z1 / 0.2e1) ^ (n - 1) * n / h / 0.2e1 + B3 * B4 * (x2 - x1) * tanh(rho * (z2 + z1) / 0.2e1) ^ (n - 1) * n * (0.1e1 - tanh(rho * (z2 + z1) / 0.2e1) ^ 2) * rho * tanh(rho * (z2 + z1) * (-x2 + x1) / h / 0.2e1) * (z2 / 0.2e1 + z1 / 0.2e1) ^ n / h / 0.2e1 - B3 * B4 * (-x2 + x1) ^ 2 * tanh(rho * (z2 + z1) / 0.2e1) ^ n * (0.1e1 - tanh(rho * (z2 + z1) * (-x2 + x1) / h / 0.2e1) ^ 2) * rho * (z2 / 0.2e1 + z1 / 0.2e1) ^ n / h ^ 2 / 0.2e1 + B3 * B4 * (x2 - x1) * tanh(rho * (z2 + z1) / 0.2e1) ^ n * tanh(rho * (z2 + z1) * (-x2 + x1) / h / 0.2e1) * (z2 / 0.2e1 + z1 / 0.2e1) ^ (n - 1) * n / h / 0.2e1;

DZDX = DGx / (1/h - DGz);

%JJ(1,1) = -alpha*K/2 - 1/2*(1-alpha)*C1 * DZDX;
JJ(1,1) = -alpha*K/2;
JJ(1,2) = 1/2*fe*cos(0.5*(t2+t1)) - 1/2*m*g*sin(0.5*(t2+t1));
JJ(2,1) = 1/2*fe*cos(t2) - 0.5*m*g*sin(t2);
JJ(2,2) = -1/2*fe*(L0+x2)*sin(t2) - 1/2*m*g*(L0+x2)*cos(t2);




end

