function [gg] = Residuum_2 (DATA, q2, q1, z1, z2)

x2 = q2(1); x1=q1(1);


B1 = DATA.B1;
B2 = DATA.B2;
B3 = DATA.B3;
B4 = DATA.B4;
alpha = DATA.Alpha;
rho = DATA.Rho;
n  =  DATA.n;
h  = DATA.h;

gg = (z2 - z1) / h - B1 * B4 * (x2 - x1) / h + B2 * B4 * (x2 - x1) * tanh(rho * (z2 + z1) / 0.2e1) ^ n * (z2 / 0.2e1 + z1 / 0.2e1) ^ n / h + B3 * B4 * (x2 - x1) * tanh(rho * (z2 + z1) / 0.2e1) ^ n * tanh(rho * (z2 + z1) * (x2 - x1) / h / 0.2e1) * (z2 / 0.2e1 + z1 / 0.2e1) ^ n / h;




end

