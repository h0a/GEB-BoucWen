
% Developed by Dr. Bruno A. Roccia
% Geophysical Institute, UiB
% Contact: Bruno.roccia@uib.no

%%
clear all
close all
clc

DATA.fe     = 1;
DATA.W      = pi;
DATA.Cycles = 4;
DATA.m     = 1;
DATA.g     = 0;
DATA.K     = 3.33;
DATA.L0    = 1;
DATA.C1    = 0.4;
DATA.B1    = 1;
DATA.B2    = 0.5;
DATA.B3    = 0.5;
DATA.B4    = 1/0.12;
DATA.Alpha = 0.13;
DATA.Rho   = 100;
DATA.n     = 2;
DATA.h     = 0.01;                                                          % Time step
DATA.TOL   = 1e-8;                                                          % Tolerance for all the three loops (Inner and outer)
DATA.Iter  = 100;                                                           % Max number of iterations for NR solver
DATA.Nstep = 1/(DATA.W/(2*pi)) * DATA.Cycles / DATA.h;                      % Number of time steps
DATA.Omega = 0.1;                                                           % Relaxation parameter

DATA.CI    = [0; pi/2];
DATA.Z0    = 0;

%% 

q(:,1) = DATA.CI;
Z(1)   = DATA.Z0;
f(1)   = DATA.Alpha*DATA.K*q(1,1) + (1-DATA.Alpha)*DATA.C1*Z(1);

Time   = 0:DATA.h:DATA.h*DATA.Nstep;
NStep  = DATA.Nstep;

fprintf ('-------------------------------------------------------------\n')
fprintf ('\t Bouc-Wen hysteretic model\n')
fprintf ('-------------------------------------------------------------\n')
fprintf ('Time Step h = %8.6f \n', DATA.h)
fprintf ('Number of Time steps = %i \n', NStep)
fprintf ('Initial Conditions [Eta, Theta, z] = %8.6f %8.6f %8.6f \n', DATA.CI.')
fprintf ('-------------------------------------------------------------\n')
fprintf ('\n')

Aux1 = 0; Aux2 =0;

for i = 1:NStep
    
    %fei = i * DATA.fe/NStep;
    %gi  = i * DATA.g/NStep;
    
    fei  = DATA.fe * sin(DATA.W * (Time(i+1)+Time(i))/2); 
    gi   = DATA.g;
        
    q1 = q(:,i);
    q2 = q1;
    z1 = Z(i);
    z2 = z1;
    
    kk = 1;
    ErrorTotal = 1;
    
    while ErrorTotal >= DATA.TOL && kk <= DATA.Iter
        
        Error = 1;
        k     = 1;
        q_aux = q2;
        while (Error >= DATA.TOL && k <= DATA.Iter)                         % First Newton-Raphson
            
            [JJ] = Jacobian_1 (DATA, q2, q1, z1, z2, fei, gi);
            gg   = Residuum_1 (DATA, q2, q1, z1, z2, fei, gi);
            Dq   = - JJ \ gg;
            
            q2   = q2 + Dq;
            
            Error = norm (Dq, inf);
            
            k = k + 1;
        end
        
        Error = 1;
        
        kz = 1;
        z_aux = z2;
        while (Error >= DATA.TOL && kz <= DATA.Iter)                        % Second Newton-Raphson

            [JJ] = Jacobian_2 (DATA, q2, q1, z1, z2);
            gg   = Residuum_2 (DATA, q2, q1, z1, z2);
            Dz   = - gg / JJ;            
            
            z2   = z2 + Dz;
            
            Error = norm (Dz, inf);
            
            kz = kz + 1;            
            
        end
        
        z2 = (1 - DATA.Omega)*z_aux + DATA.Omega * z2;                      % Relaxation step
        
        ErrorTotal = max(norm(q2-q_aux,inf), abs(z2-z_aux));                % Outer loop Error
        kk = kk + 1;
        
    end
    
    q(:,i+1) = q2;
    Z(i+1)   = z2;
    
    RES.Iter_E(i) = k;
    RES.Iter_Z(i) = kz;
    RES.Iter_Total(i) = kk;
    RES.ErrorT(i)     = ErrorTotal;
    RES.Error_E(i)    = norm(Dq);
    RES.Error_Z(i)    = abs(Dz);
    
    f(i+1)   = DATA.Alpha*DATA.K*q(1,i+1) + (1-DATA.Alpha)*DATA.C1*Z(i+1);

    
    fprintf(repmat('\b', 1, Aux1));
    fprintf(repmat('\b', 1, Aux2));
    msg1 = sprintf('Time Step = %i\n', i);
    msg2 = sprintf('Iterations = %i\n', kk);
    Aux1 = fprintf('%s', msg1);
    Aux2 = fprintf('%s', msg2);
    
end

fprintf ('Max number of iteration through the simulation - Outer Loop = %i \n', max(RES.Iter_Total))
fprintf ('Max ERROR through the simulation - Outer Loop = %i \n', max(RES.ErrorT))

%% PLOTS

subplot(1,2,1)

plot (Time, q(1,:), 'color', 'b', 'linewidth', 2)
hold on
plot (Time, q(2,:), 'color', 'r', 'linewidth', 2)
legend ('Eta', 'Theta')

subplot (1,2,2)

plot (q(1,:), f, 'color', 'b', 'linewidth', 2)
hold on
plot ([0 0], ylim, 'color', 'k')
plot (xlim, [0 0], 'color', 'k')

figure (2)

plot (Time, Z, 'color', 'r', 'linewidth', 2)
title ('z(t)')
xlabel('Time')

