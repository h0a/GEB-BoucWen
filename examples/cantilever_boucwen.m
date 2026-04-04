clear
close all
clc

% static analysis of GEB examples: a cantilever straight beam subjected to
% harmonic point load at the free end, without gravity, and
% using the Bouc-Wen hysteretic model for 2 bending responses



%% INPUT

beam.L = 4;
beam.EA = 100;
beam.GA1 = 5e10;
beam.GA2 = 5e10;
beam.EI1 = 200;
beam.EI2 = 200;
beam.GIt = 2e10;

% external forces
beam.Fend = [-1 0 0];              % load amplitude
beam.bodyFext = @(t) [0; 0; 0];     % f(t): parametrized body force function
                                    % (if straight beam: t = arc length s)
                                    % if considering gravity, it goes here.

% boundary conditions
beam.BCs = {'clamped', 'free'};     % boundary conditions at the left and right ends - 'clamped', 'fixed', 'free', 'custom'
beam.customFixedDofs = {[], []};    % if 'custom': array of local dof number of fixed dofs at the left and right ends

% discretization
mesh.nelms = 4;

% parameters to estimate number of time (load) steps and load factors = sin(wF * t)
beam.loadOmega = pi;        % wF
beam.Ncycles = 3;           % number of load cycles
beam.timeStep = 0.01;       % time step



% parameters alpha and beta for the stress function sigma_i(e,z) of each stress component
beam.alpha4Sigma = ones(6,1);
beam.alpha4Sigma(5) = 0.13;   % bending component(s)

beam.beta4Sigma = ones(6,1);
beam.beta4Sigma(5) = 0.4;


% initial conditions of strain components and z of each element
beam.e0 = zeros(6*mesh.nelms,1);
beam.z0 = zeros(6*mesh.nelms,1);


% parameters for the Bouc-Wen hysteretic model
beam.BWomega = zeros(6,1);              % relaxation parameter for each z-component (same for each element)
beam.BWomega(5) = 0.1;

beam.BWdata.B1 = 1;
beam.BWdata.B2 = 0.5;
beam.BWdata.B3 = 0.5;
beam.BWdata.B4 = 1/0.12;
beam.BWdata.Alpha = beam.alpha4Sigma(5);
beam.BWdata.Rho = 100;
beam.BWdata.n = 2;
beam.BWdata.h = beam.timeStep;



%% PRE-PROCESSING

[mesh, beam] = preprocessing(mesh, beam);

beam = useBoucWenStresses(beam);


% node vector of the reference configuration
beam.nodalTheta = linspace(0, beam.L, mesh.numNodes);
beam.nodeVec = [zeros(1,mesh.numNodes); zeros(1,mesh.numNodes); beam.nodalTheta];

% nodal directors of the reference configuration (orthonormal)
beam.nodalD1 = [ones(1,mesh.nelms+1); zeros(1,mesh.nelms+1); zeros(1,mesh.nelms+1)];
beam.nodalD3 = [zeros(1,mesh.nelms+1); zeros(1,mesh.nelms+1); ones(1,mesh.nelms+1)];
beam.nodalD2 = [zeros(1,mesh.nelms+1); ones(1,mesh.nelms+1); zeros(1,mesh.nelms+1)];


beam.dsdtJacobianFunc = @(t) 1;     % function of the jacobian ds/dt
                                    % where t is the parametric coorcinates and s the arc length
                                    % (if straight beam: s = t and hence ds/dt = @(t) 1)


% nodal external force vectors (point force at each node)
beam.nodalFext = zeros(3,mesh.numNodes);
beam.nodalFext(:,end) = beam.Fend;

% number of time/load steps
beam.numLoadSteps = 1/(beam.loadOmega/(2*pi)) * beam.Ncycles / beam.timeStep;

% load factors =  sin(wF * t)
beam.T = beam.timeStep * beam.numLoadSteps;
beam.ttVec = 0:beam.timeStep:beam.T;
beam.loadFactors = sin(beam.loadOmega .* beam.ttVec(2:end));



%% SOLVING

[beam, mesh] = solveStaticBWGEB(beam, mesh);



%% PLOTS

figure('Color',[0 0 0],'Position',[10 50 1200 600]);

% snapshots
subplot(2,2,1)

for i = unique([0:80:beam.numLoadSteps,beam.numLoadSteps])    
    [nodalphi0, ~, ~, ~] = deassembleGlobalQ(beam.q(:,i+1), mesh);

    if ~(all(nodalphi0(2,:)==zeros(1,mesh.numNodes)))
        fprintf('The beam deformed out of the X1-X3 plane.\n')
    end    
    
    hold on
    plot(nodalphi0(3,:),nodalphi0(1,:),'LineWidth',2)
end

hold off
axis equal

xlabel('X1','Interpreter','latex'); 
ylabel('X3','Interpreter','latex');
grid on; box on;
title('Snapshots','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','FontSize',16)

clear i nodalphi0

% load function
subplot(2,2,2)
plot(beam.ttVec,beam.Fend(1) .* [0,beam.loadFactors],'LineWidth',2)

xlabel('$t$ [s]','Interpreter','latex'); 
ylabel('$F_{z,end}$ [N]','Interpreter','latex');
title('Load function $F(t)$ [N]','Interpreter','latex')
grid on; box on;
set(gca,'TickLabelInterpreter','latex','FontSize',16)


% stress-strain curve of 1 chosen component (of the last element and over
% load steps)
compId = 5;        % index of the chosen stress/strain component
sh = beam.s(compId:6:end,:); sh = sh(end,:);
eh = beam.e(compId:6:end,:); eh = eh(end,:);

subplot(2,2,3)
plot(eh,sh,'LineWidth',2)
xlabel('$e_{i,h}$ [-]','Interpreter','latex'); 
ylabel('$s_{i,h}$','Interpreter','latex');

grid on; box on;
title(sprintf('%.1d-th stress-strain component',compId),'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','FontSize',16)

clear compId eh sh


% tip displacement
subplot(2,2,4)
u3end = zeros(beam.numLoadSteps+1,1);

for i = 1:beam.numLoadSteps
    [nodalphi0, ~, ~, ~] = deassembleGlobalQ(beam.q(:,i+1), mesh);

    u3end(i+1) = nodalphi0(1,end);  % Store the tip displacement for each load step
end

plot(beam.ttVec,u3end,'LineWidth',2)
xlabel('$t$ [s]','Interpreter','latex'); 
ylabel('$u_{z,end}$','Interpreter','latex');

grid on; box on;
title('Tip displacement','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','FontSize',16)

clear i nodalphi0 u3End


