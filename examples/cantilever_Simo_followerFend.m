clear
close all
clc

% static analysis of GEB examples
% benchmark: cantilever straight beam subjected to a following end force
%            from [Example 7.2, Simo, 1986, CMAME]
% ne >= 32 = convergence to Simo's results (visual)


%% INPUT

beam.L = 100;
beam.EA = 42e7;
beam.GA1 = 1.61538e8;
beam.GA2 = 1.61538e8;
beam.EI1 = 3.5e7;
beam.EI2 = 3.5e7;
beam.Req = sqrt(20/pi);
beam.Jeq = pi * (beam.Req)^4/2;
beam.GIt = (2.1/2/1.3*beam.Jeq)*1e7;

% scaling factor for conditioning
beam.condScFac = 1e-3;

% external forces
beam.Fmax = 130e3 * beam.condScFac;

% number of load steps
beam.numLoadSteps = 1000;

% boundary conditions
beam.BCs = {'clamped', 'free'};     % boundary conditions at the left and right ends - 'clamped', 'fixed', 'free', 'custom'
beam.customFixedDofs = {[], []};    % if 'custom': array of local dof number of fixed dofs at the left and right ends

% discretization
mesh.nelms = 32;



%% PRE-PROCESSING

[mesh, beam] = preprocessing(mesh, beam);

% node vector of the reference configuration
beam.nodalTheta = linspace(0, beam.L, mesh.numNodes);
beam.nodeVec = [beam.nodalTheta; zeros(1,mesh.numNodes); zeros(1,mesh.numNodes)];

% nodal directors of the reference configuration (orthonormal)
beam.nodalD3 = [ones(1,mesh.nelms+1); zeros(1,mesh.nelms+1); zeros(1,mesh.nelms+1)];
beam.nodalD1 = [zeros(1,mesh.nelms+1); zeros(1,mesh.nelms+1); ones(1,mesh.nelms+1)];
beam.nodalD2 = [zeros(1,mesh.nelms+1); ones(1,mesh.nelms+1); zeros(1,mesh.nelms+1)];


beam.dsdtJacobianFunc = @(t) 1;     % function of the jacobian ds/dt
                                    % where t is the parametric coorcinates and s the arc length
                                    % (if straight beam: s = t and hence ds/dt = @(t) 1)


% zero generalized force vector (required for the routine)
beam.fext = zeros(mesh.ndofs,1);

% load factors
beam.loadFactors = linspace(1/beam.numLoadSteps,1,beam.numLoadSteps);



%% SOLVING

[beam, mesh] = solveStaticFollowingFendGEB(beam,mesh);



%% PLOTS

colors = [  [0.5, 0.5, 0.5];
            [0, 0.4470, 0.7410];
            [0.4660, 0.6740, 0.1880];
            [0.8500, 0.3250, 0.0980];
            [0.4940, 0.1840, 0.5560];
            [0.9290, 0.6940, 0.1250];
            [0.3010, 0.7450, 0.9330];
            [0.6350, 0.0780, 0.1840]];

% load-deflection curve of the free end
uEnd = zeros(beam.numLoadSteps+1,3);

for i = 0:beam.numLoadSteps
    [nodalphi0, ~, ~, ~] = deassembleGlobalQ(beam.q(:,i+1), mesh);
    uEnd(i+1,:) = nodalphi0(:,end) - beam.nodeVec(:,end);
end

figure('Color',[1 1 1],'Position',[10 50 900 600]);

plot((beam.Fmax/beam.condScFac/1e3/beam.numLoadSteps) .* (0:beam.numLoadSteps), uEnd(:,3),'b','LineWidth',2,'DisplayName','vertical $u_z$')
hold on;
plot((beam.Fmax/beam.condScFac/1e3/beam.numLoadSteps) .*(0:beam.numLoadSteps), -uEnd(:,1),'b-.','LineWidth',2,'DisplayName','horizontal $-u_x$')
hold on;
plot((beam.Fmax/beam.condScFac/1e3/beam.numLoadSteps) .*(0:beam.numLoadSteps), uEnd(:,2),'r','LineWidth',2,'DisplayName','$u_y$')

xlabel('$F_{end}$ [kN]','Interpreter','latex'); 
ylabel('$u_{End}$ [m]','Interpreter','latex');
legend('Interpreter','latex');
grid on; box on;
set(gca,'TickLabelInterpreter','latex','FontSize',16,'XLim',[0,130],'YLim',[-50,120])

clear uEnd i nodalphi0 



% snapshots
figure('Color',[1 1 1],'Position',[10 50 800 600]);

cc = 1;
for i = [0,39,154,385,539,693,846,beam.numLoadSteps]
    hold on
    [nodalphi0, ~, ~, ~] = deassembleGlobalQ(beam.q(:,i+1), mesh);

    if ~(all(nodalphi0(2,:)==zeros(1,mesh.numNodes)))
        fprintf('The beam deformed out of the X1-X3 plane.\n')
    end
    
    if cc >8
        plot(nodalphi0(1,:),nodalphi0(3,:),'LineWidth',2,'DisplayName',"GEB")    
        hold on
    else
        plot(nodalphi0(1,:),nodalphi0(3,:),'Color',colors(cc,:),'LineWidth',2,'DisplayName',"GEB")    
        hold on
    end
    cc = cc + 1;
end

axis equal

xlabel('X1','Interpreter','latex'); 
ylabel('X3','Interpreter','latex');
grid on; box on;
set(gca,'TickLabelInterpreter','latex','FontSize',16)

clear cc i nodalphi0 colors

