clear
close all
clc

% static analysis of GEB examples
% benchmark: planar rollup of a cantilever straight beam
%            subjected to a moment at free end Mmax = 2*EI * pi / L;
%            the beam is rolled up in X1-X3 plane


%% INPUT

beam.L = 40;
beam.EA = 100;
beam.GA1 = 5e10;
beam.GA2 = 5e10;
beam.EI1 = 200;
beam.EI2 = 200;
beam.GIt = 2e10;

% scaling factor for conditioning
beam.condScFac = 1e-2;

% external forces
beam.bodyFext = @(t) [0; 0; 0] .* beam.condScFac;                       % config.-independent body force
beam.MmaxEnd = [0, 2*beam.EI2 * pi / beam.L * beam.condScFac, 0];       % rollup moments: M = [m1, m2, m3] = m1*d1 + m2*d2 + m3*d3
%                  m2*(-d2)

% number of load steps
beam.numLoadSteps = 55;

% boundary conditions
beam.BCs = {'clamped', 'free'};     % boundary conditions at the left and right ends - 'clamped', 'fixed', 'free', 'custom'
beam.customFixedDofs = {[], []};    % if 'custom': array of local dof number of fixed dofs at the left and right ends

% discretization
mesh.nelms = 16;


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

% nodal constant external force vectors (point force at each node)
beam.nodalFext = zeros(3,mesh.numNodes);

% load factors
beam.loadFactors = linspace(1/beam.numLoadSteps,1,beam.numLoadSteps);



%% SOLVING

[beam, mesh] = solveStaticRollUpGEB(beam,mesh);


%% POSTPROCESSING

colors = [  [0.5, 0.5, 0.5];
            [0, 0.4470, 0.7410];
            [0.4660, 0.6740, 0.1880];
            [0.8500, 0.3250, 0.0980];
            [0.4940, 0.1840, 0.5560];
            [0.9290, 0.6940, 0.1250];
            [0.3010, 0.7450, 0.9330];
            [0.6350, 0.0780, 0.1840]];

figure('Color',[1 1 1],'Position',[10 50 800 600]);

cc = 1;
for i = 0:11:beam.numLoadSteps
    hold on
    [nodalphi0, ~, ~, ~] = deassembleGlobalQ(beam.q(:,i+1), mesh);

    if ~(all(nodalphi0(2,:)==zeros(1,mesh.numNodes)))
        fprintf('The beam deformed out of the X1-X3 plane.\n')
    end
    
    if cc >8
        plot(nodalphi0(1,:),nodalphi0(3,:),'LineWidth',2,'DisplayName',"GEB")    
        hold on
        scatter(nodalphi0(1,:),nodalphi0(3,:),15,'o')
    else
        plot(nodalphi0(1,:),nodalphi0(3,:),'Color',colors(cc,:),'LineWidth',2,'DisplayName',"GEB")    
        hold on
        scatter(nodalphi0(1,:),nodalphi0(3,:),15,'o','Color',colors(cc,:))   
    end
    cc = cc + 1;
end

xlim([-10,45])
ylim([-10,45])
axis equal

xlabel('X1','Interpreter','latex'); 
ylabel('X3','Interpreter','latex');
grid on; box on;
set(gca,'TickLabelInterpreter','latex','FontSize',16)

clear cc i nodalphi0 colors