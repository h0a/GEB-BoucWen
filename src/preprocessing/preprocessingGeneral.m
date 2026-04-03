function [mesh, beam] = preprocessingGeneral(mesh, beam)

%% default parameters

% define default scaling factor for conditioning
if ~isfield(beam,'condScFac')
    beam.condScFac = 1.0;
end

% default option to compute linearization numerically
if ~isfield(beam,'computeAnumerically')
    beam.computeAnumerically = false;
end

% default number of Gauss points per element
if ~isfield(mesh, 'numGaussPtsPerEle')
    mesh.numGaussPtsPerEle = 1;
end

% default parameters for the Newton-Raphson scheme
if ~isfield(beam,'maxNumIter')
    beam.maxNumIter = 50;       % maximum iteration step
end

if ~isfield(beam,'NRtol')
    beam.NRtol = 1e-10;         % tolerance for the newton-raphson scheme
end



%% elasticity matrix
beam.elasMatrix = elasticityMatrix(beam) .* beam.condScFac;


%% NUMBER OF DOFS AND LOCATION MATRIX OF THE PRIMAL AND SECONDARY FIELDS

% number of nodes
mesh.numNodes = mesh.nelms+1;

% ndofs
mesh.numDofsPerNode = 12;
mesh.numStrainsPerEle = 6;
mesh.numStressesPerEle = 6;

% number of (both active and non-active) dofs
mesh.ndofs = mesh.numDofsPerNode * mesh.numNodes;
mesh.ndofse = mesh.numStrainsPerEle * mesh.nelms;
mesh.ndofss = mesh.numStressesPerEle * mesh.nelms;

% location matrix
mesh = locationMatrix(mesh);
mesh = locationMatrixStrains(mesh);
mesh = locationMatrixStresses(mesh);


%% NUMBER OF LAGRANGE MULTIPLIERS AND LOCATION MATRIX
mesh.numLagrangeMulPerNode = 6;
mesh.numLagrangeMul = mesh.numLagrangeMulPerNode * mesh.numNodes;

% location matrix
mesh = locationMatrixLagrangeMul(mesh);


%% TOTAL DOFS
mesh.ndofsTotal = mesh.ndofs + mesh.ndofse + mesh.ndofss + mesh.numLagrangeMul;

