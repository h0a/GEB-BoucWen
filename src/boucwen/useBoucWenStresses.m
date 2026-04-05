function beam = useBoucWenStresses(beam)
% add element stress function sigma_i(e,z) using the Bouc-Wen hysteretic model 
% and related matrices for solving steps


% check and set default parameters alpha and beta for the stress function
% sigma_i(e,z) if missing
if ~isfield(beam, 'alpha4Sigma')
    beam.alpha4Sigma = ones(6,1);
end

if ~isfield(beam, 'beta4Sigma')
    beam.beta4Sigma = ones(6,1);
end

% related matrices
beam.AmatSigma = spdiags(beam.alpha4Sigma, 0, 6, 6);
beam.BmatSigma = spdiags(beam.beta4Sigma,  0, 6, 6);

% stress function
beam.eleSigmaFunc = @(e,z) beam.AmatSigma * beam.elasMatrix * e + (beam.condScFac .* beam.BmatSigma) * (speye(6)-beam.AmatSigma) * z;

% default parameters regarding the solver for z when using the Bouc-Wen hysteretic model
beam.BWflag = 0;     % flag to stop the computation when the NR scheme does not converge
                     % 0: converged      1: non-converged

if ~isfield(beam,'maxNumBWIter')
    beam.maxNumBWIter = 50;       % maximum iteration step
end

if ~isfield(beam,'BWNRtol')
    beam.BWNRtol = 1e-10;         % tolerance for the newton-raphson scheme
end

% maximal number of iterations for the Uzawa scheme
if ~isfield(beam,'maxUzawaIter')
    beam.maxUzawaIter = 100;       % maximum iteration step
end

if ~isfield(beam,'UzawaTol')
    beam.UzawaTol = 1e-6;         % tolerance for the Uzawa scheme
end