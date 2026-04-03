function fext = assembleFextRollingM(beam, mesh, alpha, x)

% allocation
fext = zeros(mesh.ndofs,1);

% rollup moment at the current load step
m = beam.MmaxEnd .* alpha;

% nodal variables at the right end
eleId = mesh.nelms;
dofsInd = mesh.LM(:,eleId);

% q of the last element
q = x(dofsInd);

% external force vector
% based on [Betsch, 2002, IJNME, Eq. 7]
% i-th block: fi = mi (delta dj * dk)
for i = 1:3
    if i == 1
        rowIds = dofsInd(19:21); % delta d2
        dk = q(end-2:end);       % d3
    elseif i == 2
        rowIds = dofsInd(22:24); % delta d3
        dk = q(16:18);           % d1
    else
        rowIds = dofsInd(16:18); % delta d1
        dk = q(19:21);           % d2
    end
    fext(rowIds) = m(i) .* dk;
end

% based on [Gebhardt, 2020, Journal of Nonlinear Science] - same results as
% based on [Betsch, 2002, IJNME]
% dofsInd = mesh.LM(:,eleId);
% ids1 = dofsInd(16:18);
% ids2 = dofsInd(end-2:end);
% d1 = q(16:18);
% d3 = q(22:24);
% fext(ids1) = 0.5*m .* d3;
% fext(ids2) = -0.5*m .* d1;

% based on force pair acting on the last 2 nodes - did not work
% dofsInd = mesh.LM(:,eleId);
% ids1 = dofsInd(1:3);
% ids2 = dofsInd(13:15);
% 
% dphi0 = q(13:15) - q(1:3);
% d1 = q(4:6);
% fext(ids1) = -m/norm(cross(d1,dphi0)) .* d1;
% 
% d1 = q(16:18);
% fext(ids2) = m/norm(cross(d1,dphi0)) .* d1;

