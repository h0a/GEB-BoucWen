function fext = assembleFextFollowingF(beam, mesh, alpha, x)

% allocation
fext = zeros(mesh.ndofs,1);

% vector of the constant external moment (at the current load step)
f = beam.Fmax .* alpha;

% nodal variables at the right end
eleId = mesh.nelms;

% q of the last element
q = x(mesh.LM(:,eleId));

% external force vector
dofsInd = mesh.LM(:,eleId);
ids2 = dofsInd(13:15);
d1 = q(16:18);
fext(ids2) = f .* (d1);

