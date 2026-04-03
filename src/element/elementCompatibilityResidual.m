function r = elementCompatibilityResidual(beam, mesh, eleId, x)
% compute the compatibility residual
% input:  beam:      struct incl. relevant information of the beam structure
%         mesh:      struct incl. relevant information of the discretization
%         eleId:     index of the current element
%         x:         vector of the solution of the current configuration
% output: r:         compatibility residual vector 6x1 (6 strain/stress components per element)


% (nodal) q and (element) e of the current element
qe = x(mesh.LM(:,eleId));
ee = x(mesh.LMstrains(:,eleId));

% (nodal) Q of the current element in the reference configuration
Qe = beam.q(mesh.LM(:,eleId),1);

% node id of the current elements
n1 = mesh.IEN(1,eleId);
n2 = mesh.IEN(2,eleId);

% parametric coordinates of the nodes of the current elements
t1 = beam.nodalTheta(n1);
t2 = beam.nodalTheta(n2);

% shape functions on [-1,1]
L = linearLagrangePolynomial(-1,1);
L1 = L{1}{1};
dL1 = L{2}{1};
L2 = L{1}{2};
dL2 = L{2}{2};

% gauss points and weights on [-1,1]
[xi, w] = GaussPoints(mesh.numGaussPtsPerEle);

% gauss points in parametric coordinates of the current element
t = linearMapping(t1, t2, xi);

% jacobian ds/dt (not necessarily constant)
J1 = beam.dsdtJacobianFunc(t);

% jacobian dt/dxi (constant - 1D - linear mapping)
J2 = (t2-t1)/2;

% allocation element compatibility residual vector
r = zeros(mesh.numStrainsPerEle,1);

% integration
for i = 1:mesh.numGaussPtsPerEle
    % interpolate q of the current and reference configuration
    qh =    L1(xi(i)) .* qe(1:12) +  L2(xi(i)) .* qe(13:end);
    dqh = (dL1(xi(i)) .* qe(1:12) + dL2(xi(i)) .* qe(13:end)) ./ (J1(i) * J2);
    
    Qh =    L1(xi(i)) .* Qe(1:12) +  L2(xi(i)) .* Qe(13:end);
    dQh = (dL1(xi(i)) .* Qe(1:12) + dL2(xi(i)) .* Qe(13:end))./ (J1(i) * J2);

    % displacement-based strains
    e = displacementBasedStrains(qh, dqh, Qh, dQh);    
    
    % compatibility residual
    r = r + (e - ee) .* (w(i) * J1(i) * J2);    
end
