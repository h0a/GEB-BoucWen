function r = elementBalanceResidual(beam, mesh, eleId, x)
% compute the balance residual incl. the contribution of the constraints
% input:  beam:      struct incl. relevant information of the beam structure
%         mesh:      struct incl. relevant information of the discretization
%         eleId:     index of the current element
%         x:         vector of the solution of the current configuration
% output: r:         balance residual vector 24x1 (24 dofs of the primal fields per element)


% (nodal) q, (element) s, and langrange multipliers of the current element
qe = x(mesh.LM(:,eleId));
se = x(mesh.LMstresses(:,eleId));
lagrangeMule = x(mesh.LMlamu(:,eleId));

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

% allocation element balance residual vector
r = zeros(2*mesh.numDofsPerNode,1);

% integration
for i = 1:mesh.numGaussPtsPerEle
    B = computeBmatrix(qe, [L1(xi(i)),L2(xi(i))], [dL1(xi(i)),dL2(xi(i))] ./ (J1(i) * J2));    
    
    r = r + (B' * se) .* (w(i) * J1(i) * J2);    
end


%% contribution of the constraints (nodal value)

qleft = qe(1:mesh.numDofsPerNode);
qright = qe(mesh.numDofsPerNode+1:end);

Hleft = computeHmatrix(qleft);
Hright = computeHmatrix(qright);

lagrangeMulLeft = lagrangeMule(1:mesh.numLagrangeMulPerNode);
lagrangeMulRight = lagrangeMule(mesh.numLagrangeMulPerNode+1:end);

% add a factor of 0.5 due to repeating at each node during the assembly routine - except the first and last node
if mesh.nelms == 1
    r = r + [Hleft' * lagrangeMulLeft; Hright' * lagrangeMulRight];
else
    if eleId == 1
        r = r + [Hleft' * lagrangeMulLeft; 0.5 .* Hright' * lagrangeMulRight];
    elseif eleId == mesh.nelms
        r = r + [0.5 .* Hleft' * lagrangeMulLeft; Hright' * lagrangeMulRight];
    else
        r = r + 0.5 .* [Hleft' * lagrangeMulLeft; Hright' * lagrangeMulRight];
    end
end    
