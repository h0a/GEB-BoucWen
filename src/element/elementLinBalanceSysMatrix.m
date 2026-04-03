function A11 = elementLinBalanceSysMatrix(beam, mesh, eleId, x)
% compute the element linearized system matrix associated to the balance residual incl. the contribution of the constraints
% input:  beam:      struct incl. relevant information of the beam structure
%         mesh:      struct incl. relevant information of the discretization
%         eleId:     index of the current element
%         x:         vector of the solution of the current configuration
% output: A11:       element linearized system matrix 24x24 (24 dofs of the primal fields per element)


% (element) s and langrange multipliers of the current element
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
A11 = zeros(2*mesh.numDofsPerNode);

% integration
for i = 1:mesh.numGaussPtsPerEle
    A11e = linearizeBtransposeS(se, [L1(xi(i)),L2(xi(i))], [dL1(xi(i)),dL2(xi(i))] ./ (J1(i) * J2));    
    
    A11 = A11 + A11e .* (w(i) * J1(i) * J2);    
end


%% contribution of the constraints (nodal value)

lagrangeMulLeft = lagrangeMule(1:mesh.numLagrangeMulPerNode);
lagrangeMulRight = lagrangeMule(mesh.numLagrangeMulPerNode+1:end);

A11hLeft = linearizeHtransposeV(lagrangeMulLeft);
A11hRight = linearizeHtransposeV(lagrangeMulRight);

% add a factor of 0.5 due to repeating at each node during the assembly routine - except the first and last node
if mesh.nelms == 1
    A11 = A11 + blkdiag(A11hLeft, A11hRight);
else
    if eleId == 1
        A11 = A11 + blkdiag(A11hLeft, 0.5 .* A11hRight);
    elseif eleId == mesh.nelms
        A11 = A11 + blkdiag(0.5 .* A11hLeft, A11hRight);
    else
        A11 = A11 + 0.5 .* blkdiag(A11hLeft, A11hRight);
    end
end
