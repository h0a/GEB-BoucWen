function B = elementBmatrix(beam, mesh, eleId, x)
% compute the element B matrix
% input:  beam:      struct incl. relevant information of the beam structure
%         mesh:      struct incl. relevant information of the discretization
%         eleId:     index of the current element
%         x:         vector of the solution of the current configuration
% output: B:         element B matrix 6x24 (6 straint components x 24 dofs of the primal fields per element)


% (nodal) q of the current element
qe = x(mesh.LM(:,eleId));

% node id of the current elements
n1 = mesh.IEN(1,eleId);
n2 = mesh.IEN(2,eleId);

% parametric coordinates of the nodes of the current elements
t1 = beam.nodalTheta(n1);
t2 = beam.nodalTheta(n2);

% shape functions on [-1,1]
% L = linearLagrangePolynomial(-1,1);
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
B = zeros(mesh.numStrainsPerEle, 2*mesh.numDofsPerNode);

% integration
for i = 1:mesh.numGaussPtsPerEle
    Be = computeBmatrix(qe, [L1(xi(i)),L2(xi(i))], [dL1(xi(i)),dL2(xi(i))] ./ (J1(i) * J2));
    
    B = B + Be .* (w(i) * J1(i) * J2);    
end