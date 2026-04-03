function h = elementArcLength(beam, mesh, eleId)
% compute the arc length of an element
% input:  beam:      struct incl. relevant information of the beam structure
%         mesh:      struct incl. relevant information of the discretization
%         eleId:     index of the current element
% output: h:       arc length of the current element


% node id of the current elements
n1 = mesh.IEN(1,eleId);
n2 = mesh.IEN(2,eleId);

% parametric coordinates of the nodes of the current elements
t1 = beam.nodalTheta(n1);
t2 = beam.nodalTheta(n2);

% gauss points and weights on [-1,1]
[xi, w] = GaussPoints(mesh.numGaussPtsPerEle);

% gauss points in parametric coordinates of the current element
t = linearMapping(t1, t2, xi);

% jacobian ds/dt (not necessarily constant)
J1 = beam.dsdtJacobianFunc(t);

% jacobian dt/dxi (constant - 1D - linear mapping)
J2 = (t2-t1)/2;

% allocation element balance residual vector
h = 0;

% integration
for i = 1:mesh.numGaussPtsPerEle
    h = h + (w(i) * J1(i) * J2);    
end



