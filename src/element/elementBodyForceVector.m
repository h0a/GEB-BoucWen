function fe = elementBodyForceVector(beam, mesh, eleId)

    % node id of the current elements
    n1 = mesh.IEN(1,eleId);
    n2 = mesh.IEN(2,eleId);

    % parametric coordinates of the nodes of the current elements
    t1 = beam.nodalTheta(n1);
    t2 = beam.nodalTheta(n2);

    % shape functions on [-1,1]
    L = linearLagrangePolynomial(-1,1);
    L1 = L{1}{1};
    L2 = L{1}{2};

    % gauss points and weights on [-1,1]
    [xi, w] = GaussPoints(mesh.numGaussPtsPerEle);

    % gauss points in parametric coordinates of the current element
    t = linearMapping(t1, t2, xi);

    % jacobian ds/dt (not necessarily constant)
    J1 = beam.dsdtJacobianFunc(t);

    % jacobian dt/dxi (constant - 1D - linear mapping)
    J2 = (t2-t1)/2;

    % allocation element force vector
    fe = zeros(2*mesh.numDofsPerNode,1);

    % integration
    for i = 1:mesh.numGaussPtsPerEle
        fe(1:3) = fe(1:3) + ( L1(xi(i)) .* beam.bodyFext(t(i)) ) .* (w(i) * J1(i) * J2);
        fe(mesh.numDofsPerNode+1:mesh.numDofsPerNode+3) = fe(mesh.numDofsPerNode+1:mesh.numDofsPerNode+3) + ( L2(xi(i)) .* beam.bodyFext(t(i)) ) .* (w(i) * J1(i) * J2);
    end
end