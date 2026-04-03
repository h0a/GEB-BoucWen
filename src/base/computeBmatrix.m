function B = computeBmatrix(qe, shapeFuncs, shapeFuncs1stDeriv)
% compute the strain-displacement matrix B
% input:  qe:                 nodal q vector 24x1 of the current element
%         dqh:                interpolated vector of the 1st derivative of the dofs qh
%         shapeFuncs:         array of the 2 evaluated shape functions L1, L2 (functions associated to qh)
%         shapeFuncs1stDeriv: array of the evaluated 1st derivative of 2 shape functions L1, L2 (functions associated to qh)
% output: B: matrix B of the size 6x24 (6 strain components x 24 dofs per element) evaluated at 
%            the same point as shapeFuncs and shapeFuncs1stDeriv
   

% interpolate variable fields in q
qh = shapeFuncs(1) .* qe(1:12) + shapeFuncs(2) .* qe(13:end);
dqh = shapeFuncs1stDeriv(1) .* qe(1:12) + shapeFuncs1stDeriv(2) .* qe(13:end);

% extract relevant fields
d1h   = qh(4:6);
d2h   = qh(7:9);
d3h   = qh(10:12);

dphi0h = dqh(1:3);
dd1h   = dqh(4:6);
dd2h   = dqh(7:9);
dd3h   = dqh(10:12);

% compute B-matrix
Bi = zeros(6,12,2);

for i = 1:2         % loop over shape functions per element
    L = shapeFuncs(i);
    dL = shapeFuncs1stDeriv(i);

    Bi(1,1:3,i)   = d1h.*dL;
    Bi(1,4:6,i)   = dphi0h.*L;

    Bi(2,1:3,i)   = d2h.*dL;
    Bi(2,7:9,i)   = dphi0h.*L;

    Bi(3,1:3,i)   = d3h.*dL;
    Bi(3,10:12,i) = dphi0h.*L;

    Bi(4,7:9,i)   = 0.5 .* ( d3h.*dL - dd3h.*L );
    Bi(4,10:12,i) = 0.5 .* ( dd2h.*L - d2h.*dL );

    Bi(5,4:6,i)   = 0.5 .* ( dd3h.*L - d3h.*dL );
    Bi(5,10:12,i) = 0.5 .* ( d1h.*dL - dd1h.*L );

    Bi(6,4:6,i)   = 0.5 .* ( d2h.*dL - dd2h.*L );
    Bi(6,7:9,i)   = 0.5 .* ( dd1h.*L - d1h.*dL );
end

B = sparse([Bi(:,:,1)  Bi(:,:,2)]);