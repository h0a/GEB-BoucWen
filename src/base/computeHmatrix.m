function H = computeHmatrix(q)
% compute the jacobian matrix associated to the vector h of the orthonomal NODAL director constraints
% input:  q: nodal dofs
% output: H: nodal matrix H of the size 6x12 (6 constraints per node x 12 dofs per node)

d1 = q(4:6);
d2 = q(7:9);
d3 = q(10:12);

H = sparse(6,12);

H(1,4:6)   = d1;
H(2,7:9)   = d2;
H(3,10:12) = d3;
H(4,7:12)  = [d3', d2'];
H(5,4:6)   = d3;
H(5,10:12) = d1;
H(6,4:9)   = [d2', d1'];
