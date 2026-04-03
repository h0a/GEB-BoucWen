function h = computeHVector(q)
% compute the vector h of the orthonomal NODAL director constraints
% input:  q: nodal dofs
% output: h: nodal vector h of the size 6x1 (6 constraints per node)

d1 = q(4:6);
d2 = q(7:9);
d3 = q(10:12);

h = 0.5 .* [dot(d1,d1)-1;
            dot(d2,d2)-1;
            dot(d3,d3)-1;
            2*dot(d2,d3);
            2*dot(d1,d3);
            2*dot(d1,d2)];