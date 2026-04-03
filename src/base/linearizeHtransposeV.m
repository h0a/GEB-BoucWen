function V = linearizeHtransposeV(v)
% compute the linearization of (H^T(q) * v) that is a 12x12 matrix V (12 constraints per node)
% note that V then associates to the node where vector v is evaluated

V = sparse(12,12);

V(4:end,4:end) = [v(1).*eye(3) v(6).*eye(3) v(5).*eye(3);
                  v(6).*eye(3) v(2).*eye(3) v(4).*eye(3);
                  v(5).*eye(3) v(4).*eye(3) v(3).*eye(3)];