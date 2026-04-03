function H = elementHmatrix(mesh, eleId, x)
% compute the element H matrix
% input:  mesh:      struct incl. relevant information of the discretization
%         eleId:     index of the current element
%         x:         vector of the solution of the current configuration
% output: H:         constraint residual vector 12x24 (12 constraint per element x 24 dofs per element)


% (nodal) q of the current element
qe = x(mesh.LM(:,eleId));

% q at left and right node
qleft = qe(1:mesh.numDofsPerNode);
qright = qe(mesh.numDofsPerNode+1:end);

% H-matrix at each node
Hleft = computeHmatrix(qleft);
Hright = computeHmatrix(qright);

% element H-matrix - factor 0.5 due to repeating at each node during the assembly routine - except the first and last node
if mesh.nelms == 1
    H = blkdiag(Hleft, Hright);
else
    if eleId == 1
        H = blkdiag(Hleft, 0.5 .* Hright);
    elseif eleId == mesh.nelms
        H = blkdiag(0.5 .* Hleft, Hright);
    else
        H = 0.5 .* blkdiag(Hleft, Hright);
    end
end
