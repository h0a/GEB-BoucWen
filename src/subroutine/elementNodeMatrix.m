function IEN = elementNodeMatrix(mesh)

% IEN(nc,e) = ng; 
% where: nc is the local node number of element number e
%        ng is the global node number

% standard numbering
IEN = [1:mesh.nelms; [1:mesh.nelms]+1];