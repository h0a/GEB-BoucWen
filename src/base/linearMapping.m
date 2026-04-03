function x = linearMapping(x0, x1, xi)
% Maps natural coordinate xi to the global coordinate x.
% x0: Coordinate of element left node.
% x1: Coordinate of element right node.
% xi: Vector of natural coordinates.
% x:  Vector of points in the global x domain. x = Q(xi).

x = (1-xi)/2 * x0 + (1+xi)/2 * x1;