function e = displacementBasedStrains(qh, dqh, Qh, dQh)
% compute the displacement-based strains
% input:  qh:       interpolated vector of dofs qh
%         dqh:      interpolated vector of the 1st derivative of the dofs qh
%         Qh:       interpolated vector of dofs Qh in the reference configuration
%         dQh:      interpolated vector of the 1st derivative of the dofs Qh in the reference configuration
% output: e:        6x1 vector of strain components: shear strains in 2 directions, elongation, bending in 2 directions, and torsion

D1h   = Qh(4:6);
D2h   = Qh(7:9);
D3h   = Qh(10:12);

dPhi0h = dQh(1:3);
dD1h   = dQh(4:6);
dD2h   = dQh(7:9);
dD3h   = dQh(10:12);

d1h   = qh(4:6);
d2h   = qh(7:9);
d3h   = qh(10:12);

dphi0h = dqh(1:3);
dd1h   = dqh(4:6);
dd2h   = dqh(7:9);
dd3h   = dqh(10:12);

e = [dot(d1h,dphi0h)-dot(D1h,dPhi0h);
     dot(d2h,dphi0h)-dot(D2h,dPhi0h);
     dot(d3h,dphi0h)-dot(D3h,dPhi0h);
     0.5*( dot(d3h,dd2h) - dot(d2h,dd3h) - dot(D3h,dD2h) + dot(D2h,dD3h) );
     0.5*( dot(d1h,dd3h) - dot(d3h,dd1h) - dot(D1h,dD3h) + dot(D3h,dD1h) );
     0.5*( dot(d2h,dd1h) - dot(d1h,dd2h) - dot(D2h,dD1h) + dot(D1h,dD2h) )];