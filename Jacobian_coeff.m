function [J11,J12,J21,J22] = Jacobian_coeff(vertices,pts_h)

x1=vertices(1,1);
y1=vertices(2,1);
x2=vertices(1,2);
y2=vertices(2,2);
x3=vertices(1,3);
y3=vertices(2,3);

J11 = (x2-x1)*ones(1,size(pts_h,2));
J12 = (x3-x1)*ones(1,size(pts_h,2));
J21 = (y2-y1)*ones(1,size(pts_h,2));
J22 = (y3-y1)*ones(1,size(pts_h,2));