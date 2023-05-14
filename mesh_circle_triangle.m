function [T,P] = mesh_circle_triangle(x0,y0,r,h0)

fd = @(p)dcircle(p,x0,y0,r);
[P,T] = distmesh2d(fd,@huniform,h0,[x0-r,y0-r;x0+r,y0+r],[]);
P = P';T = T';