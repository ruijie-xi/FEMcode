function [T,P] = mesh_square_triangle_unstructured(left,right,bottom,top,h0)

fd = @(p)dpoly(p,[left,bottom;right,bottom;right,top;left,top;left,bottom]);
[P,T] = distmesh2d(fd,@huniform,h0,[left,bottom;right,top],[]);
P = P';
T = T';