function [T,P] = mesh_square_triangle(left,right,bottom,top,Nx,Ny)
% generate mesh data of structural triangular mesh on a rectangular area. 

map_element = @(r,c,idx) ((r-1)*Nx+c)*2-2+idx;
map_node = @(r,c) (r-1)*(Nx+1)+c;

P = zeros(2,(Nx+1)*(Ny+1));
for r = 1:Ny+1
    for c = 1:Nx+1
        P(:,map_node(r,c)) = [(c-1)*(right-left)/Nx+left;(r-1)*(top-bottom)/Ny+bottom];
    end
end
T = zeros(3,Nx*Ny*2);
for r = 1:Ny
    for c = 1:Nx
        T(:,map_element(r,c,1)) = [map_node(r,c);map_node(r,c+1);map_node(r+1,c)];
        T(:,map_element(r,c,2)) = [map_node(r+1,c+1);map_node(r+1,c);map_node(r,c+1)];
    end
end


