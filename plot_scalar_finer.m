function plot_scalar_finer(vector,mesh,mesh_FE)
% plot a scalar function. create 3 polygons in each triangle element.

clf;
X = zeros(4,mesh.N_elem*3);
Y = zeros(4,mesh.N_elem*3);
Z = zeros(4,mesh.N_elem*3);

for i_elem = 1:mesh.N_elem
    vertices = mesh.P(:,mesh.T(:,i_elem));
    x1 = vertices(1,1);
    x2 = vertices(1,2);
    x3 = vertices(1,3);
    y1 = vertices(2,1);
    y2 = vertices(2,2);
    y3 = vertices(2,3);

    x0 = (x1+x2+x3)/3;
    y0 = (y1+y2+y3)/3;

    x12 = (x1+x2)/2;
    y12 = (y1+y2)/2;
    x31 = (x3+x1)/2;
    y31 = (y3+y1)/2;
    x23 = (x2+x3)/2;
    y23 = (y2+y3)/2;

    pts = [x1,x2,x3,x12,x23,x31,x0;y1,y2,y3,y12,y23,y31,y0];
    val = FE_function_local_2D(pts,vector,mesh,mesh_FE,i_elem,1,0,0);
    X(:,3*i_elem-2) = [x1,x12,x0,x31];
    X(:,3*i_elem-1) = [x2,x23,x0,x12];
    X(:,3*i_elem) = [x3,x31,x0,x23];

    Y(:,3*i_elem-2) = [y1,y12,y0,y31];
    Y(:,3*i_elem-1) = [y2,y23,y0,y12];
    Y(:,3*i_elem) = [y3,y31,y0,y23];
    
    Z(:,3*i_elem-2) = [val(1),val(4),val(7),val(6)];
    Z(:,3*i_elem-1) = [val(2),val(5),val(7),val(4)];
    Z(:,3*i_elem) = [val(3),val(6),val(7),val(5)];
    
end
patch(X,Y,Z,'LineStyle','none');
colorbar