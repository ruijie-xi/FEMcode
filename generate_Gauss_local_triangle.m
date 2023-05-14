function [Gauss_weights,Gauss_pts]=generate_Gauss_local_triangle(vertices,Gauss_type)
%Generate the Gauss coefficients and Gauss points on the local triangle by using affine tranformation
%Some details are in my 2 dimensional bilinear IFE notes(part 8)
%Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle: the Gauss coefficients and Gauss points on the reference triangle
%vertices_triangle: the coordinates of the vertices of the local triangle 
%Gauss_coefficient_local_triangle,Gauss_point_local_triangle:the Gauss coefficients and Gauss points on the local triangle.

x1=vertices(1,1);
y1=vertices(2,1);
x2=vertices(1,2);
y2=vertices(2,2);
x3=vertices(1,3);
y3=vertices(2,3);
Jacobi=abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));

[Gauss_weights_ref,Gauss_pts_ref] = generate_Gauss_reference_triangle(Gauss_type);
Gauss_weights=Gauss_weights_ref*Jacobi;
Gauss_pts = Gauss_pts_ref;
Gauss_pts(1,:)=x1+(x2-x1)*Gauss_pts_ref(1,:)+(x3-x1)*Gauss_pts_ref(2,:);
Gauss_pts(2,:)=y1+(y2-y1)*Gauss_pts_ref(1,:)+(y3-y1)*Gauss_pts_ref(2,:);