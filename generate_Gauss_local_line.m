function [weights,pts] = generate_Gauss_local_line(vertices,Gauss_type)
% compute Gaussian pts and weights on edges.

[ref_weights,ref_pts] = generate_Gauss_reference_line(Gauss_type);

pt1 = vertices(:,1);pt2 = vertices(:,2);
mid = (pt1+pt2)/2;
arclength = sqrt(sum((pt2-pt1).^2));

pts = mid + ref_pts.*(pt2-pt1)/2;
weights = ref_weights*arclength/2;