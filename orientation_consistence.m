function [T,P] = orientation_consistence(T,P)
% ensure that all the triangles have a positive orientation.

for i_elem = 1:size(T,2)
    vertices = P(:,T(:,i_elem));
    J = [vertices(:,2)-vertices(:,1),vertices(:,3)-vertices(:,1)];
    if det(J)<0
        tmp = T(1,i_elem);
        T(1,i_elem) = T(2,i_elem);
        T(2,i_elem) = tmp;
    end
end