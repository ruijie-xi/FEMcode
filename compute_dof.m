function result = compute_dof(fun,mesh,mesh_FE,i)
% compute the degrees of freedom of a certain function.

if strcmp(mesh_FE.basis_type(1),'P') || strcmp(mesh_FE.basis_type,'bubbleP1') || strcmp(mesh_FE.basis_type,'DGP1') || strcmp(mesh_FE.basis_type,'DG-P2-quad')
    result = fun(mesh_FE.P(:,i));
    result = result(mesh_FE.idx(i));
elseif strcmp(mesh_FE.basis_type,'CR') || strcmp(mesh_FE.basis_type,'CR-P0') || strcmp(mesh_FE.basis_type,'CR-RT0') 
    result = fun(mesh_FE.P(:,i));
    result = result(mesh_FE.idx(i));
elseif strcmp(mesh_FE.basis_type,'BR') || strcmp(mesh_FE.basis_type,'BR-RT0')
    if i<mesh.N_node+1
        result = fun(mesh_FE.P(:,i));
        result = result(1);
    elseif i<mesh.N_node*2+1
        result = fun(mesh_FE.P(:,i));
        result = result(2);
    else
        i_edge = mesh_FE.P(1,i);
        endpts = [mesh.P(:,mesh.E(3,i_edge)),mesh.P(:,mesh.E(4,i_edge))];
        [weights,pts] = generate_Gauss_local_line(endpts,5);
        val = fun(pts);
        result = (val(1,:)*mesh.N(1,i_edge)+val(2,:)*mesh.N(2,i_edge))*weights';
        result = result/mesh.s(i_edge);
    end
elseif strcmp(mesh_FE.basis_type,'RT0')
    i_edge = mesh_FE.P(1,i);
    endpts = [mesh.P(:,mesh.E(3,i_edge)),mesh.P(:,mesh.E(4,i_edge))];
    [weights,pts] = generate_Gauss_local_line(endpts,5);
    val = fun(pts);
    result = (val(1,:)*mesh.N(1,i_edge)+val(2,:)*mesh.N(2,i_edge))*weights';
    result = result/mesh.s(i_edge);
elseif strcmp(mesh_FE.basis_type,'MTW')
    if i<=mesh.N_edge
        i_edge = mesh_FE.P(i);
        endpts = [mesh.P(:,mesh.E(3,i_edge)),mesh.P(:,mesh.E(4,i_edge))];
        [weights,pts] = generate_Gauss_local_line(endpts,5);
        val = fun(pts);
        s = sqrt(sum((pts-endpts(:,1)).^2))./sqrt(sum((endpts(:,2)-endpts(:,1)).^2));
        result = (val(1,:)*mesh.N(1,i_edge)+val(2,:)*mesh.N(2,i_edge)).*s*weights';
        result = result/mesh.s(i_edge);
    elseif i<=2*mesh.N_edge
        i_edge = mesh_FE.P(i);
        endpts = [mesh.P(:,mesh.E(3,i_edge)),mesh.P(:,mesh.E(4,i_edge))];
        [weights,pts] = generate_Gauss_local_line(endpts,5);
        val = fun(pts);
        s = sqrt(sum((pts-endpts(:,1)).^2))./sqrt(sum((endpts(:,2)-endpts(:,1)).^2));
        result = (val(1,:)*mesh.N(1,i_edge)+val(2,:)*mesh.N(2,i_edge)).*(1-s)*weights';
        result = result/mesh.s(i_edge);
    elseif i<=3*mesh.N_edge
        i_edge = mesh_FE.P(i);
        endpts = [mesh.P(:,mesh.E(3,i_edge)),mesh.P(:,mesh.E(4,i_edge))];
        [weights,pts] = generate_Gauss_local_line(endpts,5);
        val = fun(pts);
        result = (val(1,:)*mesh.tau(1,i_edge)+val(2,:)*mesh.tau(2,i_edge))*weights';
        result = result/mesh.s(i_edge);
    end
elseif strcmp(mesh_FE.basis_type,'RT0-MTW')
    if i<=mesh.N_edge
        i_edge = mesh_FE.P(i);
        endpts = [mesh.P(:,mesh.E(3,i_edge)),mesh.P(:,mesh.E(4,i_edge))];
        [weights,pts] = generate_Gauss_local_line(endpts,5);
        val = fun(pts);
        result = (val(1,:)*mesh.N(1,i_edge)+val(2,:)*mesh.N(2,i_edge))*weights';
        result = result/mesh.s(i_edge);
    elseif i<=2*mesh.N_edge
        i_edge = mesh_FE.P(i);
        endpts = [mesh.P(:,mesh.E(3,i_edge)),mesh.P(:,mesh.E(4,i_edge))];
        [weights,pts] = generate_Gauss_local_line(endpts,5);
        val = fun(pts);
        result = (val(1,:)*mesh.tau(1,i_edge)+val(2,:)*mesh.tau(2,i_edge))*weights';
        result = result/mesh.s(i_edge);
    end
elseif strcmp(mesh_FE.basis_type,'Ned1-1')
    i_edge = mesh_FE.P(i);
    endpts = [mesh.P(:,mesh.E(3,i_edge)),mesh.P(:,mesh.E(4,i_edge))];
    [weights,pts] = generate_Gauss_local_line(endpts,5);
    val = fun(pts);
    result = (val(1,:)*mesh.tau(1,i_edge)+val(2,:)*mesh.tau(2,i_edge))*weights';
    result = result/mesh.s(i_edge);

elseif strcmp(mesh_FE.basis_type,'Ned1-2')
    if i<=mesh.N_edge
        i_edge = mesh_FE.P(i);
        endpts = [mesh.P(:,mesh.E(3,i_edge)),mesh.P(:,mesh.E(4,i_edge))];
        [weights,pts] = generate_Gauss_local_line(endpts,5);
        val = fun(pts);
        s = sqrt(sum((pts-endpts(:,1)).^2))./sqrt(sum((endpts(:,2)-endpts(:,1)).^2));
        result = (val(1,:)*mesh.tau(1,i_edge)+val(2,:)*mesh.tau(2,i_edge)).*s*weights';
        result = result/mesh.s(i_edge);
    elseif i<=mesh.N_edge*2
        i_edge = mesh_FE.P(i);
        endpts = [mesh.P(:,mesh.E(3,i_edge)),mesh.P(:,mesh.E(4,i_edge))];
        [weights,pts] = generate_Gauss_local_line(endpts,5);
        val = fun(pts);
        s = sqrt(sum((pts-endpts(:,1)).^2))./sqrt(sum((endpts(:,2)-endpts(:,1)).^2));
        result = (val(1,:)*mesh.tau(1,i_edge)+val(2,:)*mesh.tau(2,i_edge)).*(1-s)*weights';
        result = result/mesh.s(i_edge);
    elseif i<=2*mesh.N_edge+mesh.N_elem
        i_elem = mesh_FE.P(i);
        vertices = mesh.P(:,mesh.T(:,i_elem));
        y1 = vertices(2,1);y2 = vertices(2,2);y3 = vertices(2,3);
        x1 = vertices(1,1);x2 = vertices(1,2);x3 = vertices(1,3);
        J = [x2-x1,x3-x1;y2-y1,y3-y1];
        J0 = abs(det(J));
        p = (J*[1;0]/J0^(3/2))';
        [weights,pts] = generate_Gauss_local_triangle(vertices,9);
        val = fun(pts);
        result = p*val*weights';
    elseif i<=2*mesh.N_edge+2*mesh.N_elem
        i_elem = mesh_FE.P(i);
        vertices = mesh.P(:,mesh.T(:,i_elem));
        y1 = vertices(2,1);y2 = vertices(2,2);y3 = vertices(2,3);
        x1 = vertices(1,1);x2 = vertices(1,2);x3 = vertices(1,3);
        J = [x2-x1,x3-x1;y2-y1,y3-y1];
        J0 = abs(det(J));
        p = (J*[0;1]/J0^(3/2))';
        [weights,pts] = generate_Gauss_local_triangle(vertices,9);
        val = fun(pts);
        result = p*val*weights';
    end

elseif strcmp(mesh_FE.basis_type,'Ned2-1')
    if i<=mesh.N_edge
        i_edge = mesh_FE.P(i);
        endpts = [mesh.P(:,mesh.E(3,i_edge)),mesh.P(:,mesh.E(4,i_edge))];
        [weights,pts] = generate_Gauss_local_line(endpts,5);
        val = fun(pts);
        s = sqrt(sum((pts-endpts(:,1)).^2))./sqrt(sum((endpts(:,2)-endpts(:,1)).^2));
        result = (val(1,:)*mesh.tau(1,i_edge)+val(2,:)*mesh.tau(2,i_edge)).*s*weights';
        result = result/mesh.s(i_edge);
    elseif i<=mesh.N_edge*2
        i_edge = mesh_FE.P(i);
        endpts = [mesh.P(:,mesh.E(3,i_edge)),mesh.P(:,mesh.E(4,i_edge))];
        [weights,pts] = generate_Gauss_local_line(endpts,5);
        val = fun(pts);
        s = sqrt(sum((pts-endpts(:,1)).^2))./sqrt(sum((endpts(:,2)-endpts(:,1)).^2));
        result = (val(1,:)*mesh.tau(1,i_edge)+val(2,:)*mesh.tau(2,i_edge)).*(1-s)*weights';
        result = result/mesh.s(i_edge);
    end

else
    fprintf('cannot handle this type of basis!');
end
