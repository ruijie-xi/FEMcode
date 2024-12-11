function result = compute_dof(fun,mesh,mesh_FE,i)
% compute the degrees of freedom of a certain function.

if strcmp(mesh_FE.basis_type(1),'P') || strcmp(mesh_FE.basis_type,'bubbleP1') || strcmp(mesh_FE.basis_type,'DGP1') || strcmp(mesh_FE.basis_type,'DG-P2-quad')
    result = fun(mesh_FE.P(:,i));
    result = result(mesh_FE.idx(i));
elseif strcmp(mesh_FE.basis_type,'CR') || strcmp(mesh_FE.basis_type,'CR-P0') || strcmp(mesh_FE.basis_type,'CR-RT0') 
    result = fun(mesh_FE.P(:,i));
    result = result(mesh_FE.idx(i));
elseif strcmp(mesh_FE.basis_type,'BR0') || strcmp(mesh_FE.basis_type,'BR-RT0')
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
elseif strcmp(mesh_FE.basis_type,'BR')
    if i<=mesh_FE.N_dof_node*mesh.N_node
        i_node = ceil(i/mesh_FE.N_dof_node);
        j = i-(i_node-1)*mesh_FE.N_dof_node;
        val = fun(mesh.P(:,i_node));
        result = val(j);
    else
        i_edge = ceil((i-mesh_FE.N_dof_node*mesh.N_node)/mesh_FE.N_dof_edge);
        i_edge_in_elem = mesh.ET(2,i_edge);
        i_elem = mesh.E(2,i_edge);
        vertices = mesh.P(:,mesh.T(:,i_elem));
        if(i_edge_in_elem==1)
            endpts_h = [0,1;0,0];
            normal_hat = [0;-1];
        elseif(i_edge_in_elem==2)
            endpts_h = [1,0;0,1];
            normal_hat = [1/sqrt(2);1/sqrt(2)];
        elseif(i_edge_in_elem==3)
            endpts_h = [0,0;1,0];
            normal_hat = [-1;0];
        end

        if(mesh.TE(i_edge_in_elem,i_elem)<0)
            endpts_h(:,[2,1]) = endpts_h(:,[1,2]);
            normal_hat = -normal_hat;
        end

        [~,pts_line_h] = generate_Gauss_reference_line(5);
        pts_line_h = 0.5+0.5*pts_line_h;
        pts_ref = endpts_h(:,1) + (endpts_h(:,2)-endpts_h(:,1)).*pts_line_h;

        [J11,J12,J21,J22] = Jacobian_coeff(vertices,pts_ref);

        arc_ref = [1,sqrt(2),1];
        normal = zeros(2,size(pts_ref,2));
        normal(1,:) = (J22*normal_hat(1)-J21*normal_hat(2));
        normal(2,:) = (-J12*normal_hat(1)+J11*normal_hat(2));
        normal = normal/mesh.s(i_edge)*arc_ref(i_edge_in_elem);

        endpts = mesh.P(:,mesh.E([3,4],i_edge));
        [weights,pts] = generate_Gauss_local_line(endpts,5);
        val = fun(pts);
        result = (val(1,:).*normal(1,:)+val(2,:).*normal(2,:))*weights';

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
