function result = basis_local_2D(coords,mesh,mesh_FE,i_elem,dx,dy)
% evaluate local basis functions (or its derivatives) on a certain element.

vertices = mesh.P(:,mesh.T(:,i_elem));

x = coords(1,:);
y = coords(2,:);
y1 = vertices(2,1);y2 = vertices(2,2);y3 = vertices(2,3);
x1 = vertices(1,1);x2 = vertices(1,2);x3 = vertices(1,3);
basis_type = mesh_FE.basis_type;
num_pts = size(coords,2);
num_basis = mesh_FE.N_lb;
dim = mesh_FE.dim;
zero = zeros(num_basis,num_pts);

if strcmp(basis_type(1),'P') || strcmp(basis_type,'CR') || strcmp(basis_type,'CR-P0') || strcmp(basis_type,'bubbleP1') || strcmp(basis_type,'DGP1') || strcmp(basis_type,'DG-P2-quad')
    % obtain reference coordinates [xh,yh]
    % Note that [x;y] = [x1;y1]+J*[xh;yh].
    J = [x2-x1,x3-x1;y2-y1,y3-y1];
    coords_h = J\([x;y]-vertices(:,1));
    J0 = det(J);

    % calculate local basis functions. (through basis_reference_2D)
    if dx==0 && dy==0
        result = basis_reference_2D(coords_h,basis_type,dx,dy);
    elseif dx==1 && dy==0
        dxh = basis_reference_2D(coords_h,basis_type,1,0);
        dyh = basis_reference_2D(coords_h,basis_type,0,1);
        result = dxh*(y3-y1)/J0+dyh*(y1-y2)/J0;
    elseif dx==0 && dy==1
        dxh = basis_reference_2D(coords_h,basis_type,1,0);
        dyh = basis_reference_2D(coords_h,basis_type,0,1);
        result = dxh*(x1-x3)/J0+dyh*(x2-x1)/J0;
    elseif dx==2 && dy==0
        dxxh = basis_reference_2D(coords_h,basis_type,2,0);
        dxyh = basis_reference_2D(coords_h,basis_type,1,1);
        dyyh = basis_reference_2D(coords_h,basis_type,0,2);
        result = dxxh*(y3-y1)^2/J0^2 + 2*dxyh*(y3-y1)*(y1-y2)/J0^2 + dyyh*(y1-y2)^2/J0^2;
    elseif dx==0 && dy==2
        dxxh = basis_reference_2D(coords_h,basis_type,2,0);
        dxyh = basis_reference_2D(coords_h,basis_type,1,1);
        dyyh = basis_reference_2D(coords_h,basis_type,0,2);
        result = dxxh*(x1-x3)^2/J0^2 + 2*dxyh*(x1-x3)*(x2-x1)/J0^2 + dyyh*(x2-x1)^2/J0^2;
    elseif dx==1 && dy==1
        dxxh = basis_reference_2D(coords_h,basis_type,2,0);
        dxyh = basis_reference_2D(coords_h,basis_type,1,1);
        dyyh = basis_reference_2D(coords_h,basis_type,0,2);
        result = dxxh*(x1-x3)*(y3-y1)/J0^2 + dxyh*(x1-x3)*(y1-y2)/J0^2 ...
            + dyyh*(x2-x1)*(y1-y2)/J0^2 + dxyh*(x2-x1)*(y3-y1)/J0^2;
    else
        result = zero;
    end
    result = repmat(result,dim,1,1);

elseif strcmp(basis_type,'CR-RT0')
    J = [x2-x1,x3-x1;y2-y1,y3-y1];
    coords_h = J\([x;y]-vertices(:,1));
    J0 = det(J);

    s = sign(mesh.TE(:,i_elem));
    s = [s;s];
    arc = mesh.s(abs(mesh.TE(:,i_elem)));
    arc = [arc,arc];
    
    n = mesh.N(:,abs(mesh.TE(:,i_elem)));
    N = [n(1,1),n(1,2),n(1,3),n(2,1),n(2,2),n(2,3)];
    A = diag(s.*arc'.*N');
    if dx==0 && dy==0
        basis_ref = basis_reference_2D(coords_h,basis_type,dx,dy);

        basis_ref = reshape(basis_ref,2,num_basis*num_pts);
        basis = (1/J0)*J*basis_ref;
        basis = reshape(basis,2,num_basis,num_pts);

        result(1,:,:) = A*reshape(basis(1,:,:),num_basis,num_pts);
        result(2,:,:) = A*reshape(basis(2,:,:),num_basis,num_pts);

    elseif dx==1 && dy==0
        basis_ref_dxh = basis_reference_2D(coords_h,basis_type,1,0);
        basis_ref_dyh = basis_reference_2D(coords_h,basis_type,0,1);
        basis_ref_dx = basis_ref_dxh*(y3-y1)/J0+basis_ref_dyh*(y1-y2)/J0;

        basis_ref_dx = reshape(basis_ref_dx,2,num_basis*num_pts);
        basis_dx = (1/J0)*J*basis_ref_dx;
        basis_dx = reshape(basis_dx,2,num_basis,num_pts);

        result(1,:,:) = A*reshape(basis_dx(1,:,:),num_basis,num_pts);
        result(2,:,:) = A*reshape(basis_dx(2,:,:),num_basis,num_pts);

    elseif dx==0 && dy==1
        basis_ref_dxh = basis_reference_2D(coords_h,basis_type,1,0);
        basis_ref_dyh = basis_reference_2D(coords_h,basis_type,0,1);
        basis_ref_dy =  basis_ref_dxh*(x1-x3)/J0+basis_ref_dyh*(x2-x1)/J0;

        basis_ref_dy = reshape(basis_ref_dy,2,num_basis*num_pts);
        basis_dy = (1/J0)*J*basis_ref_dy;
        basis_dy = reshape(basis_dy,2,num_basis,num_pts);

        result(1,:,:) = A*reshape(basis_dy(1,:,:),num_basis,num_pts);
        result(2,:,:) = A*reshape(basis_dy(2,:,:),num_basis,num_pts);
    else
        fprintf('cannot handle this type of derivative!');
    end

elseif strcmp(basis_type,'BR0') || strcmp(basis_type,'BR-RT0')
    % obtain reference coordinates [xh,yh]
    % Note that [x;y] = [x1;y1]+J*[xh;yh].
    J = [x2-x1,x3-x1;y2-y1,y3-y1];
    coords_h = J\([x;y]-vertices(:,1));
    J0 = det(J);

    n1 = [y2-y1;x1-x2];
    n2 = [y3-y2;x2-x3];
    n3 = [y1-y3;x3-x1];

    s = sign(mesh.TE(:,i_elem));

    A = [n1(1),0,0,n3(1),0,0,0,0,0;
        0,n1(1),0,0,n2(1),0,0,0,0;
        0,0,n3(1),0,0,n2(1),0,0,0;
        n1(2),0,0,n3(2),0,0,0,0,0;
        0,n1(2),0,0,n2(2),0,0,0,0;
        0,0,n3(2),0,0,n2(2),0,0,0;
        0,0,0,0,0,0,s(1)*norm(n1),0,0;
        0,0,0,0,0,0,0,s(2)*norm(n2),0;
        0,0,0,0,0,0,0,0,s(3)*norm(n3);];

    if dx==0 && dy==0
        basis_ref = basis_reference_2D(coords_h,basis_type,dx,dy);

        basis_ref = reshape(basis_ref,2,num_basis*num_pts);
        basis = (1/J0)*J*basis_ref;
        basis = reshape(basis,2,num_basis,num_pts);

        result(1,:,:) = A*reshape(basis(1,:,:),num_basis,num_pts);
        result(2,:,:) = A*reshape(basis(2,:,:),num_basis,num_pts);

    elseif dx==1 && dy==0
        basis_ref_dxh = basis_reference_2D(coords_h,basis_type,1,0);
        basis_ref_dyh = basis_reference_2D(coords_h,basis_type,0,1);
        basis_ref_dx = basis_ref_dxh*(y3-y1)/J0+basis_ref_dyh*(y1-y2)/J0;

        basis_ref_dx = reshape(basis_ref_dx,2,num_basis*num_pts);
        basis_dx = (1/J0)*J*basis_ref_dx;
        basis_dx = reshape(basis_dx,2,num_basis,num_pts);

        result(1,:,:) = A*reshape(basis_dx(1,:,:),num_basis,num_pts);
        result(2,:,:) = A*reshape(basis_dx(2,:,:),num_basis,num_pts);

    elseif dx==0 && dy==1
        basis_ref_dxh = basis_reference_2D(coords_h,basis_type,1,0);
        basis_ref_dyh = basis_reference_2D(coords_h,basis_type,0,1);
        basis_ref_dy =  basis_ref_dxh*(x1-x3)/J0+basis_ref_dyh*(x2-x1)/J0;

        basis_ref_dy = reshape(basis_ref_dy,2,num_basis*num_pts);
        basis_dy = (1/J0)*J*basis_ref_dy;
        basis_dy = reshape(basis_dy,2,num_basis,num_pts);

        result(1,:,:) = A*reshape(basis_dy(1,:,:),num_basis,num_pts);
        result(2,:,:) = A*reshape(basis_dy(2,:,:),num_basis,num_pts);
    else
        fprintf('cannot handle this type of derivative!');
    end

elseif strcmp(basis_type,'BR')
    % obtain reference coordinates [xh,yh]
    % Note that [x;y] = [x1;y1]+J*[xh;yh].
    J = [x2-x1,x3-x1;y2-y1,y3-y1];
    coords_h = J\([x;y]-vertices(:,1));
    J0 = det(J);

    s = sign(mesh.TE(:,i_elem));

    n1_hat = [1;0];
    n2_hat = [0;1];
    n1 = (J')\n1_hat*J0;
    n2 = (J')\n2_hat*J0;


    if dx==0 && dy==0
        basis_ref = basis_reference_2D(coords_h,basis_type,dx,dy);

        basis_ref = reshape(basis_ref,2,num_basis*num_pts);
        basis = (1/J0)*J*basis_ref;
        basis = reshape(basis,2,num_basis,num_pts);

        result = basis;

        for i = 1:3
            result(:,2*i-1,:) = n1(1)*basis(:,2*i-1,:) + n2(1)*basis(:,2*i,:);
            result(:,2*i,:) = n1(2)*basis(:,2*i-1,:) + n2(2)*basis(:,2*i,:);
            result(:,6+i,:) = s(i)*basis(:,6+i,:);
        end

    elseif dx==1 && dy==0
        basis_ref_dxh = basis_reference_2D(coords_h,basis_type,1,0);
        basis_ref_dyh = basis_reference_2D(coords_h,basis_type,0,1);
        basis_ref_dx = basis_ref_dxh*(y3-y1)/J0+basis_ref_dyh*(y1-y2)/J0;

        basis_ref_dx = reshape(basis_ref_dx,2,num_basis*num_pts);
        basis_dx = (1/J0)*J*basis_ref_dx;
        basis_dx = reshape(basis_dx,2,num_basis,num_pts);

        result = basis_dx;

        for i = 1:3
            result(:,2*i-1,:) = n1(1)*basis_dx(:,2*i-1,:) + n2(1)*basis_dx(:,2*i,:);
            result(:,2*i,:) = n1(2)*basis_dx(:,2*i-1,:) + n2(2)*basis_dx(:,2*i,:);
            result(:,6+i,:) = s(i)*basis_dx(:,6+i,:);
        end

    elseif dx==0 && dy==1
        basis_ref_dxh = basis_reference_2D(coords_h,basis_type,1,0);
        basis_ref_dyh = basis_reference_2D(coords_h,basis_type,0,1);
        basis_ref_dy =  basis_ref_dxh*(x1-x3)/J0+basis_ref_dyh*(x2-x1)/J0;

        basis_ref_dy = reshape(basis_ref_dy,2,num_basis*num_pts);
        basis_dy = (1/J0)*J*basis_ref_dy;
        basis_dy = reshape(basis_dy,2,num_basis,num_pts);

        result = basis_dy;

        for i = 1:3
            result(:,2*i-1,:) = n1(1)*basis_dy(:,2*i-1,:) + n2(1)*basis_dy(:,2*i,:);
            result(:,2*i,:) = n1(2)*basis_dy(:,2*i-1,:) + n2(2)*basis_dy(:,2*i,:);
            result(:,6+i,:) = s(i)*basis_dy(:,6+i,:);
        end
        
    else
        fprintf('cannot handle this type of derivative!');
    end

elseif strcmp(basis_type,'RT0')

    J = [x2-x1,x3-x1;y2-y1,y3-y1];
    coords_h = J\([x;y]-vertices(:,1));
    J0 = det(J);

    s = sign(mesh.TE(:,i_elem));
    arc = mesh.s(abs(mesh.TE(:,i_elem)));
    A = diag(s.*arc');
    if dx==0 && dy==0
        basis_ref = basis_reference_2D(coords_h,basis_type,dx,dy);

        basis_ref = reshape(basis_ref,2,num_basis*num_pts);
        basis = (1/J0)*J*basis_ref;
        basis = reshape(basis,2,num_basis,num_pts);

        result(1,:,:) = A*reshape(basis(1,:,:),num_basis,num_pts);
        result(2,:,:) = A*reshape(basis(2,:,:),num_basis,num_pts);

    elseif dx==1 && dy==0
        basis_ref_dxh = basis_reference_2D(coords_h,basis_type,1,0);
        basis_ref_dyh = basis_reference_2D(coords_h,basis_type,0,1);
        basis_ref_dx = basis_ref_dxh*(y3-y1)/J0+basis_ref_dyh*(y1-y2)/J0;

        basis_ref_dx = reshape(basis_ref_dx,2,num_basis*num_pts);
        basis_dx = (1/J0)*J*basis_ref_dx;
        basis_dx = reshape(basis_dx,2,num_basis,num_pts);

        result(1,:,:) = A*reshape(basis_dx(1,:,:),num_basis,num_pts);
        result(2,:,:) = A*reshape(basis_dx(2,:,:),num_basis,num_pts);

    elseif dx==0 && dy==1
        basis_ref_dxh = basis_reference_2D(coords_h,basis_type,1,0);
        basis_ref_dyh = basis_reference_2D(coords_h,basis_type,0,1);
        basis_ref_dy =  basis_ref_dxh*(x1-x3)/J0+basis_ref_dyh*(x2-x1)/J0;

        basis_ref_dy = reshape(basis_ref_dy,2,num_basis*num_pts);
        basis_dy = (1/J0)*J*basis_ref_dy;
        basis_dy = reshape(basis_dy,2,num_basis,num_pts);

        result(1,:,:) = A*reshape(basis_dy(1,:,:),num_basis,num_pts);
        result(2,:,:) = A*reshape(basis_dy(2,:,:),num_basis,num_pts);
    else
        fprintf('cannot handle this type of derivative!');
    end

elseif strcmp(basis_type,'MTW')
    J = [x2-x1,x3-x1;y2-y1,y3-y1];
    coords_h = J\([x;y]-vertices(:,1));
    J0 = det(J);

    s = sign(mesh.TE(:,i_elem));
    i_edge = abs(mesh.TE(:,i_elem));
    S = diag([s(1),s(1),s(2),s(2),s(3),s(3),s(1),s(2),s(3)]);

    G1 = [mesh.N(:,i_edge(1)),mesh.tau(:,i_edge(1))];
    G2 = [mesh.N(:,i_edge(2)),mesh.tau(:,i_edge(2))];
    G3 = [mesh.N(:,i_edge(3)),mesh.tau(:,i_edge(3))];

    ab1 = G1'*(J'*J)*mesh.tau(:,i_edge(1))/J0;
    ab2 = G2'*(J'*J)*mesh.tau(:,i_edge(2))/J0;
    ab3 = G3'*(J'*J)*mesh.tau(:,i_edge(3))/J0;

    A = [1,0,0,0,0,0,0,0,0;
        0,1,0,0,0,0,0,0,0;
        0,0,1,0,0,0,0,0,0;
        0,0,0,1,0,0,0,0,0;
        0,0,0,0,1,0,0,0,0;
        0,0,0,0,0,1,0,0,0;
        -ab1(1)/ab1(2),-ab1(1)/ab1(2),0,0,0,0,1/ab1(2),0,0;
        0,0,-ab2(1)/ab2(2),-ab2(1)/ab2(2),0,0,0,1/ab2(2),0;
        0,0,0,0,-ab3(1)/ab3(2),-ab3(1)/ab3(2),0,0,1/ab3(2)];
    arc = mesh.s(abs(mesh.TE(:,i_elem)));
    ARC = diag([arc(1),arc(1),arc(2),arc(2),arc(3),arc(3),arc(1),arc(2),arc(3)]);
    A = ARC*A'*S;

    if dx==0 && dy==0
        basis_ref = basis_reference_2D(coords_h,basis_type,dx,dy);

        basis_ref = reshape(basis_ref,2,num_basis*num_pts);
        basis = (1/J0)*J*basis_ref;
        basis = reshape(basis,2,num_basis,num_pts);

        result(1,:,:) = A*reshape(basis(1,:,:),num_basis,num_pts);
        result(2,:,:) = A*reshape(basis(2,:,:),num_basis,num_pts);

    elseif dx==1 && dy==0
        basis_ref_dxh = basis_reference_2D(coords_h,basis_type,1,0);
        basis_ref_dyh = basis_reference_2D(coords_h,basis_type,0,1);
        basis_ref_dx = basis_ref_dxh*(y3-y1)/J0+basis_ref_dyh*(y1-y2)/J0;

        basis_ref_dx = reshape(basis_ref_dx,2,num_basis*num_pts);
        basis_dx = (1/J0)*J*basis_ref_dx;
        basis_dx = reshape(basis_dx,2,num_basis,num_pts);

        result(1,:,:) = A*reshape(basis_dx(1,:,:),num_basis,num_pts);
        result(2,:,:) = A*reshape(basis_dx(2,:,:),num_basis,num_pts);

    elseif dx==0 && dy==1
        basis_ref_dxh = basis_reference_2D(coords_h,basis_type,1,0);
        basis_ref_dyh = basis_reference_2D(coords_h,basis_type,0,1);
        basis_ref_dy =  basis_ref_dxh*(x1-x3)/J0+basis_ref_dyh*(x2-x1)/J0;

        basis_ref_dy = reshape(basis_ref_dy,2,num_basis*num_pts);
        basis_dy = (1/J0)*J*basis_ref_dy;
        basis_dy = reshape(basis_dy,2,num_basis,num_pts);

        result(1,:,:) = A*reshape(basis_dy(1,:,:),num_basis,num_pts);
        result(2,:,:) = A*reshape(basis_dy(2,:,:),num_basis,num_pts);
    else
        fprintf('cannot handle this type of derivative!');
    end

elseif strcmp(basis_type,'RT0-MTW')
    J = [x2-x1,x3-x1;y2-y1,y3-y1];
    coords_h = J\([x;y]-vertices(:,1));
    J0 = det(J);

    s = sign(mesh.TE(:,i_elem));
    i_edge = abs(mesh.TE(:,i_elem));
    S = diag([s(1),s(2),s(3),s(1),s(2),s(3)]);

    G1 = [mesh.N(:,i_edge(1)),mesh.tau(:,i_edge(1))];
    G2 = [mesh.N(:,i_edge(2)),mesh.tau(:,i_edge(2))];
    G3 = [mesh.N(:,i_edge(3)),mesh.tau(:,i_edge(3))];

    ab1 = G1'*(J'*J)*mesh.tau(:,i_edge(1))/J0;
    ab2 = G2'*(J'*J)*mesh.tau(:,i_edge(2))/J0;
    ab3 = G3'*(J'*J)*mesh.tau(:,i_edge(3))/J0;

    A = [1,0,0,0,0,0;
        0,1,0,0,0,0;
        0,0,1,0,0,0;
        -ab1(1)/ab1(2),0,0,1/ab1(2),0,0;
        0,-ab2(1)/ab2(2),0,0,1/ab2(2),0;
        0,0,-ab3(1)/ab3(2),0,0,1/ab3(2)];
    arc = mesh.s(abs(mesh.TE(:,i_elem)));
    ARC = diag([arc(1),arc(2),arc(3),arc(1),arc(2),arc(3)]);
    A = ARC*A'*S;

    if dx==0 && dy==0
        basis_ref = basis_reference_2D(coords_h,basis_type,dx,dy);

        basis_ref = reshape(basis_ref,2,num_basis*num_pts);
        basis = (1/J0)*J*basis_ref;
        basis = reshape(basis,2,num_basis,num_pts);

        result(1,:,:) = A*reshape(basis(1,:,:),num_basis,num_pts);
        result(2,:,:) = A*reshape(basis(2,:,:),num_basis,num_pts);

    elseif dx==1 && dy==0
        basis_ref_dxh = basis_reference_2D(coords_h,basis_type,1,0);
        basis_ref_dyh = basis_reference_2D(coords_h,basis_type,0,1);
        basis_ref_dx = basis_ref_dxh*(y3-y1)/J0+basis_ref_dyh*(y1-y2)/J0;

        basis_ref_dx = reshape(basis_ref_dx,2,num_basis*num_pts);
        basis_dx = (1/J0)*J*basis_ref_dx;
        basis_dx = reshape(basis_dx,2,num_basis,num_pts);

        result(1,:,:) = A*reshape(basis_dx(1,:,:),num_basis,num_pts);
        result(2,:,:) = A*reshape(basis_dx(2,:,:),num_basis,num_pts);

    elseif dx==0 && dy==1
        basis_ref_dxh = basis_reference_2D(coords_h,basis_type,1,0);
        basis_ref_dyh = basis_reference_2D(coords_h,basis_type,0,1);
        basis_ref_dy =  basis_ref_dxh*(x1-x3)/J0+basis_ref_dyh*(x2-x1)/J0;

        basis_ref_dy = reshape(basis_ref_dy,2,num_basis*num_pts);
        basis_dy = (1/J0)*J*basis_ref_dy;
        basis_dy = reshape(basis_dy,2,num_basis,num_pts);

        result(1,:,:) = A*reshape(basis_dy(1,:,:),num_basis,num_pts);
        result(2,:,:) = A*reshape(basis_dy(2,:,:),num_basis,num_pts);
    else
        fprintf('cannot handle this type of derivative!');
    end

elseif strcmp(basis_type,'Ned1-1')
    J = [x2-x1,x3-x1;y2-y1,y3-y1];
    coords_h = J\([x;y]-vertices(:,1));
    J0 = det(J);

    s = sign(mesh.TE(:,i_elem));
    S = diag([s(1),s(2),s(3)]);

    arc = mesh.s(abs(mesh.TE(:,i_elem)));
    ARC = diag([arc(1),arc(2),arc(3)]);
    A = ARC*S;

    if dx==0 && dy==0
        basis_ref = basis_reference_2D(coords_h,basis_type,dx,dy);

        basis_ref = reshape(basis_ref,2,num_basis*num_pts);
        basis = J'\basis_ref;
        basis = reshape(basis,2,num_basis,num_pts);

        result(1,:,:) = A*reshape(basis(1,:,:),num_basis,num_pts);
        result(2,:,:) = A*reshape(basis(2,:,:),num_basis,num_pts);

    elseif dx==1 && dy==0
        basis_ref_dxh = basis_reference_2D(coords_h,basis_type,1,0);
        basis_ref_dyh = basis_reference_2D(coords_h,basis_type,0,1);
        basis_ref_dx = basis_ref_dxh*(y3-y1)/J0+basis_ref_dyh*(y1-y2)/J0;

        basis_ref_dx = reshape(basis_ref_dx,2,num_basis*num_pts);
        basis_dx = J'\basis_ref_dx;
        basis_dx = reshape(basis_dx,2,num_basis,num_pts);

        result(1,:,:) = A*reshape(basis_dx(1,:,:),num_basis,num_pts);
        result(2,:,:) = A*reshape(basis_dx(2,:,:),num_basis,num_pts);

    elseif dx==0 && dy==1
        basis_ref_dxh = basis_reference_2D(coords_h,basis_type,1,0);
        basis_ref_dyh = basis_reference_2D(coords_h,basis_type,0,1);
        basis_ref_dy =  basis_ref_dxh*(x1-x3)/J0+basis_ref_dyh*(x2-x1)/J0;

        basis_ref_dy = reshape(basis_ref_dy,2,num_basis*num_pts);
        basis_dy = J'\basis_ref_dy;
        basis_dy = reshape(basis_dy,2,num_basis,num_pts);

        result(1,:,:) = A*reshape(basis_dy(1,:,:),num_basis,num_pts);
        result(2,:,:) = A*reshape(basis_dy(2,:,:),num_basis,num_pts);
    else
        fprintf('cannot handle this type of derivative!');
    end
elseif strcmp(basis_type,'Ned1-2')
    J = [x2-x1,x3-x1;y2-y1,y3-y1];
    coords_h = J\([x;y]-vertices(:,1));
    J0 = det(J);

    s = sign(mesh.TE(:,i_elem));
    S = diag([s(1),s(1),s(2),s(2),s(3),s(3),sqrt(J0),sqrt(J0)]);

    arc = mesh.s(abs(mesh.TE(:,i_elem)));
    ARC = diag([arc(1),arc(1),arc(2),arc(2),arc(3),arc(3),1,1]);
    A = ARC*S;

    if dx==0 && dy==0
        basis_ref = basis_reference_2D(coords_h,basis_type,dx,dy);

        basis_ref = reshape(basis_ref,2,num_basis*num_pts);
        basis = J'\basis_ref;
        basis = reshape(basis,2,num_basis,num_pts);

        result(1,:,:) = A*reshape(basis(1,:,:),num_basis,num_pts);
        result(2,:,:) = A*reshape(basis(2,:,:),num_basis,num_pts);

    elseif dx==1 && dy==0
        basis_ref_dxh = basis_reference_2D(coords_h,basis_type,1,0);
        basis_ref_dyh = basis_reference_2D(coords_h,basis_type,0,1);
        basis_ref_dx = basis_ref_dxh*(y3-y1)/J0+basis_ref_dyh*(y1-y2)/J0;

        basis_ref_dx = reshape(basis_ref_dx,2,num_basis*num_pts);
        basis_dx = J'\basis_ref_dx;
        basis_dx = reshape(basis_dx,2,num_basis,num_pts);

        result(1,:,:) = A*reshape(basis_dx(1,:,:),num_basis,num_pts);
        result(2,:,:) = A*reshape(basis_dx(2,:,:),num_basis,num_pts);

    elseif dx==0 && dy==1
        basis_ref_dxh = basis_reference_2D(coords_h,basis_type,1,0);
        basis_ref_dyh = basis_reference_2D(coords_h,basis_type,0,1);
        basis_ref_dy =  basis_ref_dxh*(x1-x3)/J0+basis_ref_dyh*(x2-x1)/J0;

        basis_ref_dy = reshape(basis_ref_dy,2,num_basis*num_pts);
        basis_dy = J'\basis_ref_dy;
        basis_dy = reshape(basis_dy,2,num_basis,num_pts);

        result(1,:,:) = A*reshape(basis_dy(1,:,:),num_basis,num_pts);
        result(2,:,:) = A*reshape(basis_dy(2,:,:),num_basis,num_pts);
    else
        fprintf('cannot handle this type of derivative!');
    end
elseif strcmp(basis_type,'Ned2-1')
    J = [x2-x1,x3-x1;y2-y1,y3-y1];
    coords_h = J\([x;y]-vertices(:,1));
    J0 = det(J);

    s = sign(mesh.TE(:,i_elem));
    S = diag([s(1),s(1),s(2),s(2),s(3),s(3)]);

    arc = mesh.s(abs(mesh.TE(:,i_elem)));
    ARC = diag([arc(1),arc(1),arc(2),arc(2),arc(3),arc(3)]);
    A = ARC*S;

    if dx==0 && dy==0
        basis_ref = basis_reference_2D(coords_h,basis_type,dx,dy);

        basis_ref = reshape(basis_ref,2,num_basis*num_pts);
        basis = J'\basis_ref;
        basis = reshape(basis,2,num_basis,num_pts);

        result(1,:,:) = A*reshape(basis(1,:,:),num_basis,num_pts);
        result(2,:,:) = A*reshape(basis(2,:,:),num_basis,num_pts);

    elseif dx==1 && dy==0
        basis_ref_dxh = basis_reference_2D(coords_h,basis_type,1,0);
        basis_ref_dyh = basis_reference_2D(coords_h,basis_type,0,1);
        basis_ref_dx = basis_ref_dxh*(y3-y1)/J0+basis_ref_dyh*(y1-y2)/J0;

        basis_ref_dx = reshape(basis_ref_dx,2,num_basis*num_pts);
        basis_dx = J'\basis_ref_dx;
        basis_dx = reshape(basis_dx,2,num_basis,num_pts);

        result(1,:,:) = A*reshape(basis_dx(1,:,:),num_basis,num_pts);
        result(2,:,:) = A*reshape(basis_dx(2,:,:),num_basis,num_pts);

    elseif dx==0 && dy==1
        basis_ref_dxh = basis_reference_2D(coords_h,basis_type,1,0);
        basis_ref_dyh = basis_reference_2D(coords_h,basis_type,0,1);
        basis_ref_dy =  basis_ref_dxh*(x1-x3)/J0+basis_ref_dyh*(x2-x1)/J0;

        basis_ref_dy = reshape(basis_ref_dy,2,num_basis*num_pts);
        basis_dy = J'\basis_ref_dy;
        basis_dy = reshape(basis_dy,2,num_basis,num_pts);

        result(1,:,:) = A*reshape(basis_dy(1,:,:),num_basis,num_pts);
        result(2,:,:) = A*reshape(basis_dy(2,:,:),num_basis,num_pts);
    else
        fprintf('cannot handle this type of derivative!');
    end
end

