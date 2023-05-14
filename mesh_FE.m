function mesh_FE = mesh_FE(basis_type,dim,mesh,Gauss_type_2D,Gauss_type_1D,bndy_number)
% generate data for a certain FE space.
% 
% input:
%   basis_type: type of the FE space, 'P0', 'P1', 'P1dc', 'DGP1', ...
%   dim: dimension of the FE space
%   mesh: mesh data
%   Gauss_type_2D: type of the 2D Gauss quadrature rule
%   Gauss_type_1D: type of the 1D Gauss quadrature rule
%   bndy_number: number of boundary edges to be used in the computation of the boundary integral, [] means all edges
%
% output:
%   mesh_FE: mesh data for the FE space
%
% Attributes:
%   mesh_FE.dim: dimension of the FE space
%   mesh_FE.N_lb: number of local basis functions
%   mesh_FE.P: data of the DOFs, coordinates of the pts or number of edges, etc.
%   mesh_FE.N_node: number of DOFS
%   mesh_FE.basis_type: type of the FE space
%   mesh_FE.T: data of the elements, DOFs of the elements
%   mesh_FE.Dbndynodes: DOFs on the Dirichlet boundary
%   mesh_FE.basis: cached basis functions
%   mesh_FE.basis_dx: cached basis functions' derivatives
%   mesh_FE.basis_dy: cached basis functions' derivatives


if strcmp(basis_type,'P0')

    mesh_FE.dim = dim;
    mesh_FE.N_lb = 1;
    P = zeros(2,mesh.N_elem);
    for i = 1:mesh.N_elem
        vertices = mesh.P(:,mesh.T(:,i));
        P(:,i) = mean(vertices,2);
    end
    mesh_FE.P = repmat(P,1,dim);
    mesh_FE.N_node = mesh.N_elem*dim;
    mesh_FE.basis_type = basis_type;
    mesh_FE.idx = kron(1:dim,ones(1,mesh.N_elem));
    mesh_FE.T = repmat(1:mesh.N_elem,dim,1);

    for i = 1:dim
        mesh_FE.T(i,:) =  (1:mesh.N_elem)+(i-1)*mesh.N_elem;
    end

    flag = mesh.E(1,:)==0;
    Dbndynodes = unique(mesh.E(2,flag));
    N_bndy = numel(Dbndynodes);
    mesh_FE.Dbndynodes = zeros(1,dim*N_bndy);
    for i = 1:dim
        mesh_FE.Dbndynodes(1+(i-1)*N_bndy:i*N_bndy) = Dbndynodes + (i-1)*mesh.N_elem;
    end

elseif strcmp(basis_type,'P1')

    mesh_FE.dim = dim;
    mesh_FE.N_lb = 3;
    mesh_FE.P = repmat(mesh.P,1,dim);
    mesh_FE.N_node = mesh.N_node*dim;
    mesh_FE.basis_type = basis_type;
    mesh_FE.idx = kron(1:dim,ones(1,mesh.N_node));
    mesh_FE.T = zeros(3*dim,mesh.N_elem);
    for i = 1:dim
        mesh_FE.T(3*(i-1)+1:3*i,:) = mesh.T + (i-1)*mesh.N_node;
    end

    flag = mesh.E(1,:)==0;
    Dbndynodes = unique([mesh.E(3,flag),mesh.E(4,flag)]);
    N_bndy = numel(Dbndynodes);
    mesh_FE.Dbndynodes = zeros(1,dim*N_bndy);
    for i = 1:dim
        mesh_FE.Dbndynodes(1+(i-1)*N_bndy:i*N_bndy) = Dbndynodes + (i-1)*mesh.N_node;
    end

elseif strcmp(basis_type,'P1dc')

    mesh_FE.dim = dim;
    mesh_FE.N_lb = 3;
    mesh_FE.N_node = 3*mesh.N_elem*dim;
    P = zeros(2,3*mesh.N_elem);
    for i_elem = 1:mesh.N_elem
        P(:,3*(i_elem-1)+1:3*i_elem) = mesh.P(:,mesh.T(:,i_elem));
    end
    mesh_FE.P = repmat(P,1,dim);
    mesh_FE.basis_type = basis_type;
    mesh_FE.idx = kron(1:dim,ones(1,3*mesh.N_elem));
    mesh_FE.T = zeros(3*dim,mesh.N_elem);
    for i = 1:dim
        mesh_FE.T(3*i,:) = 3*(1:mesh.N_elem)+(i-1)*3*mesh.N_elem;
        mesh_FE.T(3*i-1,:) = 3*(1:mesh.N_elem)+(i-1)*3*mesh.N_elem-1;
        mesh_FE.T(3*i-2,:) = 3*(1:mesh.N_elem)+(i-1)*3*mesh.N_elem-2;
    end

    % boundary dofs
    flag = mesh.E(1,:)==0;
    tmp = 1:mesh.N_edge;
    Dbndy_edges = tmp(flag);
    N_Dbndyedges = sum(flag);
    Dbndynodes = zeros(1,2*N_Dbndyedges*dim);
    count = 0;
    for i_bndy = 1:N_Dbndyedges
        i_edge = Dbndy_edges(i_bndy);
        i_elem = mesh.E(2,i_edge);
        for i_dim = 1:dim
            count = count + 1;
            Dbndynodes(1,count) = 3*(i_elem-1) + (i_dim-1)*3*mesh.N_elem + mesh.ET(2,i_edge);
            count = count + 1;
            Dbndynodes(1,count) = 3*(i_elem-1) + (i_dim-1)*3*mesh.N_elem + mod(mesh.ET(2,i_edge),3)+1;
        end
    end
    mesh_FE.Dbndynodes = unique(Dbndynodes);

elseif strcmp(basis_type,'DGP1')
    % DGP1 is a discontinuous P1 element, with degrees of freedoms on 3 edges of each triangle

    mesh_FE.dim = dim;
    mesh_FE.N_lb = 3;
    mesh_FE.N_node = 3*mesh.N_elem*dim;
    P = zeros(2,3*mesh.N_elem);
    for i_elem = 1:mesh.N_elem
        vertices = mesh.P(:,mesh.T(:,i_elem));
        % calculate the midpoints of the edges
        P(:,3*(i_elem-1)+1) = mean(vertices(:,[1,2]),2);
        P(:,3*(i_elem-1)+2) = mean(vertices(:,[2,3]),2);
        P(:,3*(i_elem-1)+3) = mean(vertices(:,[3,1]),2);
    end
    mesh_FE.P = repmat(P,1,dim);
    mesh_FE.basis_type = basis_type;
    mesh_FE.idx = kron(1:dim,ones(1,3*mesh.N_elem));
    mesh_FE.T = zeros(3*dim,mesh.N_elem);
    for i_dim = 1:dim
        mesh_FE.T(3*i_dim-2,:) = 3*(1:mesh.N_elem)+(i_dim-1)*3*mesh.N_elem-2;
        mesh_FE.T(3*i_dim-1,:) = 3*(1:mesh.N_elem)+(i_dim-1)*3*mesh.N_elem-1;
        mesh_FE.T(3*i_dim,:) = 3*(1:mesh.N_elem)+(i_dim-1)*3*mesh.N_elem;
    end

    % boundary dofs
    flag = mesh.E(1,:)==0;
    tmp = 1:mesh.N_edge;
    Dbndy_edges = tmp(flag);
    N_Dbndyedges = sum(flag);
    Dbndynodes = zeros(1,N_Dbndyedges*dim);
    count = 0;
    for i_bndy = 1:N_Dbndyedges
        i_edge = Dbndy_edges(i_bndy);
        i_elem = mesh.E(2,i_edge);
        for i_dim = 1:dim
            count = count + 1;
            Dbndynodes(1,count) = 3*(i_elem-1) + (i_dim-1)*3*mesh.N_elem + mesh.ET(2,i_edge);
        end
    end
    mesh_FE.Dbndynodes = unique(Dbndynodes);

 
elseif strcmp(basis_type,'P2')

    mesh_FE.dim = dim;
    mesh_FE.N_lb = 6;
    mesh_FE.basis_type = basis_type;

    N_node = mesh.N_node + mesh.N_edge;
    mesh_FE.N_node = N_node*dim;
    mesh_FE.idx = kron(1:dim,ones(1,N_node));

    P = zeros(2,N_node);
    P(:,1:mesh.N_node) = mesh.P;

    for i = 1:mesh.N_edge
        edge = mesh.E([3,4],i);
        vertices = mesh.P(:,edge);
        P(:,i+mesh.N_node) = mean(vertices,2);
    end
    mesh_FE.P = repmat(P,1,dim);

    T = zeros(6,mesh.N_elem);
    T(1:3,:) = mesh.T;
    T([4,5,6],:) = mesh.N_node+abs(mesh.TE([1,2,3],:));
    mesh_FE.T = zeros(6*dim,mesh.N_elem);
    for i = 1:dim
        mesh_FE.T(6*(i-1)+1:6*i,:) = T + (i-1)*N_node;
    end

    flag = mesh.E(1,:)==0;
    newnodes = mesh.N_node+1:mesh.N_node+mesh.N_edge;
    Dbndynodes = unique([mesh.E(3,flag),mesh.E(4,flag),newnodes(flag)]);
    N_bndy = numel(Dbndynodes);
    mesh_FE.Dbndynodes = zeros(1,dim*N_bndy);
    for i = 1:dim
        mesh_FE.Dbndynodes(1+(i-1)*N_bndy:i*N_bndy) = Dbndynodes + (i-1)*N_node;
    end

elseif strcmp(basis_type,'bubbleP1')

    mesh_FE.dim = dim;
    mesh_FE.N_lb = 4;
    P_center = zeros(2,mesh.N_elem);
    for i_elem = 1:mesh.N_elem
        P_center(:,i_elem) = mean(mesh.P(:,mesh.T(:,i_elem)),2);
    end
    P = [mesh.P P_center];
    mesh_FE.P = repmat(P,1,dim);
    mesh_FE.N_node = (mesh.N_node+mesh.N_elem)*dim;
    mesh_FE.basis_type = basis_type;
    mesh_FE.idx = kron(1:dim,ones(1,mesh.N_node+mesh.N_elem));
    mesh_FE.T = zeros(4*dim,mesh.N_elem);
    for i = 1:dim
        mesh_FE.T(4*(i-1)+1:4*i,:) = [mesh.T;(1:mesh.N_elem)+mesh.N_node] + (i-1)*(mesh.N_node+mesh.N_elem);
    end

    flag = mesh.E(1,:)==0;
    Dbndynodes = unique([mesh.E(3,flag),mesh.E(4,flag)]);
    N_bndy = numel(Dbndynodes);
    mesh_FE.Dbndynodes = zeros(1,dim*N_bndy);
    for i = 1:dim
        mesh_FE.Dbndynodes(1+(i-1)*N_bndy:i*N_bndy) = Dbndynodes + (i-1)*(mesh.N_node+mesh.N_elem);
    end

elseif strcmp(basis_type,'CR') || strcmp(basis_type,'CR-P0')
    mesh_FE.dim = dim;
    mesh_FE.N_lb = 3;

    P = (mesh.P(:,mesh.E(3,:))+mesh.P(:,mesh.E(4,:)))/2;
    mesh_FE.P = repmat(P,1,dim);
    mesh_FE.N_node = mesh.N_edge*dim;
    mesh_FE.basis_type = basis_type;
    mesh_FE.idx = kron(1:dim,ones(1,mesh.N_edge));
    mesh_FE.T = zeros(3*dim,mesh.N_elem);
    for i = 1:dim
        mesh_FE.T(3*(i-1)+1:3*i,:) = abs(mesh.TE) + (i-1)*mesh.N_edge;
    end

    flag = mesh.E(1,:)==0;
    tmp = 1:mesh.N_edge;
    Dbndynodes = tmp(flag);
    N_bndy = numel(Dbndynodes);
    mesh_FE.Dbndynodes = zeros(1,dim*N_bndy);
    for i = 1:dim
        mesh_FE.Dbndynodes(1+(i-1)*N_bndy:i*N_bndy) = Dbndynodes + (i-1)*mesh.N_edge;
    end

%     flag = mesh.E(1,:)==-1;
%     bndynodes1 = tmp(flag);
%     mesh_FE.Dbndynodes = [mesh_FE.Dbndynodes bndynodes1];
% 
%     flag = mesh.E(1,:)==-2;
%     bndynodes2 = tmp(flag);
%     mesh_FE.Dbndynodes = [mesh_FE.Dbndynodes mesh.N_edge+bndynodes2];
    

elseif strcmp(basis_type,'CR-RT0')
    mesh_FE.dim = 2;
    mesh_FE.N_lb = 6;

    P = (mesh.P(:,mesh.E(3,:))+mesh.P(:,mesh.E(4,:)))/2;
    mesh_FE.P = repmat(P,1,dim);
    mesh_FE.N_node = mesh.N_edge*dim;
    mesh_FE.basis_type = basis_type;
    mesh_FE.idx = kron(1:dim,ones(1,mesh.N_edge));
    mesh_FE.T = zeros(6*dim,mesh.N_elem);
    for i = 1:2
        mesh_FE.T(3*(i-1)+1:3*i,:) = abs(mesh.TE) + (i-1)*mesh.N_edge;
    end
    mesh_FE.T(7:12,:) = mesh_FE.T(1:6,:);

    flag = mesh.E(1,:)==0;
    tmp = 1:mesh.N_edge;
    Dbndynodes = tmp(flag);
    N_bndy = numel(Dbndynodes);
    mesh_FE.Dbndynodes = zeros(1,dim*N_bndy);
    for i = 1:dim
        mesh_FE.Dbndynodes(1+(i-1)*N_bndy:i*N_bndy) = Dbndynodes + (i-1)*mesh.N_edge;
    end

elseif strcmp(basis_type,'RT0')

    mesh_FE.dim = dim;
    mesh_FE.N_node = mesh.N_edge;
    mesh_FE.P = repmat(1:mesh.N_edge,2,1);

    flag = mesh.E(1,:)==0;
    tmp = 1:mesh.N_edge;
    mesh_FE.Dbndynodes = tmp(flag);
    mesh_FE.basis_type = basis_type;
    mesh_FE.N_lb = 3;
    mesh_FE.T = repmat(abs(mesh.TE),2,1);

elseif strcmp(basis_type,'BR') || strcmp(basis_type,'BR-RT0')
    mesh_FE.N_lb = 9;
    mesh_FE.dim = dim;
    mesh_FE.N_node = 2*mesh.N_node+mesh.N_edge;
    P = [mesh.P,mesh.P,repmat(1:mesh.N_edge,2,1)];
    mesh_FE.P = P;
    T = zeros(9,mesh.N_elem);
    T(1:3,:) = mesh.T;
    T(4:6,:) = mesh.T+mesh.N_node;
    T(7:9,:) = abs(mesh.TE)+2*mesh.N_node;
    mesh_FE.T = repmat(T,dim,1);
    mesh_FE.basis_type = basis_type;

    flag = mesh.E(1,:)==0;
    tmp = 1:mesh.N_edge;
    Dbndynodes = unique([mesh.E(3,flag),mesh.E(4,flag),mesh.E(3,flag)+mesh.N_node,...
        mesh.E(4,flag)+mesh.N_node, tmp(flag)+2*mesh.N_node]);
    mesh_FE.Dbndynodes = Dbndynodes;

elseif strcmp(basis_type,'MTW')
    mesh_FE.N_lb = 9;
    mesh_FE.dim = dim;
    mesh_FE.N_node = 3*mesh.N_edge;
    T = zeros(9,mesh.N_elem);
    S = sign(mesh.TE);
    T([1,3,5],:) = (1/2-S/2)*mesh.N_edge+abs(mesh.TE);
    T([2,4,6],:) = (1/2+S/2)*mesh.N_edge+abs(mesh.TE);
    T([7,8,9],:) = 2*mesh.N_edge+abs(mesh.TE);
    mesh_FE.T = repmat(T,dim,1);
    P = [1:mesh.N_edge,1:mesh.N_edge,1:mesh.N_edge];
    mesh_FE.P = P;

    mesh_FE.basis_type = basis_type;
    flag = mesh.E(1,:)==0;
    tmp = 1:mesh.N_edge;
    Dbndynodes = [tmp(flag),tmp(flag)+mesh.N_edge,tmp(flag)+mesh.N_edge*2];
    mesh_FE.Dbndynodes = Dbndynodes;
elseif strcmp(basis_type,'RT0-MTW')
    mesh_FE.N_lb = 6;
    mesh_FE.dim = dim;
    mesh_FE.N_node = 2*mesh.N_edge;
    T = zeros(6,mesh.N_elem);
    T([1,2,3],:) = abs(mesh.TE);
    T([4,5,6],:) = mesh.N_edge+abs(mesh.TE);
    mesh_FE.T = repmat(T,dim,1);
    P = [1:mesh.N_edge,1:mesh.N_edge];
    mesh_FE.P = P;

    mesh_FE.basis_type = basis_type;
    flag = mesh.E(1,:)==0;
    tmp = 1:mesh.N_edge;
    Dbndynodes = [tmp(flag),tmp(flag)+mesh.N_edge];
    mesh_FE.Dbndynodes = Dbndynodes;

elseif strcmp(basis_type,'Ned1-1')
    mesh_FE.N_lb = 3;
    mesh_FE.dim = dim;
    mesh_FE.N_node = mesh.N_edge;
    T = abs(mesh.TE);
    mesh_FE.T = repmat(T,dim,1);
    P = 1:mesh.N_edge;
    mesh_FE.P = P;
    mesh_FE.basis_type = basis_type;
    flag = mesh.E(1,:)==0;
    tmp = 1:mesh.N_edge;
    Dbndynodes = tmp(flag);
    mesh_FE.Dbndynodes = Dbndynodes;

elseif strcmp(basis_type,'Ned1-2')
    mesh_FE.N_lb = 8;
    mesh_FE.dim = dim;
    mesh_FE.N_node = 2*mesh.N_edge+2*mesh.N_elem;

    T = zeros(8,mesh.N_elem);
    S = sign(mesh.TE);
    T([1,3,5],:) = (1/2-S/2)*mesh.N_edge+abs(mesh.TE);
    T([2,4,6],:) = (1/2+S/2)*mesh.N_edge+abs(mesh.TE);
    T(7,:) = 2*mesh.N_edge+(1:mesh.N_elem);
    T(8,:) = 2*mesh.N_edge+mesh.N_elem+(1:mesh.N_elem);
    mesh_FE.T = repmat(T,dim,1);

    P = [1:mesh.N_edge,1:mesh.N_edge,1:mesh.N_elem,1:mesh.N_elem];
    mesh_FE.P = P;
    mesh_FE.basis_type = basis_type;

    flag = mesh.E(1,:)==0;
    tmp = 1:mesh.N_edge;
    Dbndynodes = [tmp(flag),mesh.N_edge+tmp(flag)];
    mesh_FE.Dbndynodes = Dbndynodes;

elseif strcmp(basis_type,'Ned2-1')
    mesh_FE.N_lb = 6;
    mesh_FE.dim = dim;
    mesh_FE.N_node = 2*mesh.N_edge;

    T = zeros(6,mesh.N_elem);
    S = sign(mesh.TE);
    T([1,3,5],:) = (1/2-S/2)*mesh.N_edge+abs(mesh.TE);
    T([2,4,6],:) = (1/2+S/2)*mesh.N_edge+abs(mesh.TE);
    mesh_FE.T = repmat(T,dim,1);

    P = [1:mesh.N_edge,1:mesh.N_edge];
    mesh_FE.P = P;
    mesh_FE.basis_type = basis_type;

    flag = mesh.E(1,:)==0;
    tmp = 1:mesh.N_edge;
    Dbndynodes = [tmp(flag),mesh.N_edge+tmp(flag)];
    mesh_FE.Dbndynodes = Dbndynodes;

else
    fprintf('cannot handle the basis type!\n');
end

[basis,basis_dx,basis_dy] = generate_basis_data_triangle(mesh,mesh_FE,Gauss_type_2D);
mesh_FE.basis = basis;
mesh_FE.basis_dx = basis_dx;
mesh_FE.basis_dy = basis_dy;

if Gauss_type_1D~=0
    [basis,basis_dx,basis_dy] = generate_basis_data_line(mesh,mesh_FE,Gauss_type_1D,bndy_number);
    mesh_FE.basis_line = basis;
    mesh_FE.basis_line_dx = basis_dx;
    mesh_FE.basis_line_dy = basis_dy;
end
end
