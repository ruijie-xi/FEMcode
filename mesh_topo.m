function mesh = mesh_topo(T,P,varargin)
% calculate all mesh information from given T,P matrix, including:
%   E: Edges. 
% It's a 4*N_edge matrix, where N_edge is the number of all edges.
% The first and second rows stores the index of elements lying on both
% side of the edge. If edge j is on the boundary, then the first row of 
% the j th column will store the type of this boundary (0 for Dirichlet, -1
% for Neumann and -2 for Robin). The third and the fourth row stores the
% index of nodes on both ends of the edge. Note that each edge has an
% orientation, from the node in 3rd row to 4th row.
%   TE: Edges included in an element
% It's a 3(4)*N_elem matrix, where N_elem is the number of all elements. 
% For the i th elements, indices of all the three edges will be stored in
% the i th column. If the orientation of an edge in this element is not its
% intrinsic orientation, then the index will be negative.
%   N: unit normal vector of each edge.
% It's a 2*N_edge matrix.
% For the i th edges, its normal vector is stored in the i th column. 
% Direction of the normal vector is defined as turning the orientation of edge
% 90 degree clockwise.

%   AT: Adjacent elements for each element
% It's a 3(4)*N_elem matrix, where N_elem is the number of all elements.
% For the i th elements, indices of all the three adjacent elements will be
% stored in the i th column.

%   tau: unit tangential vector of each edge.
% It's a 2*N_edge matrix.

%   ET: Edges included in an element
% It's a 2*N_edge matrix, where N_edge is the number of all edges.
% For the i th edges, indices of the two elements lying on both side of the
% edge will be stored in the i th column. If the edge is on the boundary,
% then the second element will be 0.


N_bndytype = nargin-2;

[T,P] = orientation_consistence(T,P);

N_elem = size(T,2);
N_node = size(P,2);
N_le = size(T,1); % number of local edges, 3 for triangle, 4 for quadrilateral.

TE = zeros(N_le,N_elem);

% create matrix E
index_i = zeros(N_le*N_elem,1);
index_j = zeros(N_le*N_elem,1);
index_v = ones(N_le*N_elem,1);
count = 0;
for n = 1:N_elem
    for i = 1:N_le
        count = count+1;
        local_edge = T([i,mod(i,N_le)+1],n);
        if local_edge(1)>local_edge(2)
            local_edge = flipud(local_edge);
        end
        index_i(count) = local_edge(1);
        index_j(count) = local_edge(2);
    end
end
A = sparse(index_i,index_j,index_v,N_node,N_node);
[index_i,index_j,~] = find(A);
N_edge = nnz(A);
E = zeros(4,N_edge);
E([3,4],:) = [index_i';index_j'];

B = accumarray([index_i,index_j],(1:N_edge)',[],[],[],1);
ET = zeros(2,N_edge);

for n = 1:N_elem
    for i = 1:N_le
        local_edge = T([i,mod(i,N_le)+1],n);
        orientation = 1;
        if local_edge(1)>local_edge(2)
            local_edge = flipud(local_edge);
            orientation = -1;
        end
        idx = full(B(local_edge(1),local_edge(2)));
        TE(i,n) = orientation*idx;
        if E(2,idx)==0
            E(2,idx) = n;
            ET(2,idx) = i;
        else
            E(1,idx) = n;
            ET(1,idx) = i;
            if orientation==-1
                tmp = E(2,idx);
                E(2,idx) = E(1,idx);
                E(1,idx) = tmp;
                tmp = ET(2,idx);
                ET(2,idx) = ET(1,idx);
                ET(1,idx) = tmp;
            end
        end
    end
end

% adjacent elements for each element
AT = zeros(N_le,N_elem);
for i = 1:N_elem
    for j = 1:N_le
        idx = TE(j,i);
        if idx>0 % intrinsic orientation
            AT(j,i) = E(2,idx);
        else % non-intrinsic orientation
            AT(j,i) = E(1,-idx);
        end
    end
end

% normal vectors and tangential vectors.
N = zeros(2,N_edge);
tau = zeros(2,N_edge);
for i = 1:N_edge
    tau_e = P(:,E(4,i))-P(:,E(3,i));
    tau_e = tau_e/norm(tau_e);
    N_e = [tau_e(2);-tau_e(1)];
    N(:,i) = N_e;
    tau(:,i) = tau_e;
end

% Boundary type
bndy_flag = E(1,:)==0;
N_bndy = sum(bndy_flag);
tmp = 1:N_edge;
bndy_idx = tmp(bndy_flag);
ep = 1e-8;
for i = 1:N_bndy
    pt1 = E(3,bndy_idx(i));
    pt2 = E(4,bndy_idx(i));
    for i_bndytype = 1:N_bndytype
        fun = varargin{N_bndytype-i_bndytype+1};
        if abs(fun(P(:,pt1)))<=ep && abs(fun(P(:,pt2)))<=ep
            E(1,bndy_idx(i)) = -N_bndytype+i_bndytype;
        end
    end
end

%arclength
s = zeros(1,N_edge);
for i = 1:N_edge
    pt1 = E(3,i);
    pt2 = E(4,i);
    s(i) = norm(P(:,pt1)-P(:,pt2));
end

%hmax for triangulation.
if N_le==3
    hmax = 0;
    for i_elem = 1:N_elem
        vertices = P(:,T(:,i_elem));
        hmax_elem = max([norm(vertices(:,1)-vertices(:,2)),norm(vertices(:,2)-vertices(:,3)),norm(vertices(:,3)-vertices(:,1))]);
        hmax = max([hmax,hmax_elem]);
    end
end

% mesh area
area = 0;
min_area = 999999999999;
max_area = 0;
for i_elem = 1:N_elem
    vertices = P(:,T(:,i_elem));
    mat = [vertices(:,2)-vertices(:,1),vertices(:,3)-vertices(:,1)];
    elem_area = det(mat)/2;
    min_area = min([min_area,elem_area]);
    max_area = max([max_area,elem_area]);
    area = area + elem_area;
end

% generate boxes to make it easier to find element containing a point
xmin = min(P(1,:));
ymin = min(P(2,:));
xmax = max(P(1,:));
ymax = max(P(2,:));

Nx = ceil((xmax-xmin)/hmax);
Ny = ceil((ymax-ymin)/hmax);
hx = (xmax-xmin)/Nx;
hy = (ymax-ymin)/Ny;
boxes = cell(Nx,Ny);
for i_elem=1:N_elem
    vertices = P(:,T(:,i_elem));
    elemxmin = min(vertices(1,:));
    elemxmax = max(vertices(1,:));
    elemymin = min(vertices(2,:));
    elemymax = max(vertices(2,:));

    boxleft = ceil((elemxmin-xmin)/hx)-1;
    boxright = ceil((elemxmax-xmin)/hx)+1;
    boxbottom = ceil((elemymin-ymin)/hy)-1;
    boxtop = ceil((elemymax-ymin)/hy)+1;
    boxleft = min([max([1,boxleft]),Nx]);
    boxright = min([max([1,boxright]),Nx]);
    boxtop = min([max([1,boxtop]),Ny]);
    boxbottom = min([max([1,boxbottom]),Ny]);

    for i_box = boxleft:boxright
        for j_box = boxbottom:boxtop
            boxes{i_box,j_box} = [boxes{i_box,j_box} i_elem];
        end
    end
end

% neighbor matrix from element to node (N_elem x N_node)
neigh_i = zeros(N_le*N_elem,1);
neigh_j = zeros(N_le*N_elem,1);
neigh_v = ones(N_le*N_elem,1);
for i_elem = 1:N_elem
    neigh_i((i_elem-1)*N_le+1:i_elem*N_le) = i_elem;
    neigh_j((i_elem-1)*N_le+1:i_elem*N_le) = T(:,i_elem);
end
Mat_elem_node = sparse(neigh_i,neigh_j,neigh_v);

% neighbor matrix from element to edge (N_elem x N_edge)
neigh_i = zeros(N_le*N_elem,1);
neigh_j = zeros(N_le*N_elem,1);
neigh_v = ones(N_le*N_elem,1);
for i_elem = 1:N_elem
    neigh_i((i_elem-1)*N_le+1:i_elem*N_le) = i_elem;
    neigh_j((i_elem-1)*N_le+1:i_elem*N_le) = abs(TE(:,i_elem));
end
Mat_elem_edge = sparse(neigh_i,neigh_j,neigh_v);

mesh.T = T;
mesh.P = P;
mesh.E = E;
mesh.TE = TE;
mesh.ET = ET;
mesh.AT = AT;
mesh.N = N;
mesh.tau = tau;
mesh.N_elem = N_elem;
mesh.N_edge = N_edge;
mesh.N_node = N_node;
mesh.type = N_le;
mesh.s = s;
mesh.hmax = hmax;
mesh.area = area;
mesh.min_area = min_area;
mesh.max_area = max_area;
mesh.xmin = xmin;
mesh.xmax = xmax;
mesh.ymin = ymin;
mesh.ymax = ymax;
mesh.boxes = boxes;
mesh.Mat_elem_node = Mat_elem_node;
mesh.Mat_elem_edge = Mat_elem_edge;
end


