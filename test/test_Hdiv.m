function test_Hdiv
% Regression checks on the unit square; run with the repository on the path.
[T,P] = mesh_square_triangle(0,1,0,1,2,2);
mesh = mesh_topo(T,P);
space = mesh_FE('P1',2,mesh,3,0,[]);
tol = 1e-12;
zero = zeros(space.N_node,1);

% Opposite diagonal derivatives cancel in the divergence.
u = interpolate(@(x) [x(1,:);-x(2,:)],mesh,space);
assert(compute_norm(u,mesh,space,'Hdiv-semi',3) < tol);
assert(compute_error(u,@exact_field,mesh,space,'Hdiv',3) < tol);
assert(abs(compute_error(zero,@exact_field,mesh,space,'Hdiv',3)-sqrt(2/3)) < tol);

% Nonzero divergence tests the cross term and the full error norm.
v = interpolate(@(x) [x(1,:);x(2,:)],mesh,space);
assert(abs(compute_norm(v,mesh,space,'Hdiv-semi',3)-2) < tol);
assert(abs(compute_error(v,@zero_field,mesh,space,'Hdiv',3)-sqrt(14/3)) < tol);

% A constant error must retain its L2 contribution.
c = interpolate(@(x) repmat([1;2],1,size(x,2)),mesh,space);
assert(abs(compute_error(c,@zero_field,mesh,space,'Hdiv',3)-sqrt(5)) < tol);
fprintf('Hdiv regression checks passed.\n');
end

function value = exact_field(x,dx,dy)
value = zeros(2,size(x,2));
if dx==0 && dy==0
    value = [x(1,:);-x(2,:)];
elseif dx==1 && dy==0
    value(1,:) = 1;
elseif dx==0 && dy==1
    value(2,:) = -1;
end
end

function value = zero_field(x,~,~)
value = zeros(2,size(x,2));
end
