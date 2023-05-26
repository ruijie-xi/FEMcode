function result = basis_reference_2D(coords,basis_type,dx,dy)
% evaluate basis functions (or its derivatives) on reference element.
% coords: coordinates of the pt on which the basis function is evaluated.

x = coords(1,:);y = coords(2,:);
num_pts = size(coords,2);
one = ones(1,num_pts);
zero = zeros(1,num_pts);

if strcmp(basis_type,'P0')
    result = zeros(1,1,num_pts);
    if dx==0 && dy==0
        basis_ref = one;
    else
        basis_ref = zero;
    end
    result(1,:,:) = basis_ref;

elseif strcmp(basis_type,'P1') || strcmp(basis_type,'P1dc')
    result = zeros(1,3,num_pts);
    A = [1,-1,-1;0,1,0;0,0,1];
    if dx==0 && dy==0
        b = [one;x;y];
    elseif dx==1 && dy==0
        b = [zero;one;zero];
    elseif dx==0 && dy==1
        b = [zero;zero;one];
    else
        b = zeros(3,num_pts);
    end
    basis_ref = A*b;
    result(1,:,:) = basis_ref;

elseif strcmp(basis_type,'P2')
    result = zeros(1,6,num_pts);
    A = [1,-3,-3,2,4,2;
        0,-1,0,2,0,0;
        0,0,-1,0,0,2;
        0,4,0,-4,-4,0;
        0,0,0,0,4,0;
        0,0,4,0,-4,-4];
    if dx==0 && dy==0
        b = [one;x;y;x.^2;x.*y;y.^2];
    elseif dx==1 && dy==0
        b = [zero;one;zero;2*x;y;zero];
    elseif dx==0 && dy==1
        b = [zero;zero;one;zero;x;2*y];
    elseif dx==2 && dy==0
        b = [zero;zero;zero;2*one;zero;zero];
    elseif dx==1 && dy==1
        b = [zero;zero;zero;zero;one;zero];
    elseif dx==0 && dy==2
        b = [zero;zero;zero;zero;zero;2*one];
    else
        b = zeros(6,num_pts);
    end
    basis_ref = A*b;
    result(1,:,:) = basis_ref;

elseif strcmp(basis_type,'DG-P2-quad')
    result = zeros(1,6,num_pts);
    A = [-0.638559587411938,0.408594905200259,7.99630368687820,0.356305408700074,-7.99630368687822,-7.99630368687820;
        -0.638559587411938,7.99630368687821,0.408594905200258,-7.99630368687820,-7.99630368687822,0.356305408700075;
        0.126340726488395,-1.12120572260041,-1.12120572260041,0.356305408700077,8.70891450427837,0.356305408700075;
        0.138559587411936,-1.98973373528443,0.165973973290164,3.72483342138410,-0.165973973290164,-0.165973973290164;
        0.138559587411936,0.165973973290163,-1.98973373528444,-0.165973973290163,-0.165973973290163,3.72483342138410;
        1.87365927351161,-5.45993310748379,-5.45993310748379,3.72483342138411,7.61564081605839,3.72483342138411];
    if dx==0 && dy==0
        b = [one;x;y;x.^2;x.*y;y.^2];
    elseif dx==1 && dy==0
        b = [zero;one;zero;2*x;y;zero];
    elseif dx==0 && dy==1
        b = [zero;zero;one;zero;x;2*y];
    elseif dx==2 && dy==0
        b = [zero;zero;zero;2*one;zero;zero];
    elseif dx==1 && dy==1
        b = [zero;zero;zero;zero;one;zero];
    elseif dx==0 && dy==2
        b = [zero;zero;zero;zero;zero;2*one];
    else
        b = zeros(6,num_pts);
    end
    basis_ref = A*b;
    result(1,:,:) = basis_ref;

elseif strcmp(basis_type,'bubbleP1')
    result = zeros(1,4,num_pts);
    A = [1,-1,-1,0,-9,0,0,9,9,0;
        0,1,0,0,-9,0,0,9,9,0;
        0,0,1,0,-9,0,0,9,9,0;
        0,0,0,0,27,0,0,-27,-27,0];
    if dx==0 && dy==0
        b = [one;x;y;x.^2;x.*y;y.^2;x.^3;x.^2.*y;x.*y.^2;y.^3];
    elseif dx==1 && dy==0
        b = [zero;one;zero;2*x;y;zero;3*x.^2;2.*x.*y;y.^2;zero];
    elseif dx==0 && dy==1
        b = [zero;zero;one;zero;x;2*y;zero;x.^2;2*x.*y;3*y.^2];
    elseif dx==2 && dy==0
        b = [zero;zero;zero;2*one;zero;zero;6*x;2*y;zero;zero];
    elseif dx==1 && dy==1
        b = [zero;zero;zero;zero;one;zero;zero;2*x;2*y;zero];
    elseif dx==0 && dy==2
        b = [zero;zero;zero;zero;zero;2*one;zero;zero;2*x;6*y];
    end
    basis_ref = A*b;
    result(1,:,:) = basis_ref;

elseif strcmp(basis_type,'CR') || strcmp(basis_type,'DGP1')
    result = zeros(1,3,num_pts);
    A = [1,0,-2;-1,2,2;1,-2,0];
    if dx==0 && dy==0
        b = [one;x;y];
    elseif dx==1 && dy==0
        b = [zero;one;zero];
    elseif dx==0 && dy==1
        b = [zero;zero;one];
    else
        b = zeros(3,num_pts);
    end
    basis_ref = A*b;
    result(1,:,:) = basis_ref;
elseif strcmp(basis_type,'CR-P0')
    result = zeros(1,3,num_pts);
    A = [1/3,0,0;1/3,0,0;1/3,0,0];
    if dx==0 && dy==0
        b = [one;x;y];
    elseif dx==1 && dy==0
        b = [zero;one;zero];
    elseif dx==0 && dy==1
        b = [zero;zero;one];
    else
        b = zeros(3,num_pts);
    end
    basis_ref = A*b;
    result(1,:,:) = basis_ref;

elseif strcmp(basis_type,'CR-RT0')
    result = zeros(2,6,num_pts);
    A1 = [0,1,0;
        0,1,0;
        -1,1,0;
        0,1,0;
        0,1,0;
        -1,1,0];
    A2 = [-1,0,1;
        0,0,1;
        0,0,1;
        -1,0,1;
        0,0,1;
        0,0,1];
    if dx==0 && dy==0
        b = [one;x;y];
    elseif dx==1 && dy==0
        b = [zero;one;zero];
    elseif dx==0 && dy==1
        b = [zero;zero;one];
    else
        b = zeros(3,num_pts);
    end
    result(1,:,:) = A1*b;
    result(2,:,:) = A2*b;

elseif strcmp(basis_type,'RT0')
    result = zeros(2,3,num_pts);
    A1 = [0,1,0;
        0,1,0;
        -1,1,0];
    A2 = [-1,0,1;
        0,0,1;
        0,0,1];
    if dx==0 && dy==0
        b = [one;x;y];
    elseif dx==1 && dy==0
        b = [zero;one;zero];
    elseif dx==0 && dy==1
        b = [zero;zero;one];
    else
        b = zeros(3,num_pts);
    end
    result(1,:,:) = A1*b;
    result(2,:,:) = A2*b;

elseif strcmp(basis_type,'BR')
    result = zeros(2,9,num_pts);
    A1 = [0,0,0,0,0,0;
        0,1,0,0,0,0;
        0,0,2,0,-3,-3;
        -1,1,4,0,-3,-3;
        0,1,0,0,-3/2,0;
        0,0,0,0,-3/2,0;
        0,0,0,0,0,0;
        0,0,0,0,3,0;
        0,0,-6,0,6,6];
    A2 = [-1,4,1,-3,-3,0;
        0,2,0,-3,-3,0;
        0,0,1,0,0,0;
        0,0,0,0,0,0;
        0,0,0,0,-3/2,0;
        0,0,1,0,-3/2,0;
        0,-6,0,6,6,0;
        0,0,0,0,3,0;
        0,0,0,0,0,0];
    if dx==0 && dy==0
        b = [one;x;y;x.^2;x.*y;y.^2];
    elseif dx==1 && dy==0
        b = [zero;one;zero;2*x;y;zero];
    elseif dx==0 && dy==1
        b = [zero;zero;one;zero;x;2*y];
    elseif dx==2 && dy==0
        b = [zero;zero;zero;2*one;zero;zero];
    elseif dx==1 && dy==1
        b = [zero;zero;zero;zero;one;zero];
    elseif dx==0 && dy==2
        b = [zero;zero;zero;zero;zero;2*one];
    else
        b = zeros(6,num_pts);
    end

    result(1,:,:) = A1*b;
    result(2,:,:) = A2*b;
elseif strcmp(basis_type,'BR-RT0')
    result = zeros(2,9,num_pts);
    A1 = [0,0,0,0,0,0;
        0,0,0,0,0,0;
        0,0,0,0,0,0;
        0,0,0,0,0,0;
        0,0,0,0,0,0;
        0,0,0,0,0,0;
        0,1,0,0,0,0;
        0,1,0,0,0,0;
        -1,1,0,0,0,0];
    A2 = [0,0,0,0,0,0;
        0,0,0,0,0,0;
        0,0,0,0,0,0;
        0,0,0,0,0,0;
        0,0,0,0,0,0;
        0,0,0,0,0,0;
        -1,0,1,0,0,0;
        0,0,1,0,0,0;
        0,0,1,0,0,0];
    if dx==0 && dy==0
        b = [one;x;y;x.^2;x.*y;y.^2];
    elseif dx==1 && dy==0
        b = [zero;one;zero;2*x;y;zero];
    elseif dx==0 && dy==1
        b = [zero;zero;one;zero;x;2*y];
    elseif dx==2 && dy==0
        b = [zero;zero;zero;2*one;zero;zero];
    elseif dx==1 && dy==1
        b = [zero;zero;zero;zero;one;zero];
    elseif dx==0 && dy==2
        b = [zero;zero;zero;zero;zero;2*one];
    else
        b = zeros(6,num_pts);
    end

    result(1,:,:) = A1*b;
    result(2,:,:) = A2*b;
elseif strcmp(basis_type,'MTW')
    result = zeros(2,9,num_pts);
    A1 = [0,-14,0,30,132,0,-12,-120,-144,0;
        0,-2,0,12,-48,0,-12,24,72,0;
        0,25,0,-69,-114,0,42,144,90,0;
        0,-23,0,57,138,0,-30,-144,-126,0;
        -4,4,6,24,-24,0,-24,-24,36,0;
        2,16,-6,-66,-60,0,48,120,36,0;
        0,6,0,-6,-36,0,0,24,36,0;
        0,-3,0,9,18,0,-6,-24,-18,0;
        0,6,0,-18,-12,0,12,24,0,0;];
    A2 = [2,-6,16,0,-60,-66,0,36,120,48;
        -4,6,4,0,-24,24,0,36,-24,-24;
        0,0,-23,0,138,57,0,-126,-144,-30;
        0,0,25,0,-114,-69,0,90,144,42;
        0,0,-2,0,-48,12,0,72,24,-12;
        0,0,-14,0,132,30,0,-144,-120,-12;
        0,0,-6,0,12,18,0,0,-24,-12;
        0,0,3,0,-18,-9,0,18,24,6;
        0,0,-6,0,36,6,0,-36,-24,0;];
    if dx==0 && dy==0
        b = [one;x;y;x.^2;x.*y;y.^2;x.^3;x.^2.*y;x.*y.^2;y.^3];
    elseif dx==1 && dy==0
        b = [zero;one;zero;2*x;y;zero;3*x.^2;2*x.*y;y.^2;zero];
    elseif dx==0 && dy==1
        b = [zero;zero;one;zero;x;2*y;zero;x.^2;2*x.*y;3*y.^2];
    end
    result(1,:,:) = A1*b;
    result(2,:,:) = A2*b;

elseif strcmp(basis_type,'RT0-MTW')
    result = zeros(2,6,num_pts);
    A1 = [0	-8	0	21	42	0	-12	-48	-36	0
        0	1	0	-6	12	0	6	0	-18	0
        -1	10	0	-21	-42	0	12	48	36	0
        0	6	0	-6	-36	0	0	24	36	0
        0	-3	0	9	18	0	-6	-24	-18	0
        0	6	0	-18	-12	0	12	24	0	0];
    A2 = [-1	0	10	0	-42	-21	0	36	48	12
        0	0	1	0	12	-6	0	-18	0	6
        0	0	-8	0	42	21	0	-36	-48	-12
        0	0	-6	0	12	18	0	0	-24	-12
        0	0	3	0	-18	-9	0	18	24	6
        0	0	-6	0	36	6	0	-36	-24	0];
    if dx==0 && dy==0
        b = [one;x;y;x.^2;x.*y;y.^2;x.^3;x.^2.*y;x.*y.^2;y.^3];
    elseif dx==1 && dy==0
        b = [zero;one;zero;2*x;y;zero;3*x.^2;2*x.*y;y.^2;zero];
    elseif dx==0 && dy==1
        b = [zero;zero;one;zero;x;2*y;zero;x.^2;2*x.*y;3*y.^2];
    end
    result(1,:,:) = A1*b;
    result(2,:,:) = A2*b;

elseif strcmp(basis_type,'Ned1-1')
    A1 = [1,0,-1;0,0,-1;0,0,-1];
    A2 = [0,1,0;0,1,0;-1,1,0];
    if dx==0 && dy==0
        b = [one;x;y];
    elseif dx==1 && dy==0
        b = [zero;one;zero];
    elseif dx==0 && dy==1
        b = [zero;zero;one];
    end
    result(1,:,:) = A1*b;
    result(2,:,:) = A2*b;
elseif strcmp(basis_type,'Ned1-2')
    A1 = [-2,6,2,0,-8,0;
        4,-6,-12,0,8,8;
        0,0,4,0,0,-8;
        0,0,2,0,-8,0;
        0,0,-6,0,8,8;
        0,0,4,0,0,-8;
        0,0,16,0,-8,-16;
        0,0,-8,0,16,8];
    A2 = [0,-4,0,8,0,0;
        0,6,0,-8,-8,0;
        0,-2,0,0,8,0;
        0,-4,0,8,0,0;
        -4,12,6,-8,-8,0;
        2,-2,-6,0,8,0;
        0,-8,0,8,16,0;
        0,16,0,-16,-8,0];
    if dx==0 && dy==0
        b = [one;x;y;x.^2;x.*y;y.^2];
    elseif dx==1 && dy==0
        b = [zero;one;zero;2*x;y;zero];
    elseif dx==0 && dy==1
        b = [zero;zero;one;zero;x;2*y];
    end
    result(1,:,:) = A1*b;
    result(2,:,:) = A2*b;
elseif strcmp(basis_type,'Ned2-1')
    A1 = [-2,6,2;
        4,-6,-4;
        0,0,-4;
        0,0,2;
        0,0,2;
        0,0,-4;];
    A2 = [0,4,0;
        0,-2,0;
        0,-2,0;
        0,4,0;
        -4,4,6;
        2,-2,-6];
    if dx==0 && dy==0
        b = [one;x;y];
    elseif dx==1 && dy==0
        b = [zero;one;zero];
    elseif dx==0 && dy==1
        b = [zero;zero;one];
    end
    result(1,:,:) = A1*b;
    result(2,:,:) = A2*b;
else
    fprintf('Cannot handle this type of basis.');
end


