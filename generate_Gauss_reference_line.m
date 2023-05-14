function [weights,pts] = generate_Gauss_reference_line(Gauss_type)
% compute Gaussian pts and weights on reference edges.

if Gauss_type == 1
    pts = 0;
    weights = 2;
elseif Gauss_type == 2
    pts = [-sqrt(1/3),sqrt(1/3)];
    weights = [1,1];
elseif Gauss_type == 3
    pts = [-sqrt(3/5),0,sqrt(3/5)];
    weights = [5/9,8/9,5/9];
elseif Gauss_type == 4
    pts = [-sqrt(525-70*sqrt(30))/35,sqrt(525-70*sqrt(30))/35,...
        -sqrt(525+70*sqrt(30))/35,sqrt(525+70*sqrt(30))/35];
    weights = [(18+sqrt(30))/36,(18+sqrt(30))/36,(18-sqrt(30))/36,(18-sqrt(30))/36];
elseif Gauss_type == 5
    pts = [0,-sqrt(245-14*sqrt(70))/21,sqrt(245-14*sqrt(70))/21,...
        -sqrt(245+14*sqrt(70))/21,sqrt(245+14*sqrt(70))/21];
    weights = [128/225,(322+13*sqrt(70))/900,(322+13*sqrt(70))/900,...
        (322-13*sqrt(70))/900,(322-13*sqrt(70))/900];
else
    fprintf('cannot handle this type of Gaussian quadrature!');
end