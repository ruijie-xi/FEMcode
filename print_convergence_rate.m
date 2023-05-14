function [rate_list] = print_convergence_rate(h_list,err_list)

N_test = size(h_list,2);

assert(N_test>1,'number of tests should be more than 1!');
rate_list = zeros(1,N_test - 1);

for i = 2:N_test
    rate_list(i-1) = log(err_list(i)/err_list(i-1))/log(h_list(i)/h_list(i-1));
    fprintf('h1 = %.2e, h2 = %.2e, rate = %.2f\n',h_list(i-1),h_list(i),rate_list(i-1));
end

