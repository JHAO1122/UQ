function hw2()
    rng(0);
    n = 400;
    M_list = [100, 500, 1000, 2000]; 
    results_method1 = [];
    results_method2 = [];
    results_method3 = [];
    for i = 1:length(M_list)
        M_c = M_list(i);
        [Y_all1, Y_true_all1] = mcode(M_c, n, 1);
        mean_numerical1 = mean(Y_all1, 2);  
        mean_true1 = Y_true_all1(:,1);  
        errors_mean1 = mean_numerical1 - mean_true1;  
        eL2_1 = sqrt(mean(errors_mean1.^2));  
        einf_1 = max(abs(errors_mean1));       
        results_method1(i).M = M_c;
        results_method1(i).eL2 = eL2_1;
        results_method1(i).einf = einf_1;
        [Y_all2, Y_true_all2] = mcode(M_c, n, 2);
        mean_numerical2 = mean(Y_all2, 2);
        mean_true2 = Y_true_all2(:,1);
        errors_mean2 = mean_numerical2 - mean_true2;
        eL2_2 = sqrt(mean(errors_mean2.^2));
        einf_2 = max(abs(errors_mean2));        
        results_method2(i).M = M_c;
        results_method2(i).eL2 = eL2_2;
        results_method2(i).einf = einf_2;

        [u_pde_all, u_true_all] = mcpde(M_c, n, 5, 1, 0.9);
        mean_numerical_pde = mean(u_pde_all, 3);
        errors_mean_pde = mean_numerical_pde - u_true_all;
        eL2_pde = sqrt(mean(errors_mean_pde(:).^2));
        einf_pde = max(abs(errors_mean_pde(:))); 
        results_method3(i).M = M_c;
        results_method3(i).eL2 = eL2_pde;
        results_method3(i).einf = einf_pde;     
    end
    
    create_error_table(results_method1, 'EulerMC');
    create_error_table(results_method2, 'RKMC');
    create_pde_error_table(results_method3, "Scheme3 PDE");
end


function Y =  ode(n, xa, xb, c, f, type)
Y = zeros(1, n+1) ;
X = zeros(1, n+1) ;
X(1) = xa;
Y(1) = c;
delta = (xb - xa) / n;

if type == 1
for i = 1:n
    Y(i+1) = f(X(i), Y(i)) * delta + Y(i) ;
    X(i+1) = X(i) + delta ;
end
end
if type == 2
for i = 1:n
    k1 = Y(i) + delta * f(X(i), Y(i));
    n1 = X(i) + delta;
    n_half = X(i) + delta/2;
    k2 = 3/4 * Y(i) + 1/4 * k1 + delta/4 * f(n1, k1);
    Y(i+1) = 1/3 * Y(i) + 2/3 * k2 + 2 * delta/3 * f(n_half, k2);
    X(i+1) = X(i) + delta;
end
end
end

function [Y_all, Y_true_all] = mcode(M, n, type2)
    Y_all = zeros(n+1, M);
    t = linspace(0, 2, n+1);
    true_expect = exp(t.^2/2)';
    Y_true_all = repmat(true_expect, 1, M);
    if type2 == 1
        for j = 1:M
        alpha = normrnd(0, 1);
        beta = 1;
        Y = ode(n, 0, 2, beta, @(x,y) (-1)*alpha*y, 1);
        Y_all(:, j) = Y';
        end
    end

    if type2 == 2
        for j = 1:M
        alpha = normrnd(0, 1);
        beta = 1;
        Y = ode(n, 0, 2, beta, @(x,y) (-1)*alpha*y, 2);
        Y_all(:, j) = Y';
        end
    end
end


function create_error_table(results, method)
        fprintf('\n%s误差表:\n', method);
        fprintf('M\t\teL2\t\torder\t\te∞\t\torder\n');
        fprintf('------------------------------------------------\n');
        for i = 1:length(results)
            M = results(i).M;
            eL2 = results(i).eL2;
            einf = results(i).einf;
            if i == 1
                order_L2 = '-';
                order_inf = '-';
            else
                pre_eL2 = results(i-1).eL2;
                pre_einf = results(i-1).einf;
                pre_M = results(i-1).M; 
                order_L2 = log(pre_eL2 / eL2) / log(M / pre_M);
                order_inf = log(pre_einf / einf) / log(M / pre_M);
            end
            if i == 1
                fprintf('%d\t\t%.2e\t%s\t\t%.2e\t%s\n', ...
                       M, eL2, order_L2, einf, order_inf);
            else
                fprintf('%d\t\t%.2e\t%.3f\t\t%.2e\t%.3f\n', ...
                        M, eL2, order_L2, einf, order_inf);
            end
        end
    end





function [u_pde_all, u_pde_true_all] = mcpde(M, N, T, k, txratio)
    delta_x = 2 * pi / N;
    delta_t = txratio * delta_x;
    nt = floor(T / delta_t);
    u_pde_all = zeros(N+1, nt+1, M);
    u_pde_true_all = zeros(N+1, nt+1);
    x = linspace(0, 2*pi, N+1);
    t = 0:delta_t:T;
    for i = 1:nt+1
        for j = 1:N+1
            if abs(t(i)) < 1e-10
                u_pde_true_all(j, i) = sin(k*x(j));
            else
                u_pde_true_all(j, i) = (1/(0.2*t(i))) * (cos(k*(x(j)-1.1*t(i))) - cos(k*(x(j)-0.9*t(i))));
            end
        end
    end
    for j = 1:M
        a_j = 1 + 0.1*(2*rand-1);
        u_pde = pde(N, T, k, a_j, txratio);
        u_pde_all(:, :, j) = u_pde;
    end
end


function u = pde(N, T, k, a, txratio)
    delta_x = 2 * pi / N;
    delta_t = txratio * delta_x;
    nt = floor(T / delta_t);
    u = zeros(N+1, nt+1);
    x = linspace(0, 2*pi, N+1);
    for j = 1:N+1
            u(j, 1) = sin(k * x(j));
    end
    for i = 1:nt
        for j = 2:N
        u(j, i+1) = u(j, i) + a * txratio * (u(j+1, i) - u(j, i));
        end
        u(1, i+1) = u(1, i) + a * txratio * (u(2, i) - u(1, i));
        u(N+1, i+1) = u(1, i+1);
    end
end



function create_pde_error_table(results, method_name)
    fprintf('\n%s误差表:\n', method_name);
    fprintf('M\t\teL2\t\torder\t\te∞\t\torder\n');
    fprintf('------------------------------------------------\n');
    table_results = [];
    for i = 1:length(results)
        table_results(i).M = results(i).M;
        table_results(i).eL2 = results(i).eL2;
        table_results(i).einf = results(i).einf;
    end
    for i = 1:length(table_results)
        M = table_results(i).M;
        eL2 = table_results(i).eL2;
        einf = table_results(i).einf; 
        if i == 1
            order_L2 = '-';
            order_inf = '-';
        else
            pre_eL2 = table_results(i-1).eL2;
            pre_einf = table_results(i-1).einf;
            pre_M = table_results(i-1).M; 
            order_L2 = log(pre_eL2 / eL2) / log(M / pre_M);
            order_inf = log(pre_einf / einf) / log(M / pre_M);
        end
        if i == 1
            fprintf('%d\t\t%.2e\t%s\t\t%.2e\t%s\n', ...
                   M, eL2, order_L2, einf, order_inf);
        else
            fprintf('%d\t\t%.2e\t%.3f\t\t%.2e\t%.3f\n', ...
                    M, eL2, order_L2, einf, order_inf);
        end
    end
end
