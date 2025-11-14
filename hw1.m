function hw1()
    [~] = pde(20, 5, 1, 1, 0.9, 1);
    fq = @(x, y) -5 * y;
    N_values = [20, 40, 80, 160];
    results_euler = [];
    results_midform= [];
    results_rk3 = [];
    for i = 1:length(N_values)
        N = N_values(i);
        Y1 = ode(N, 0, 2, 1, fq, 1);
        Y_true = ode_true(N, 0, 2);
        Y2 = ode(N, 0, 2, 1, fq, 2);
        Y3 = ode(N, 0, 2, 1, fq, 4);
        eL2_1 = sqrt(sum((Y1 - Y_true).^2) / (N+1));
        eL2_2 = sqrt(sum((Y2 - Y_true).^2) / (N+1));
        eL2_3 = sqrt(sum((Y3 - Y_true).^2) / (N+1));
        einf_1 = max(abs(Y1 - Y_true));
        einf_2 = max(abs(Y2 - Y_true));
        einf_3 = max(abs(Y3 - Y_true));
        results_euler(i).N = N;
        results_euler(i).eL2 = eL2_1;
        results_euler(i).einf = einf_1;
        results_midform(i).N = N;
        results_midform(i).eL2 = eL2_2;
        results_midform(i).einf = einf_2;
        results_rk3(i).N = N;
        results_rk3(i).eL2 = eL2_3;
        results_rk3(i).einf = einf_3;
    end

N_a_values = [20, 40, 80];
k_a_values = [1, 5, 10];
T_a_values = [5, 10, 20];
results_pde_a = [];
for i = 1:length(N_a_values)
    N_a = N_a_values(i);
    for j = 1:length(k_a_values)
        for l = 1:length(T_a_values)
            Y = pde(N_a, T_a_values(l), k_a_values(j), 1, 0.9, 1);
            Y_true = pde_true(N_a, T_a_values(l), k_a_values(j), 1, 0.9); 
            error = Y(:, end) - Y_true(:, end);
            eL2 = sqrt(sum(error.^2) / (N_a+1));
            einf = max(abs(error));
            results_pde_a(i,j,l).N = N_a;
            results_pde_a(i,j,l).k = k_a_values(j);  
            results_pde_a(i,j,l).T = T_a_values(l);  
            results_pde_a(i,j,l).eL2 = eL2;
            results_pde_a(i,j,l).einf = einf;
        end
    end
end

N_b_values_sch2 = [20, 40, 80, 160];
k_b_values_sch2 = [1, 5, 10];
T_b_values_sch2 = [5, 10, 20];
results_pde_b_sch2 = [];  
for i = 1:length(N_b_values_sch2)
    N_b_sch2 = N_b_values_sch2(i);
    for j = 1:length(k_b_values_sch2)
        for l = 1:length(T_b_values_sch2)
            Y = pde(N_b_sch2, T_b_values_sch2(l), k_b_values_sch2(j), 1, 0.9, 2);
            Y_true = pde_true(N_b_sch2, T_b_values_sch2(l), k_b_values_sch2(j), 1, 0.9); 
            error = Y(:, end) - Y_true(:, end);
            eL2 = sqrt(sum(error.^2) / (N_b_sch2+1));
            einf = max(abs(error));
            results_pde_b_sch2(i,j,l).N = N_b_sch2;
            results_pde_b_sch2(i,j,l).k = k_b_values_sch2(j);  
            results_pde_b_sch2(i,j,l).T = T_b_values_sch2(l);  
            results_pde_b_sch2(i,j,l).eL2 = eL2;
            results_pde_b_sch2(i,j,l).einf = einf;
        end
    end
end
    
N_b_values_sch3 = [20, 40, 80, 160];
k_b_values_sch3 = [1, 5, 10];
T_b_values_sch3 = [5, 10, 20];
results_pde_b_sch3 = [];  
for i = 1:length(N_b_values_sch3)
    N_b_sch3 = N_b_values_sch3(i);
    for j = 1:length(k_b_values_sch3)
        for l = 1:length(T_b_values_sch3)
            Y = pde(N_b_sch3, T_b_values_sch3(l), k_b_values_sch3(j), 1, 0.9, 3);
            Y_true = pde_true(N_b_sch3, T_b_values_sch3(l), k_b_values_sch3(j), 1, 0.9); 
            error = Y(:, end) - Y_true(:, end);
            eL2 = sqrt(sum(error.^2) / (N_b_sch3+1));
            einf = max(abs(error));
            results_pde_b_sch3(i,j,l).N = N_b_sch3;
            results_pde_b_sch3(i,j,l).k = k_b_values_sch3(j);  
            results_pde_b_sch3(i,j,l).T = T_b_values_sch3(l);  
            results_pde_b_sch3(i,j,l).eL2 = eL2;
            results_pde_b_sch3(i,j,l).einf = einf;
        end
    end
end

N_c_values = [10, 20, 40];
k_c_values = [1, 5, 10];
T_c_values = [5, 10, 20];
results_pde_c = [];  
for i = 1:length(N_c_values)
    N_c = N_c_values(i);
    for j = 1:length(k_c_values)
        for l = 1:length(T_c_values)
            Y = pde(N_c, T_c_values(l), k_c_values(j), -1, 0.9, 3);
            Y_true = pde_true(N_c, T_c_values(l), k_c_values(j), -1, 0.9); 
            error = Y(:, end) - Y_true(:, end);
            eL2 = sqrt(sum(error.^2) / (N_c+1));
            einf = max(abs(error));
            results_pde_c(i,j,l).N = N_c;
            results_pde_c(i,j,l).k = k_c_values(j);  
            results_pde_c(i,j,l).T = T_c_values(l);  
            results_pde_c(i,j,l).eL2 = eL2;
            results_pde_c(i,j,l).einf = einf;            
        end 
    end
end
    create_error_table(results_euler, 'Euler');
    create_error_table(results_midform, 'Mid-form');
    create_error_table(results_rk3, 'Runge-Kuta3');
    create_pde_error_table(results_pde_a, 1, 1, 'scheme1 (k=1, T=5)');
    create_pde_error_table(results_pde_a, 2, 2, 'scheme1 (k=5, T=10)');
    create_pde_error_table(results_pde_a, 3, 3, 'scheme1 (k=10, T=20)');
    create_pde_error_table(results_pde_b_sch2, 1, 1, 'scheme2 (k=1, T=5)');
    create_pde_error_table(results_pde_b_sch2, 2, 2, 'scheme2 (k=5, T=10)');
    create_pde_error_table(results_pde_b_sch2, 3, 3, 'scheme2 (k=10, T=20)');
    create_pde_error_table(results_pde_b_sch3, 1, 1, 'scheme3 (k=1, T=5)');
    create_pde_error_table(results_pde_b_sch3, 2, 2, 'scheme3 (k=5, T=10)');
    create_pde_error_table(results_pde_b_sch3, 3, 3, 'scheme3 (k=10, T=20)');
    create_pde_error_table(results_pde_c, 1, 1, 'scheme3 (k=1, T=5, a=-1)');
    create_pde_error_table(results_pde_c, 2, 2, 'scheme3 (k=5, T=10, a=-1)');
    create_pde_error_table(results_pde_c, 3, 3, 'scheme3 (k=10, T=20, a=-1)');
end

function Y = ode(n, xa, xb, c, f, type)
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
X(2) = X(1) + delta;
Y(2) = Y(1) + delta * f(X(1), Y(1));
for i = 2:n
    Y(i+1) = 2 * delta * f(X(i), Y(i)) + Y(i-1); 
    X(i+1) = X(i) + delta;
end
end

if type ==3
for i = 1:n
    k1 = f(X(i), Y(i));
    k2 = f(X(i) + delta/2, Y(i) + delta * k1/2);
    k3 = f(X(i) + delta, Y(i) - delta * k1 + 2 * delta * k2);
    Y(i+1) = Y(i) + delta/6 * (k1 + 4 * k2 + k3);
    X(i+1) = X(i) + delta;
end
end

if type == 4
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
function Y_true = ode_true(n, xa, xb)
Y_true = zeros(1, n+1);
X = zeros(1, n+1);
for i = 1:n+1
    X(i) = (i-1) * (xb - xa) / n + xa;
    Y_true(i) = exp(-5 * X(i));
end
end

function create_error_table(results, method)
        fprintf('\n%s误差表:\n', method);
        fprintf('N\t\teL2\t\torder\t\te∞\t\torder\n');
        fprintf('------------------------------------------------\n');
        for i = 1:length(results)
            N = results(i).N;
            eL2 = results(i).eL2;
            einf = results(i).einf;
            if i == 1
                order_L2 = '-';
                order_inf = '-';
            else
                pre_eL2 = results(i-1).eL2;
                pre_einf = results(i-1).einf;
                pre_N = results(i-1).N; 
                order_L2 = log(pre_eL2 / eL2) / log(N / pre_N);
                order_inf = log(pre_einf / einf) / log(N / pre_N);
            end
            if i == 1
                fprintf('%d\t\t%.2e\t%s\t\t%.2e\t%s\n', ...
                       N, eL2, order_L2, einf, order_inf);
            else
                fprintf('%d\t\t%.2e\t%.3f\t\t%.2e\t%.3f\n', ...
                        N, eL2, order_L2, einf, order_inf);
            end
        end
end

function create_pde_error_table(results, k_idx, T_idx, method_name)
    fprintf('\n%s误差表:\n', method_name);
    fprintf('N\t\teL2\t\torder\t\te∞\t\torder\n');
    fprintf('------------------------------------------------\n');
    table_results = [];
    for i = 1:size(results, 1)
        table_results(i).N = results(i, k_idx, T_idx).N;
        table_results(i).eL2 = results(i, k_idx, T_idx).eL2;
        table_results(i).einf = results(i, k_idx, T_idx).einf;
    end
    for i = 1:length(table_results)
        N = table_results(i).N;
        eL2 = table_results(i).eL2;
        einf = table_results(i).einf; 
        if i == 1
            order_L2 = '-';
            order_inf = '-';
        else
            pre_eL2 = table_results(i-1).eL2;
            pre_einf = table_results(i-1).einf;
            pre_N = table_results(i-1).N; 
            order_L2 = log(pre_eL2 / eL2) / log(N / pre_N);
            order_inf = log(pre_einf / einf) / log(N / pre_N);
        end
        if i == 1
            fprintf('%d\t\t%.2e\t%s\t\t%.2e\t%s\n', ...
                   N, eL2, order_L2, einf, order_inf);
        else
            fprintf('%d\t\t%.2e\t%.3f\t\t%.2e\t%.3f\n', ...
                    N, eL2, order_L2, einf, order_inf);
        end
    end
end




function u = pde(N, T, k, a, txratio, type)
    delta_x = 2 * pi / N;
    delta_t = txratio * delta_x;
    nt = floor(T / delta_t);
    u = zeros(N+1, nt+1);
    x = linspace(0, 2*pi, N+1);
    for j = 1:N+1
            u(j, 1) = sin(k * x(j));
    end
    if type == 1
    for i = 1: nt
        for j = 2: N 
            u(j, i+1) = u(j, i) + a / 2 * txratio * (u(j+1, i) - u(j-1, i));
        end
        u(1, i+1) = u(1, i) + a/2 * txratio * (u(2, i) - u(N, i));
        u(N+1, i+1) = u(1, i+1);
    end
    end

    if type == 2
    for i = 1: nt
        for j = 2: N
            u(j, i+1) = (u(j + 1, i) + u(j - 1, i)) / 2 + a / 2 * txratio * (u(j+1, i) - u(j-1, i));
        end
        u(1, i+1) = (u(2, i) + u(N, i))/2 + a/2 * txratio * (u(2, i) - u(N, i));
        u(N+1, i+1) = u(1, i+1);
    end
    end

    
    if type == 3 && a > 0
        for i = 1:nt
            for j = 2:N
            u(j, i+1) = u(j, i) + a * txratio * (u(j+1, i) - u(j, i));
            end
            u(1, i+1) = u(1, i) + a * txratio * (u(2, i) - u(1, i));
            u(N+1, i+1) = u(1, i+1);
        end
    end

    if type == 3 && a < 0
        for i = 1:nt
            for j = 2:N+1
            u(j, i+1) = u(j, i) + a * txratio * (u(j, i) - u(j-1, i));
            end
        u(1, i+1) = u(1, i) + a * txratio * (u(1, i) - u(N, i));
        u(N+1, i+1) = u(1, i+1);
        end
    end    
end

function u = pde_true(N, T, k, a, txratio)
    delta_x = 2 * pi / N;
    delta_t = txratio * delta_x;
    nt = floor(T / delta_t);
    u = zeros(N+1, nt+1);
    x = linspace(0, 2*pi, N+1);
    t = 0:delta_t:T;
    for i = 1:nt+1
        for j = 1:N+1
            u(j, i) = sin(k * (x(j) + a*t(i)));
        end
    end
end
