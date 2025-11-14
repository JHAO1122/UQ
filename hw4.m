function hw4()
    f1 = @(x) abs(sin(pi*x)).^3;
    f2 = @(x) abs(x);
    f3 = @(x) cos(pi*x);
    f4 = @(x) sign(x);
    max_index = 50;
    plot_polynomials(f1, 'f(x) = |sin(pi*x)|^3', max_index);
    plot_polynomials(f2, 'f(x) = |x|', max_index);
    plot_polynomials(f3, 'f(x) = cos(pi*x)', max_index);
    plot_polynomials(f4, 'f(x) = sign(x)', max_index);
    plot_error_vs_order(f1, 'f(x) = |sin(pi*x)|^3', max_index);
    plot_error_vs_order(f2, 'f(x) = |x|', max_index);
    plot_error_vs_order(f3, 'f(x) = cos(pi*x)', max_index);
    plot_error_vs_order(f4, 'f(x) = sign(x)', max_index);
end

function legendre = legendre_recursion(n, x)
    if n < 0, error("n must be non-negative"); end
    if n == 0
        legendre = ones(size(x));
    elseif n == 1
        legendre = x;
    else
        a_n0 = ones(size(x));
        a_n1 = x;
        for i = 2:n
            a_current = ((2*i - 1)*x.*a_n1 - (i - 1)*a_n0)/i;
            a_n0 = a_n1;
            a_n1 = a_current;
        end
        legendre = a_current;
    end
end

function hermite = hermite_probabilistic_recursion(n, x)

    if n < 0, error("n must be non-negative"); end
    if n == 0
        hermite = ones(size(x));
    elseif n == 1
        hermite = x;
    else
        h_n0 = ones(size(x));
        h_n1 = x;
        for i = 2:n
            h_current = x.*h_n1 - (i - 1)*h_n0;
            h_n0 = h_n1;
            h_n1 = h_current;
        end
        hermite = h_current;
    end
end

function [nodes, weights] = gauss_legendre_quadrature(n)
    i = 1:n-1;
    beta = i ./ sqrt(4*i.^2 - 1);
    T = diag(beta, 1) + diag(beta, -1);
    [V, D] = eig(T);
    nodes = diag(D)'; 
    weights = 2 * V(1, :).^2;
    [nodes, idx] = sort(nodes);
    weights = weights(idx);
end

function [nodes, weights] = gauss_hermite_quadrature(n)
    if n <= 0
        nodes = [];
        weights = [];
        return;
    end
    beta = sqrt(1:(n-1));
    T = diag(beta,1) + diag(beta,-1);
    [V, D] = eig(T);
    nodes = diag(D)';              
    weights = V(1,:).^2;          
    [nodes, idx] = sort(nodes);
    weights = weights(idx);
end
function polynomial_coefficients = polynomial_discrete(f, type, max_index, quad_points)
    polynomial_coefficients = zeros(1, max_index+1);

    for i = 0:max_index
        switch type
            case 'legendre'
                [nodes, weights] = gauss_legendre_quadrature(quad_points);
                f_vals = f(nodes);
                poly_vals = legendre_recursion(i, nodes);
                inner_product = sum(weights .* f_vals .* poly_vals);
                norm_sq = 2/(2*i + 1); 
                polynomial_coefficients(i+1) = inner_product / norm_sq;

            case 'hermite'
                [nodes, weights] = gauss_hermite_quadrature(quad_points);
                f_vals = f(nodes);
                poly_vals = hermite_probabilistic_recursion(i, nodes);
                inner_product = sum(weights .* f_vals .* poly_vals);
                norm_sq = factorial(i); 
                polynomial_coefficients(i+1) = inner_product / norm_sq;

            otherwise
                error("type must be 'legendre' or 'hermite'.");
        end
    end
end

function plot_polynomials(f, title_str, max_index)
    figure('Position', [100 100 1200 500]);
    quad_points = max(2*max_index+20, 200);
    subplot(1, 2, 1);
    coeffs_legendre = polynomial_discrete(f, 'legendre', max_index, quad_points);
    x_plot = linspace(-1, 1, 1000);
    y_original = f(x_plot);
    y_legendre = zeros(size(x_plot));
    for i = 0:max_index
        y_legendre = y_legendre + coeffs_legendre(i+1).*legendre_recursion(i, x_plot);
    end
    plot(x_plot, y_original, 'b-', 'LineWidth', 2, 'DisplayName', '原函数');
    hold on;
    plot(x_plot, y_legendre, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Legendre N=%d', max_index));
    xlabel('x'); ylabel('f(x)');
    title([title_str, ' - Legendre离散投影']);
    legend('show', 'Location', 'best');
    grid on;
    xlim([-1 1]);  
    subplot(1, 2, 2);
    coeffs_hermite = polynomial_discrete(f, 'hermite', max_index, quad_points);
    x_plot = linspace(-1, 1, 1000);  
    y_original = f(x_plot);
    y_hermite = zeros(size(x_plot));
    for i = 0:max_index
        y_hermite = y_hermite + coeffs_hermite(i+1).*hermite_probabilistic_recursion(i, x_plot);
    end
    plot(x_plot, y_original, 'b-', 'LineWidth', 2, 'DisplayName', '原函数');
    hold on;
    plot(x_plot, y_hermite, 'g--', 'LineWidth', 2, 'DisplayName', sprintf('Hermite N=%d', max_index));
    xlabel('x'); ylabel('f(x)');
    title([title_str, ' - Hermite离散投影']);
    legend('show', 'Location', 'best');
    grid on;
    xlim([-1 1]);  
    all_y = [y_original(:); y_legendre(:); y_hermite(:)];
    y_min = min(all_y);
    y_max = max(all_y);
    yrange = y_max - y_min;
    if yrange == 0
        y_min = y_min - 0.5;
        y_max = y_max + 0.5;
    else
        y_min = y_min - 0.05*yrange;
        y_max = y_max + 0.05*yrange;
    end
    subplot(1,2,1); ylim([y_min, y_max]);
    subplot(1,2,2); ylim([y_min, y_max]);
    sgtitle(title_str, 'FontSize', 14, 'FontWeight', 'bold');
end

function plot_error_vs_order(f, title_str, max_index)
    orders = 1:max_index;
    eL2_legendre = zeros(max_index,1);
    einf_legendre = zeros(max_index,1);
    eL2_hermite = zeros(max_index,1);
    einf_hermite = zeros(max_index,1);
    quad_points = max(2*max_index+20, 200);
    [nodes_large, ~] = gauss_hermite_quadrature(max(quad_points,200));
    x_limit = max(8, 1.2*max(abs(nodes_large)));
    x_plot_h = linspace(-x_limit, x_limit, 3000);
    w_h = exp(-x_plot_h.^2/2)/sqrt(2*pi); 
    for k = 1:max_index
        coeffs_legendre = polynomial_discrete(f, 'legendre', k, quad_points);
        xL = linspace(-1,1,2000);
        yL = f(xL);
        y_approx_L = zeros(size(xL));
        for i = 0:k
            y_approx_L = y_approx_L + coeffs_legendre(i+1).*legendre_recursion(i, xL);
        end
        errL = abs(y_approx_L - yL);
        eL2_legendre(k) = sqrt(trapz(xL, errL.^2));
        einf_legendre(k) = max(errL);
        coeffs_hermite = polynomial_discrete(f, 'hermite', k, quad_points);
        yH = f(x_plot_h);
        y_approx_H = zeros(size(x_plot_h));
        for i = 0:k
            y_approx_H = y_approx_H + coeffs_hermite(i+1).*hermite_probabilistic_recursion(i, x_plot_h);
        end
        errH = abs(y_approx_H - yH);
        eL2_hermite(k) = sqrt(trapz(x_plot_h, errH.^2 .* w_h));
        einf_hermite(k) = max(errH .* sqrt(w_h)); 
    end
    figure('Position',[100 100 1400 800]);
    subplot(2,2,1);
    semilogy(orders, eL2_legendre, 'ro-', 'LineWidth',2, 'DisplayName','Legendre');
    hold on;
    semilogy(orders, eL2_hermite, 'bs-', 'LineWidth',2, 'DisplayName','Hermite');
    xlabel('多项式阶数 N'); ylabel('L^2 误差');
    title('(a) L^2 误差 (semilogy)');
    legend('show'); grid on;

    subplot(2,2,2);
    semilogy(orders, einf_legendre, 'ro-', 'LineWidth',2, 'DisplayName','Legendre');
    hold on;
    semilogy(orders, einf_hermite, 'bs-', 'LineWidth',2, 'DisplayName','Hermite');
    xlabel('多项式阶数 N'); ylabel('L^\infty 误差');
    title('(b) L^\infty 误差 (semilogy)');
    legend('show'); grid on;

    subplot(2,2,3);
    loglog(orders, eL2_legendre, 'ro-', 'LineWidth',2, 'DisplayName','Legendre');
    hold on;
    loglog(orders, eL2_hermite, 'bs-', 'LineWidth',2, 'DisplayName','Hermite');
    xlabel('多项式阶数 N'); ylabel('L^2 误差');
    title('(c) L^2 误差 (loglog)');
    legend('show'); grid on;

    subplot(2,2,4);
    loglog(orders, einf_legendre, 'ro-', 'LineWidth',2, 'DisplayName','Legendre');
    hold on;
    loglog(orders, einf_hermite, 'bs-', 'LineWidth',2, 'DisplayName','Hermite');
    xlabel('多项式阶数 N'); ylabel('L^\infty 误差');
    title('(d) L^\infty 误差 (loglog)');
    legend('show'); grid on;

    sgtitle([title_str, ' - 离散投影收敛性分析'], 'FontSize',14,'FontWeight','bold');
end
