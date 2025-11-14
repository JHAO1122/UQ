function hw3()
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
    if n < 0
        printf("Error: n must be a non-negative integer.\n");
    end
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

function Hermite = hermite_recursion(n, x)
    if n < 0
        printf("Error: n must be a non-negative integer.\n");
    end
    if n == 0
        Hermite = ones(size(x));
    elseif n == 1
        Hermite = x;
    else
        a_n0 = ones(size(x));
        a_n1 = x;
        for i = 2:n
            a_current = x.*a_n1 - (i - 1)*a_n0;
            a_n0 = a_n1;
            a_n1 = a_current;
        end
        Hermite = a_current;
    end
end

function polynomial_coefficients = polynomial(f, type, max_index)
    polynomial_coefficients = zeros(1, max_index+1);
    warning off
        for i = 0:max_index
            switch type
                case 'legendre'
                    integrand = @(x) f(x).*legendre_recursion(i, x);
                    a = -1; b = 1;
                    norm = 2/(2*i + 1);
                case 'hermite'
                    integrand = @(x) f(x).*hermite_recursion(i, x) .* exp(-x.^2/2) / sqrt(2*pi);
                    a = -inf; b = inf;
                    norm = factorial(i);
                otherwise
                    printf("Error: type must be 'legendre' or 'hermite'.\n");
            end
            if isinf(a) || isinf(b)
                integral = quadgk(integrand, -8, 8);
            else
                integral = quadgk(integrand, a, b);
            end
            polynomial_coefficients(i+1) = integral / norm;
        end
end

function plot_polynomials(f, title_str, max_index)
    figure('position', [100 100 1200 500]);
    subplot(1, 2, 1);
    coeffs_legendre = polynomial(f, 'legendre', max_index);
    x_step = linspace(-1, 1, 100);
    y_original = f(x_step);
    y_legendre = zeros(size(x_step));
    for i = 0:max_index
        y_legendre = y_legendre + coeffs_legendre(i+1).*legendre_recursion(i, x_step);
    end
    plot(x_step, y_original, 'b--', 'LineWidth', 2, 'DisplayName', '原函数')
    hold on
    plot(x_step, y_legendre, 'r', 'LineWidth', 2, 'DisplayName', sprintf('Legendre多项式, 最大阶数%d', max_index))
    xlabel('x')
    ylabel('y')
    title(title_str, '-legendre多项式逼近')
    legend('show')
    grid on;
    subplot(1, 2, 2);
    coeffs_hermite = polynomial(f, 'hermite', max_index);
    x_step2 = linspace(-1, 1, 100);
    y_hermite = zeros(size(x_step2));
    for i = 0:max_index
        y_hermite = y_hermite + coeffs_hermite(i+1).*hermite_recursion(i, x_step2);
    end
    plot(x_step2, f(x_step2), 'b--', 'LineWidth', 2, 'DisplayName', '原函数')
    hold on
    plot(x_step2, y_hermite, 'r', 'LineWidth', 2, 'DisplayName', sprintf('Hermite多项式, 最大阶数%d', max_index))
    xlabel('x')
    ylabel('y')
    title(title_str, '-Hermite多项式逼近')
    legend('show')
    grid on;
end

% function results = compute_errors(f, type, orders)
%     results = [ ];
%     for k = 1:length(orders)
%         max_index = orders(k);
%         polynomial_coefficients = polynomial(f, type, max_index);
%         x_step = linspace(-1, 1, 100);
%         y_original = f(x_step);
%         y_approx = zeros(size(x_step));
%         if strcmp(type, 'legendre')
%             weights = ones(size(x_step));
%         else
%             weights = exp(-x_step.^2);
%         end
%         for i = 0:max_index
%             if strcmp(type, 'legendre')
%                 func = legendre_recursion(i, x_step);
%             else
%                 func = hermite_recursion(i, x_step);
%             end
%             y_approx = y_approx + polynomial_coefficients(i+1).*func;
%         end
%         error = abs(y_approx - y_original);
%         eL2 = sqrt(trapz(x_step, weights .* error.^2));
%         einf = max(weights .* error);
%         results(k).N = max_index;
%         results(k).eL2 = eL2;
%         results(k).einf = einf;
%     end
% end

function plot_error_vs_order(f, title_str, max_index)
    orders = 1:max_index;
    eL2_legendre = zeros(max_index, 1);
    einf_legendre = zeros(max_index, 1);
    eL2_hermite = zeros(max_index, 1);
    einf_hermite = zeros(max_index, 1);
    for k = 1:max_index
        current_order = k;  
        coeffs_legendre = polynomial(f, 'legendre', current_order);
        x_step = linspace(-1, 1, 100);
        y_original = f(x_step);
        y_approx_legendre = zeros(size(x_step));
        weights_legendre = ones(size(x_step)); 
        for i = 0:current_order
            y_approx_legendre = y_approx_legendre + coeffs_legendre(i+1).*legendre_recursion(i, x_step);
        end       
        error_legendre = abs(y_approx_legendre - y_original);
        eL2_legendre(k) = sqrt(trapz(x_step, weights_legendre .* error_legendre.^2));
        einf_legendre(k) = max(weights_legendre .* error_legendre);
        coeffs_hermite = polynomial(f, 'hermite', current_order);
        y_approx_hermite = zeros(size(x_step));
        weights_hermite = exp(-x_step.^2); 
        for i = 0:current_order
            y_approx_hermite = y_approx_hermite + coeffs_hermite(i+1).*hermite_recursion(i, x_step);
        end
        error_hermite = abs(y_approx_hermite - y_original);
        eL2_hermite(k) = sqrt(trapz(x_step, weights_hermite .* error_hermite.^2));
        einf_hermite(k) = max(weights_hermite .* error_hermite);
    end
    
    figure('position', [100 100 1400 800]);
    subplot(2, 2, 1);
    loglog(orders, eL2_legendre, 'ro-', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Legendre');
    hold on;
    loglog(orders, eL2_hermite, 'bs-', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Hermite');
    xlabel('多项式阶数 N');
    ylabel('L2误差');
    title('(a) loglog坐标 - L2误差');
    legend('show', 'Location', 'best');
    grid on;
    subplot(2, 2, 2);
    loglog(orders, einf_legendre, 'ro-', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Legendre');
    hold on;
    loglog(orders, einf_hermite, 'bs-', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Hermite');
    xlabel('多项式阶数 N');
    ylabel('L∞误差');
    title('(b) loglog坐标 - L∞误差');
    legend('show', 'Location', 'best');
    grid on;
    subplot(2, 2, 3);
    semilogy(orders, eL2_legendre, 'ro-', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Legendre');
    hold on;
    semilogy(orders, eL2_hermite, 'bs-', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Hermite');
    xlabel('多项式阶数 N');
    ylabel('L2误差');
    title('(c) semilogy坐标 - L2误差');
    legend('show', 'Location', 'best');
    grid on;
    subplot(2, 2, 4);
    semilogy(orders, einf_legendre, 'ro-', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Legendre');
    hold on;
    semilogy(orders, einf_hermite, 'bs-', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Hermite');
    xlabel('多项式阶数 N');
    ylabel('L∞误差');
    title('(d) semilogy坐标 - L∞误差');
    legend('show', 'Location', 'best');
    grid on;
    sgtitle([title_str, ' - 多项式逼近误差分析'], 'FontSize', 14, 'FontWeight', 'bold');
end



