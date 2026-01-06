% =========================================================================
% Thapar Institute of Engineering and Technology, Patiala
% School of Mathematics
% Optimization Methods (UMA-035)
% Lab Experiment - 9: Fibonacci Search Technique
% =========================================================================

function Lab9_Fibonacci_Search()
    clc;
    clear;
    
    fprintf('========================================\n');
    fprintf('Fibonacci Search Technique\n');
    fprintf('========================================\n\n');
    
    fprintf('This program implements Fibonacci search for 1D optimization\n\n');
    
    % Example problem from lab document
    fprintf('Example Problem:\n');
    fprintf('Minimize f(x) = x(x-2) on [0, 1.5]\n');
    fprintf('Target: Interval of uncertainty = 0.25*L0\n\n');
    
    % Define the objective function
    f = @(x) x.*(x - 2);
    
    % Define interval
    a = 0;
    b = 1.5;
    
    % Target measure of effectiveness
    target_measure = 0.25;
    
    % Solve using Fibonacci search
    [x_opt, f_opt, iterations, final_interval] = fibonacci_search(f, a, b, ...
                                                 target_measure, true, 'min');
    
    fprintf('\n========================================\n');
    fprintf('FINAL RESULTS\n');
    fprintf('========================================\n');
    fprintf('Optimal point: x* = %.6f\n', x_opt);
    fprintf('Optimal value: f(x*) = %.6f\n', f_opt);
    fprintf('Final interval: [%.6f, %.6f]\n', final_interval(1), final_interval(2));
    fprintf('Final interval length: %.6f\n', final_interval(2) - final_interval(1));
    fprintf('Iterations: %d\n', iterations);
    
    % Plot the function and search process
    plot_fibonacci_search(f, a, b, x_opt);
end

function [x_opt, f_opt, n, final_interval] = fibonacci_search(f, a, b, ...
                                             target_measure, verbose, opt_type)
    % Fibonacci Search Algorithm
    % Inputs:
    %   f: objective function handle
    %   a, b: initial interval [a, b]
    %   target_measure: desired final interval / initial interval
    %   verbose: display iterations
    %   opt_type: 'min' or 'max'
    
    % Step 1: Calculate measure of effectiveness
    L0 = b - a;
    measure = target_measure;
    
    if verbose
        fprintf('Initial Setup:\n');
        fprintf('--------------\n');
        fprintf('Initial interval: [%.4f, %.4f]\n', a, b);
        fprintf('Initial length L0 = %.4f\n', L0);
        fprintf('Target measure = %.4f\n', measure);
        fprintf('Target final interval = %.4f\n\n', measure * L0);
    end
    
    % Step 2: Find smallest n such that 1/F_n <= measure
    F = generate_fibonacci(50);  % Generate first 50 Fibonacci numbers
    
    n = 1;
    while (1/F(n+1)) > measure
        n = n + 1;
        if n >= length(F) - 1
            error('Need more Fibonacci numbers. Increase limit.');
        end
    end
    
    if verbose
        fprintf('Fibonacci Number Selection:\n');
        fprintf('---------------------------\n');
        fprintf('Required: 1/F_n <= %.4f\n', measure);
        fprintf('Selected: n = %d, F_%d = %d\n', n, n, F(n+1));
        fprintf('Achieved measure: 1/F_%d = %.6f\n\n', n, 1/F(n+1));
        
        fprintf('Fibonacci Sequence (first %d terms):\n', min(15, n+2));
        fprintf('n:  ');
        for i = 0:min(14, n+1)
            fprintf('%4d ', i);
        end
        fprintf('\nF_n:');
        for i = 1:min(15, n+2)
            fprintf('%4d ', F(i));
        end
        fprintf('\n\n');
    end
    
    % Step 3-6: Iterative search
    if verbose
        fprintf('========================================\n');
        fprintf('Fibonacci Search Iterations\n');
        fprintf('========================================\n\n');
    end
    
    for i = n:-1:2
        if verbose
            fprintf('Iteration %d (i=%d):\n', n-i+1, i);
            fprintf('----------------\n');
            fprintf('Current interval: [%.6f, %.6f], length = %.6f\n', ...
                    a, b, b-a);
        end
        
        % Step 4: Calculate test points
        L_current = b - a;
        x1 = a + (F(i-1) / F(i+1)) * L_current;
        x2 = a + (F(i) / F(i+1)) * L_current;
        
        % Evaluate function at test points
        f1 = f(x1);
        f2 = f(x2);
        
        if verbose
            fprintf('Test points:\n');
            fprintf('  x1 = a + F_%d/F_%d * L = %.6f, f(x1) = %.6f\n', ...
                    i-2, i, x1, f1);
            fprintf('  x2 = a + F_%d/F_%d * L = %.6f, f(x2) = %.6f\n', ...
                    i-1, i, x2, f2);
        end
        
        % Step 5: Update interval based on optimization type
        if strcmp(opt_type, 'min')
            if f1 > f2
                % Minimum is in [x1, b]
                a = x1;
                if verbose
                    fprintf('Decision: f(x1) > f(x2) → New interval: [%.6f, %.6f]\n', a, b);
                end
            else
                % Minimum is in [a, x2]
                b = x2;
                if verbose
                    fprintf('Decision: f(x1) <= f(x2) → New interval: [%.6f, %.6f]\n', a, b);
                end
            end
        else  % 'max'
            if f1 > f2
                % Maximum is in [a, x2]
                b = x2;
                if verbose
                    fprintf('Decision: f(x1) > f(x2) → New interval: [%.6f, %.6f]\n', a, b);
                end
            else
                % Maximum is in [x1, b]
                a = x1;
                if verbose
                    fprintf('Decision: f(x1) <= f(x2) → New interval: [%.6f, %.6f]\n', a, b);
                end
            end
        end
        
        if verbose
            fprintf('\n');
        end
    end
    
    % Final result
    x_opt = (a + b) / 2;
    f_opt = f(x_opt);
    final_interval = [a, b];
    
    if verbose
        fprintf('========================================\n');
        fprintf('Search Complete!\n');
        fprintf('========================================\n\n');
    end
end

function F = generate_fibonacci(n)
    % Generate first n Fibonacci numbers
    % F(1) = F_0 = 1, F(2) = F_1 = 1, ...
    
    F = zeros(n, 1);
    F(1) = 1;  % F_0
    F(2) = 1;  % F_1
    
    for i = 3:n
        F(i) = F(i-1) + F(i-2);
    end
end

function plot_fibonacci_search(f, a, b, x_opt)
    % Plot the objective function and optimal point
    
    figure('Position', [100, 100, 800, 500]);
    
    % Create fine grid for plotting
    x = linspace(a - 0.2, b + 0.2, 1000);
    y = arrayfun(f, x);
    
    % Plot function
    plot(x, y, 'b-', 'LineWidth', 2);
    hold on;
    grid on;
    
    % Mark the initial interval
    y_min = min(y);
    y_max = max(y);
    plot([a a], [y_min, y_max], 'g--', 'LineWidth', 1.5, ...
         'DisplayName', 'Initial interval bounds');
    plot([b b], [y_min, y_max], 'g--', 'LineWidth', 1.5, ...
         'HandleVisibility', 'off');
    
    % Mark the optimal point
    f_opt = f(x_opt);
    plot(x_opt, f_opt, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r', ...
         'DisplayName', sprintf('Optimal point (%.4f, %.4f)', x_opt, f_opt));
    
    % Add vertical line at optimal point
    plot([x_opt x_opt], [y_min, f_opt], 'r--', 'LineWidth', 1, ...
         'HandleVisibility', 'off');
    
    xlabel('x', 'FontSize', 12);
    ylabel('f(x)', 'FontSize', 12);
    title('Fibonacci Search - Objective Function', 'FontSize', 14);
    legend('Location', 'best');
    
    % Add text annotation
    text(x_opt, f_opt + 0.1*(y_max-y_min), ...
         sprintf('  x* = %.4f\n  f(x*) = %.4f', x_opt, f_opt), ...
         'FontSize', 10, 'Color', 'red');
    
    hold off;
end
