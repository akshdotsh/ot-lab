% =========================================================================
% Thapar Institute of Engineering and Technology, Patiala
% School of Mathematics
% Optimization Methods (UMA-035)
% Lab Experiment - 8: Multi-Objective LPP (Weighted Sum Method)
% =========================================================================

function Lab8_MultiObjective_LPP()
    clc;
    clear;
    
    fprintf('========================================\n');
    fprintf('Multi-Objective LPP - Weighted Sum Method\n');
    fprintf('========================================\n\n');
    
    fprintf('Select Problem:\n');
    fprintf('1. Max (3x1+2x2+4x3) and Max (x1+5x2+3x3)\n');
    fprintf('2. Max (x1+4x2+x3) and Max (2x1+7x2+5x3)\n');
    problem = input('Enter problem number (1-2): ');
    
    switch problem
        case 1
            problem1_multiobjective();
        case 2
            problem2_multiobjective();
        otherwise
            error('Invalid problem selection');
    end
end

function problem1_multiobjective()
    % Problem 1:
    % Maximize (3x1 + 2x2 + 4x3)
    % Maximize (x1 + 5x2 + 3x3)
    % Subject to:
    % 2x1 + 4x2 + x3 <= 8
    % 3x1 + 5x2 + 4x3 >= 15
    % x1, x2, x3 >= 0
    
    fprintf('\nProblem 1: Multi-Objective Optimization\n');
    fprintf('---------------------------------------\n');
    fprintf('Objective 1: Max z1 = 3x1 + 2x2 + 4x3\n');
    fprintf('Objective 2: Max z2 = x1 + 5x2 + 3x3\n');
    fprintf('Subject to:\n');
    fprintf('  2x1 + 4x2 + x3 <= 8\n');
    fprintf('  3x1 + 5x2 + 4x3 >= 15\n');
    fprintf('  x1, x2, x3 >= 0\n\n');
    
    % Objective coefficients
    C1 = [3; 2; 4];
    C2 = [1; 5; 3];
    
    % Weighted sum: (C1 + C2) / 2
    C_combined = (C1 + C2) / 2;
    
    fprintf('Combined objective (weighted average):\n');
    fprintf('Max z = (z1 + z2)/2 = %.2fx1 + %.2fx2 + %.2fx3\n\n', ...
            C_combined(1), C_combined(2), C_combined(3));
    
    % Convert to standard form with slack/surplus/artificial variables
    % 2x1 + 4x2 + x3 + s1 = 8
    % 3x1 + 5x2 + 4x3 - s2 + R1 = 15
    
    A = [2, 4, 1, 1, 0, 0;
         3, 5, 4, 0, -1, 1];
    b = [8; 15];
    
    M = 1e6;
    c = [C_combined; 0; 0; -M];  % Include slack, surplus, artificial
    
    % Initial basis: [s1, R1] = [4, 6]
    basis = [4, 6];
    
    [x_opt, z_opt, feasible] = solve_bigM(A, b, c, basis, M, true);
    
    if feasible
        fprintf('\n========================================\n');
        fprintf('OPTIMAL SOLUTION\n');
        fprintf('========================================\n');
        fprintf('x1 = %.4f\n', x_opt(1));
        fprintf('x2 = %.4f\n', x_opt(2));
        fprintf('x3 = %.4f\n', x_opt(3));
        
        % Calculate individual objective values
        z1 = C1' * x_opt(1:3);
        z2 = C2' * x_opt(1:3);
        
        fprintf('\nObjective Values:\n');
        fprintf('z1 = %.4f\n', z1);
        fprintf('z2 = %.4f\n', z2);
        fprintf('Combined: z = %.4f\n', z_opt);
    else
        fprintf('\nProblem is INFEASIBLE!\n');
    end
end

function problem2_multiobjective()
    % Problem 2:
    % Maximize (x1 + 4x2 + x3)
    % Maximize (2x1 + 7x2 + 5x3)
    % Subject to:
    % x1 + x2 + x3 <= 8
    % x1 + 5x2 + 4x3 >= 15
    % x1, x2, x3 >= 0
    
    fprintf('\nProblem 2: Multi-Objective Optimization\n');
    fprintf('---------------------------------------\n');
    fprintf('Objective 1: Max z1 = x1 + 4x2 + x3\n');
    fprintf('Objective 2: Max z2 = 2x1 + 7x2 + 5x3\n');
    fprintf('Subject to:\n');
    fprintf('  x1 + x2 + x3 <= 8\n');
    fprintf('  x1 + 5x2 + 4x3 >= 15\n');
    fprintf('  x1, x2, x3 >= 0\n\n');
    
    C1 = [1; 4; 1];
    C2 = [2; 7; 5];
    
    C_combined = (C1 + C2) / 2;
    
    fprintf('Combined objective (weighted average):\n');
    fprintf('Max z = (z1 + z2)/2 = %.2fx1 + %.2fx2 + %.2fx3\n\n', ...
            C_combined(1), C_combined(2), C_combined(3));
    
    A = [1, 1, 1, 1, 0, 0;
         1, 5, 4, 0, -1, 1];
    b = [8; 15];
    
    M = 1e6;
    c = [C_combined; 0; 0; -M];
    
    basis = [4, 6];
    
    [x_opt, z_opt, feasible] = solve_bigM(A, b, c, basis, M, true);
    
    if feasible
        fprintf('\n========================================\n');
        fprintf('OPTIMAL SOLUTION\n');
        fprintf('========================================\n');
        fprintf('x1 = %.4f\n', x_opt(1));
        fprintf('x2 = %.4f\n', x_opt(2));
        fprintf('x3 = %.4f\n', x_opt(3));
        
        z1 = C1' * x_opt(1:3);
        z2 = C2' * x_opt(1:3);
        
        fprintf('\nObjective Values:\n');
        fprintf('z1 = %.4f\n', z1);
        fprintf('z2 = %.4f\n', z2);
        fprintf('Combined: z = %.4f\n', z_opt);
    else
        fprintf('\nProblem is INFEASIBLE!\n');
    end
end

function [x_opt, z_opt, feasible] = solve_bigM(A, b, c, basis, M, verbose)
    % Simplified Big-M solver for multi-objective problems
    
    [m, n] = size(A);
    iteration = 0;
    max_iterations = 100;
    
    artificial_vars = find(abs(c + M) < 1e-3);
    
    if verbose
        fprintf('Solving using Big-M Method...\n');
        fprintf('==============================\n\n');
    end
    
    while iteration < max_iterations
        iteration = iteration + 1;
        
        B = A(:, basis);
        B_inv = inv(B);
        XB = B_inv * b;
        
        X = zeros(n, 1);
        X(basis) = XB;
        
        CB = c(basis);
        z = CB' * XB;
        
        if verbose && mod(iteration, 5) == 1
            fprintf('Iteration %d: z = %.4e\n', iteration, z);
        end
        
        Zj_minus_cj = zeros(n, 1);
        for j = 1:n
            Aj = A(:, j);
            Zj_minus_cj(j) = CB' * (B_inv * Aj) - c(j);
        end
        
        [min_val, k] = min(Zj_minus_cj);
        
        if min_val >= -1e-10
            % Check for artificial variables in basis
            for i = 1:m
                if ismember(basis(i), artificial_vars) && XB(i) > 1e-6
                    x_opt = X;
                    z_opt = z;
                    feasible = false;
                    return;
                end
            end
            
            x_opt = X;
            z_opt = z;
            feasible = true;
            return;
        end
        
        alpha_k = B_inv * A(:, k);
        
        if all(alpha_k <= 1e-10)
            warning('Problem is unbounded');
            x_opt = X;
            z_opt = inf;
            feasible = false;
            return;
        end
        
        theta = inf;
        leaving_row = -1;
        for i = 1:m
            if alpha_k(i) > 1e-10
                ratio = XB(i) / alpha_k(i);
                if ratio < theta
                    theta = ratio;
                    leaving_row = i;
                end
            end
        end
        
        basis(leaving_row) = k;
    end
    
    x_opt = X;
    z_opt = z;
    feasible = false;
end
