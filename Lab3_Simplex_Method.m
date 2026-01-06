% =========================================================================
% Thapar Institute of Engineering and Technology, Patiala
% School of Mathematics
% Optimization Methods (UMA-035)
% Lab Experiment - 3: Simplex Method
% =========================================================================

function Lab3_Simplex_Method()
    clc;
    clear;
    
    fprintf('========================================\n');
    fprintf('Simplex Method for LPP\n');
    fprintf('========================================\n\n');
    
    fprintf('Select Problem:\n');
    fprintf('1. Max z = x1 + 2x2, s.t. -x1+x2<=1, x1+x2<=2\n');
    fprintf('2. Max z = 4x1+6x2+3x3+x4, s.t. (3 constraints)\n');
    fprintf('3. Min z = -3/4*x4+20*x5-1/2*x6+6*x7 (with equality constraints)\n');
    problem = input('Enter problem number (1-3): ');
    
    switch problem
        case 1
            problem1_simplex();
        case 2
            problem2_simplex();
        case 3
            problem3_simplex();
        otherwise
            error('Invalid problem selection');
    end
end

function problem1_simplex()
    % Max z = x1 + 2x2
    % -x1 + x2 <= 1
    %  x1 + x2 <= 2
    % x1, x2 >= 0
    
    fprintf('\nProblem 1: Max z = x1 + 2x2\n');
    fprintf('---------------------------\n');
    
    % Standard form with slack variables
    % Max z = x1 + 2x2 + 0*s1 + 0*s2
    % -x1 + x2 + s1 = 1
    %  x1 + x2 + s2 = 2
    
    A = [-1, 1, 1, 0;
          1, 1, 0, 1];
    b = [1; 2];
    c = [1; 2; 0; 0];
    
    % Initial basis: [s1, s2] = [3, 4]
    basis = [3, 4];
    
    [x_opt, z_opt, iterations] = simplex_method(A, b, c, basis, true);
    
    fprintf('\n========================================\n');
    fprintf('OPTIMAL SOLUTION\n');
    fprintf('========================================\n');
    fprintf('x1 = %.4f\n', x_opt(1));
    fprintf('x2 = %.4f\n', x_opt(2));
    fprintf('Optimal value: z = %.4f\n', z_opt);
    fprintf('Iterations: %d\n', iterations);
end

function problem2_simplex()
    % Max z = 4x1 + 6x2 + 3x3 + x4
    % x1 + 4x2 + 8x3 + 6x4 <= 11
    % 4x1 + x2 + 2x3 + x4 <= 7
    % 2x1 + 3x2 + x3 + 2x4 <= 2
    
    fprintf('\nProblem 2: Max z = 4x1 + 6x2 + 3x3 + x4\n');
    fprintf('----------------------------------------\n');
    
    A = [1, 4, 8, 6, 1, 0, 0;
         4, 1, 2, 1, 0, 1, 0;
         2, 3, 1, 2, 0, 0, 1];
    b = [11; 7; 2];
    c = [4; 6; 3; 1; 0; 0; 0];
    
    % Initial basis: [s1, s2, s3] = [5, 6, 7]
    basis = [5, 6, 7];
    
    [x_opt, z_opt, iterations] = simplex_method(A, b, c, basis, true);
    
    fprintf('\n========================================\n');
    fprintf('OPTIMAL SOLUTION\n');
    fprintf('========================================\n');
    fprintf('x1 = %.4f\n', x_opt(1));
    fprintf('x2 = %.4f\n', x_opt(2));
    fprintf('x3 = %.4f\n', x_opt(3));
    fprintf('x4 = %.4f\n', x_opt(4));
    fprintf('s1 = %.4f\n', x_opt(5));
    fprintf('s2 = %.4f\n', x_opt(6));
    fprintf('s3 = %.4f\n', x_opt(7));
    fprintf('Optimal value: z = %.4f\n', z_opt);
    fprintf('Expected: x1=1/3, x3=4/3, s2=3, z=16/3\n');
end

function problem3_simplex()
    % Min z = -3/4*x4 + 20*x5 - 1/2*x6 + 6*x7
    % Convert to Max: -z = 3/4*x4 - 20*x5 + 1/2*x6 - 6*x7
    % x1 + 1/4*x4 - 8*x5 - x6 + 9*x7 = 0
    % x2 + 1/2*x4 - 12*x5 - 1/6*x6 + 3*x7 = 0
    % x3 + x6 = 1
    
    fprintf('\nProblem 3: Min z = -3/4*x4 + 20*x5 - 1/2*x6 + 6*x7\n');
    fprintf('---------------------------------------------------\n');
    
    A = [1, 0, 0, 1/4, -8, -1, 9;
         0, 1, 0, 1/2, -12, -1/6, 3;
         0, 0, 1, 0, 0, 1, 0];
    b = [0; 0; 1];
    c = [-3/4; -20; 1/2; -6; 0; 0; 0];  % For minimization (negate for max)
    
    % Initial basis: [x1, x2, x3] = [1, 2, 3]
    basis = [1, 2, 3];
    
    [x_opt, z_opt, iterations] = simplex_method(A, b, c, basis, true);
    
    fprintf('\n========================================\n');
    fprintf('OPTIMAL SOLUTION (Minimization)\n');
    fprintf('========================================\n');
    fprintf('x1 = %.4f\n', x_opt(1));
    fprintf('x2 = %.4f\n', x_opt(2));
    fprintf('x3 = %.4f\n', x_opt(3));
    fprintf('x4 = %.4f\n', x_opt(4));
    fprintf('x5 = %.4f\n', x_opt(5));
    fprintf('x6 = %.4f\n', x_opt(6));
    fprintf('x7 = %.4f\n', x_opt(7));
    fprintf('Optimal value: z = %.4f\n', -z_opt);  % Convert back to min
    fprintf('Expected: x1=3/4, x4=1, x6=1, z=-5/4\n');
end

function [x_opt, z_opt, iterations] = simplex_method(A, b, c, basis, verbose)
    % Simplex Method Algorithm
    % Inputs:
    %   A: coefficient matrix (m x n)
    %   b: RHS vector (m x 1)
    %   c: cost vector (n x 1)
    %   basis: initial basis indices
    %   verbose: display iterations (true/false)
    
    [m, n] = size(A);
    iteration = 0;
    max_iterations = 100;
    
    if verbose
        fprintf('\nStarting Simplex Method...\n');
        fprintf('==========================\n\n');
    end
    
    while iteration < max_iterations
        iteration = iteration + 1;
        
        % Step 1: Construct basis matrix B
        B = A(:, basis);
        B_inv = inv(B);
        
        % Calculate basic solution
        XB = B_inv * b;
        
        % Create full solution
        X = zeros(n, 1);
        X(basis) = XB;
        
        % Calculate objective value
        CB = c(basis);
        z = CB' * XB;
        
        if verbose
            fprintf('Iteration %d:\n', iteration);
            fprintf('-------------\n');
            fprintf('Basis: ');
            fprintf('x%d ', basis);
            fprintf('\n');
            fprintf('Basic solution: ');
            fprintf('%.4f ', XB);
            fprintf('\n');
            fprintf('Objective value: z = %.4f\n', z);
        end
        
        % Step 2: Calculate Zj - cj for all variables
        Zj_minus_cj = zeros(n, 1);
        for j = 1:n
            Aj = A(:, j);
            Zj_minus_cj(j) = CB' * (B_inv * Aj) - c(j);
        end
        
        if verbose
            fprintf('Zj - cj: ');
            fprintf('%.4f ', Zj_minus_cj);
            fprintf('\n');
        end
        
        % Find minimum Zj - cj
        [min_val, k] = min(Zj_minus_cj);
        
        % Check optimality
        if min_val >= -1e-10
            if verbose
                fprintf('\nOptimal solution found!\n');
            end
            x_opt = X;
            z_opt = z;
            iterations = iteration;
            return;
        end
        
        if verbose
            fprintf('Entering variable: x%d (Zj-cj = %.4f)\n', k, min_val);
        end
        
        % Calculate alpha_k (column k in simplex tableau)
        alpha_k = B_inv * A(:, k);
        
        if verbose
            fprintf('Alpha_k: ');
            fprintf('%.4f ', alpha_k);
            fprintf('\n');
        end
        
        % Check for unboundedness
        if all(alpha_k <= 1e-10)
            error('Problem is unbounded!');
        end
        
        % Calculate minimum ratio (theta)
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
        
        if verbose
            fprintf('Leaving variable: x%d (row %d, theta = %.4f)\n', ...
                    basis(leaving_row), leaving_row, theta);
        end
        
        % Step 3: Update basis
        basis(leaving_row) = k;
        
        if verbose
            fprintf('\n');
        end
    end
    
    warning('Maximum iterations reached without convergence');
    x_opt = X;
    z_opt = z;
    iterations = iteration;
end
