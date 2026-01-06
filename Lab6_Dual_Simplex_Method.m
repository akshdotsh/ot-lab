% =========================================================================
% Thapar Institute of Engineering and Technology, Patiala
% School of Mathematics
% Optimization Methods (UMA-035)
% Lab Experiment - 6: Dual Simplex Method
% =========================================================================

function Lab6_Dual_Simplex_Method()
    clc;
    clear;
    
    fprintf('========================================\n');
    fprintf('Dual Simplex Method for LPP\n');
    fprintf('========================================\n\n');
    
    fprintf('Select Problem:\n');
    fprintf('1. Min z = 3x1 + 5x2, s.t. x1+3x2>=3, x1+x2>=2\n');
    fprintf('2. Min z = 12x1 + 10x2, s.t. 5x1+x2>=10, 6x1+5x2>=30, x1+4x2>=8\n');
    fprintf('3. Min z = 3x1 + 2x2, s.t. x1+x2<=1, x1+2x2>=3\n');
    problem = input('Enter problem number (1-3): ');
    
    switch problem
        case 1
            problem1_dual_simplex();
        case 2
            problem2_dual_simplex();
        case 3
            problem3_dual_simplex();
        otherwise
            error('Invalid problem selection');
    end
end

function problem1_dual_simplex()
    % Min z = 3x1 + 5x2
    % x1 + 3x2 >= 3  =>  -x1 - 3x2 + s1 = -3
    % x1 + x2 >= 2   =>  -x1 - x2 + s2 = -2
    
    fprintf('\nProblem 1: Min z = 3x1 + 5x2\n');
    fprintf('----------------------------\n');
    
    A = [-1, -3, 1, 0;
         -1, -1, 0, 1];
    b = [-3; -2];
    c = [3; 5; 0; 0];
    
    % Initial basis: [s1, s2] = [3, 4]
    basis = [3, 4];
    
    [x_opt, z_opt, feasible] = dual_simplex_method(A, b, c, basis, true);
    
    if feasible
        fprintf('\n========================================\n');
        fprintf('OPTIMAL SOLUTION (Minimization)\n');
        fprintf('========================================\n');
        fprintf('x1 = %.4f\n', x_opt(1));
        fprintf('x2 = %.4f\n', x_opt(2));
        fprintf('Optimal value: z = %.4f\n', z_opt);
    else
        fprintf('\nProblem is INFEASIBLE!\n');
    end
end

function problem2_dual_simplex()
    % Min z = 12x1 + 10x2
    % 5x1 + x2 >= 10   =>  -5x1 - x2 + s1 = -10
    % 6x1 + 5x2 >= 30  =>  -6x1 - 5x2 + s2 = -30
    % x1 + 4x2 >= 8    =>  -x1 - 4x2 + s3 = -8
    
    fprintf('\nProblem 2: Min z = 12x1 + 10x2\n');
    fprintf('-------------------------------\n');
    
    A = [-5, -1, 1, 0, 0;
         -6, -5, 0, 1, 0;
         -1, -4, 0, 0, 1];
    b = [-10; -30; -8];
    c = [12; 10; 0; 0; 0];
    
    basis = [3, 4, 5];
    
    [x_opt, z_opt, feasible] = dual_simplex_method(A, b, c, basis, true);
    
    if feasible
        fprintf('\n========================================\n');
        fprintf('OPTIMAL SOLUTION (Minimization)\n');
        fprintf('========================================\n');
        fprintf('x1 = %.4f\n', x_opt(1));
        fprintf('x2 = %.4f\n', x_opt(2));
        fprintf('Optimal value: z = %.4f\n', z_opt);
    else
        fprintf('\nProblem is INFEASIBLE!\n');
    end
end

function problem3_dual_simplex()
    % Min z = 3x1 + 2x2
    % x1 + x2 <= 1     =>  x1 + x2 + s1 = 1
    % x1 + 2x2 >= 3    =>  -x1 - 2x2 + s2 = -3
    
    fprintf('\nProblem 3: Min z = 3x1 + 2x2\n');
    fprintf('-----------------------------\n');
    
    A = [1, 1, 1, 0;
        -1, -2, 0, 1];
    b = [1; -3];
    c = [3; 2; 0; 0];
    
    basis = [3, 4];
    
    [x_opt, z_opt, feasible] = dual_simplex_method(A, b, c, basis, true);
    
    if feasible
        fprintf('\n========================================\n');
        fprintf('OPTIMAL SOLUTION (Minimization)\n');
        fprintf('========================================\n');
        fprintf('x1 = %.4f\n', x_opt(1));
        fprintf('x2 = %.4f\n', x_opt(2));
        fprintf('Optimal value: z = %.4f\n', z_opt);
    else
        fprintf('\nProblem is INFEASIBLE!\n');
    end
end

function [x_opt, z_opt, feasible] = dual_simplex_method(A, b, c, basis, verbose)
    % Dual Simplex Method Algorithm
    % Assumes: Initial solution is optimal but infeasible (some b_i < 0)
    
    [m, n] = size(A);
    iteration = 0;
    max_iterations = 100;
    
    if verbose
        fprintf('\nStarting Dual Simplex Method...\n');
        fprintf('================================\n\n');
    end
    
    while iteration < max_iterations
        iteration = iteration + 1;
        
        % Construct basis matrix
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
        
        % Calculate Zj - cj for all variables
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
        
        % Check if all Zj - cj >= 0 (optimality for minimization)
        if all(Zj_minus_cj >= -1e-10)
            if verbose
                fprintf('Optimal conditions satisfied!\n');
            end
            
            % Check feasibility: all XB >= 0
            if all(XB >= -1e-10)
                if verbose
                    fprintf('Solution is FEASIBLE!\n');
                end
                x_opt = X;
                z_opt = z;
                feasible = true;
                return;
            end
        end
        
        % Step 1: Select leaving variable (most negative XB)
        [min_XB, leaving_row] = min(XB);
        
        if min_XB >= -1e-10
            if verbose
                fprintf('All basic variables are non-negative!\n');
            end
            x_opt = X;
            z_opt = z;
            feasible = true;
            return;
        end
        
        leaving_var = basis(leaving_row);
        
        if verbose
            fprintf('Leaving variable: x%d (XB = %.4f)\n', leaving_var, min_XB);
        end
        
        % Calculate the leaving row of the tableau
        Y_r = B_inv(leaving_row, :) * A;
        
        if verbose
            fprintf('Y_r (pivot row): ');
            fprintf('%.4f ', Y_r);
            fprintf('\n');
        end
        
        % Step 2: Select entering variable using minimum ratio test
        % Ratio: |Zj - cj| / |Y_rj| for Y_rj < 0
        
        min_ratio = inf;
        entering_var = -1;
        
        for j = 1:n
            if Y_r(j) < -1e-10  % Only consider negative Y_rj
                ratio = abs(Zj_minus_cj(j)) / abs(Y_r(j));
                if ratio < min_ratio
                    min_ratio = ratio;
                    entering_var = j;
                end
            end
        end
        
        if entering_var == -1
            if verbose
                fprintf('No entering variable found!\n');
                fprintf('Problem is INFEASIBLE!\n');
            end
            x_opt = X;
            z_opt = z;
            feasible = false;
            return;
        end
        
        if verbose
            fprintf('Entering variable: x%d (ratio = %.4f)\n', entering_var, min_ratio);
        end
        
        % Step 3: Update basis
        basis(leaving_row) = entering_var;
        
        if verbose
            fprintf('\n');
        end
    end
    
    warning('Maximum iterations reached');
    x_opt = X;
    z_opt = z;
    feasible = false;
end
