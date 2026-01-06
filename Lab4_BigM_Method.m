% =========================================================================
% Thapar Institute of Engineering and Technology, Patiala
% School of Mathematics
% Optimization Methods (UMA-035)
% Lab Experiment - 4: Big-M Method
% =========================================================================

function Lab4_BigM_Method()
    clc;
    clear;
    
    fprintf('========================================\n');
    fprintf('Big-M Method for LPP\n');
    fprintf('========================================\n\n');
    
    fprintf('Select Problem:\n');
    fprintf('1. Min z = 3x1 + 5x2, s.t. x1+3x2>=3, x1+x2>=2\n');
    fprintf('2. Min z = 12x1 + 10x2, s.t. 5x1+x2>=10, 6x1+5x2>=30, x1+4x2>=8\n');
    fprintf('3. Max z = 3x1 + 2x2, s.t. x1+x2<=2, x1+3x2<=3, x1-x2=1\n');
    problem = input('Enter problem number (1-3): ');
    
    switch problem
        case 1
            problem1_bigM();
        case 2
            problem2_bigM();
        case 3
            problem3_bigM();
        otherwise
            error('Invalid problem selection');
    end
end

function problem1_bigM()
    % Min z = 3x1 + 5x2
    % x1 + 3x2 >= 3
    % x1 + x2 >= 2
    % Convert to standard form:
    % x1 + 3x2 - s1 + R1 = 3
    % x1 + x2 - s2 + R2 = 2
    
    fprintf('\nProblem 1: Min z = 3x1 + 5x2\n');
    fprintf('----------------------------\n');
    
    % Maximize -z = -3x1 - 5x2 (convert min to max)
    A = [1, 3, -1, 0, 1, 0;
         1, 1,  0, -1, 0, 1];
    b = [3; 2];
    
    M = 1e6;
    c = [-3; -5; 0; 0; -M; -M];  % For maximization
    
    % Initial basis: [R1, R2] = [5, 6]
    basis = [5, 6];
    
    [x_opt, z_opt, feasible] = bigM_method(A, b, c, basis, M, true);
    
    if feasible
        fprintf('\n========================================\n');
        fprintf('OPTIMAL SOLUTION (Minimization)\n');
        fprintf('========================================\n');
        fprintf('x1 = %.4f\n', x_opt(1));
        fprintf('x2 = %.4f\n', x_opt(2));
        fprintf('Optimal value: z = %.4f\n', -z_opt);  % Convert back to min
    else
        fprintf('\nProblem is INFEASIBLE!\n');
    end
end

function problem2_bigM()
    % Min z = 12x1 + 10x2
    % 5x1 + x2 >= 10
    % 6x1 + 5x2 >= 30
    % x1 + 4x2 >= 8
    
    fprintf('\nProblem 2: Min z = 12x1 + 10x2\n');
    fprintf('-------------------------------\n');
    
    A = [5, 1, -1, 0, 0, 1, 0, 0;
         6, 5, 0, -1, 0, 0, 1, 0;
         1, 4, 0, 0, -1, 0, 0, 1];
    b = [10; 30; 8];
    
    M = 1e6;
    c = [-12; -10; 0; 0; 0; -M; -M; -M];
    
    basis = [6, 7, 8];
    
    [x_opt, z_opt, feasible] = bigM_method(A, b, c, basis, M, true);
    
    if feasible
        fprintf('\n========================================\n');
        fprintf('OPTIMAL SOLUTION (Minimization)\n');
        fprintf('========================================\n');
        fprintf('x1 = %.4f\n', x_opt(1));
        fprintf('x2 = %.4f\n', x_opt(2));
        fprintf('Optimal value: z = %.4f\n', -z_opt);
    else
        fprintf('\nProblem is INFEASIBLE!\n');
    end
end

function problem3_bigM()
    % Max z = 3x1 + 2x2
    % x1 + x2 <= 2
    % x1 + 3x2 <= 3
    % x1 - x2 = 1
    
    fprintf('\nProblem 3: Max z = 3x1 + 2x2\n');
    fprintf('-----------------------------\n');
    
    A = [1, 1, 1, 0, 0;
         1, 3, 0, 1, 0;
         1, -1, 0, 0, 1];  % Artificial variable for equality
    b = [2; 3; 1];
    
    M = 1e6;
    c = [3; 2; 0; 0; -M];
    
    basis = [3, 4, 5];
    
    [x_opt, z_opt, feasible] = bigM_method(A, b, c, basis, M, true);
    
    if feasible
        fprintf('\n========================================\n');
        fprintf('OPTIMAL SOLUTION (Maximization)\n');
        fprintf('========================================\n');
        fprintf('x1 = %.4f\n', x_opt(1));
        fprintf('x2 = %.4f\n', x_opt(2));
        fprintf('Optimal value: z = %.4f\n', z_opt);
    else
        fprintf('\nProblem is INFEASIBLE!\n');
    end
end

function [x_opt, z_opt, feasible] = bigM_method(A, b, c, basis, M, verbose)
    % Big-M Method Algorithm
    
    [m, n] = size(A);
    iteration = 0;
    max_iterations = 100;
    
    % Identify artificial variables (those with cost -M)
    artificial_vars = find(abs(c + M) < 1e-3);
    
    if verbose
        fprintf('\nStarting Big-M Method...\n');
        fprintf('========================\n');
        fprintf('M = %.0e\n', M);
        fprintf('Artificial variables: ');
        fprintf('x%d ', artificial_vars);
        fprintf('\n\n');
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
            fprintf('Objective value: z = %.4e\n', z);
        end
        
        % Calculate Zj - cj
        Zj_minus_cj = zeros(n, 1);
        for j = 1:n
            Aj = A(:, j);
            Zj_minus_cj(j) = CB' * (B_inv * Aj) - c(j);
        end
        
        if verbose
            fprintf('Zj - cj: ');
            for j = 1:n
                fprintf('%.2e ', Zj_minus_cj(j));
            end
            fprintf('\n');
        end
        
        % Find minimum Zj - cj
        [min_val, k] = min(Zj_minus_cj);
        
        % Check optimality
        if min_val >= -1e-6
            if verbose
                fprintf('\nOptimal solution found!\n');
            end
            
            % Check if any artificial variable is in basis with positive value
            for i = 1:m
                if ismember(basis(i), artificial_vars) && XB(i) > 1e-6
                    if verbose
                        fprintf('Artificial variable x%d = %.4f in basis\n', ...
                                basis(i), XB(i));
                    end
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
        
        if verbose
            fprintf('Entering variable: x%d\n', k);
        end
        
        % Calculate alpha_k
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
        
        % Calculate minimum ratio
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
            fprintf('Leaving variable: x%d (row %d)\n\n', basis(leaving_row), leaving_row);
        end
        
        % Update basis
        basis(leaving_row) = k;
    end
    
    warning('Maximum iterations reached');
    x_opt = X;
    z_opt = z;
    feasible = false;
end
