% =========================================================================
% Thapar Institute of Engineering and Technology, Patiala
% School of Mathematics
% Optimization Methods (UMA-035)
% Lab Experiment - 5: Two-Phase Method
% =========================================================================

function Lab5_TwoPhase_Method()
    clc;
    clear;
    
    fprintf('========================================\n');
    fprintf('Two-Phase Method for LPP\n');
    fprintf('========================================\n\n');
    
    fprintf('Select Problem:\n');
    fprintf('1. Min z = 3x1 + 5x2, s.t. x1+3x2>=3, x1+x2>=2\n');
    fprintf('2. Min z = 12x1 + 10x2, s.t. 5x1+x2>=10, 6x1+5x2>=30, x1+4x2>=8\n');
    fprintf('3. Max z = 3x1 + 2x2, s.t. x1+x2<=2, x1+3x2<=3, x1-x2=1\n');
    problem = input('Enter problem number (1-3): ');
    
    switch problem
        case 1
            problem1_twophase();
        case 2
            problem2_twophase();
        case 3
            problem3_twophase();
        otherwise
            error('Invalid problem selection');
    end
end

function problem1_twophase()
    fprintf('\nProblem 1: Min z = 3x1 + 5x2\n');
    fprintf('----------------------------\n');
    
    % Standard form with surplus and artificial variables
    A = [1, 3, -1, 0, 1, 0;
         1, 1,  0, -1, 0, 1];
    b = [3; 2];
    c_original = [-3; -5; 0; 0; 0; 0];  % For max (negate for min)
    
    num_artificial = 2;
    artificial_vars = [5, 6];
    
    [x_opt, z_opt, feasible] = two_phase_method(A, b, c_original, ...
                                                 artificial_vars, true);
    
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

function problem2_twophase()
    fprintf('\nProblem 2: Min z = 12x1 + 10x2\n');
    fprintf('-------------------------------\n');
    
    A = [5, 1, -1, 0, 0, 1, 0, 0;
         6, 5, 0, -1, 0, 0, 1, 0;
         1, 4, 0, 0, -1, 0, 0, 1];
    b = [10; 30; 8];
    c_original = [-12; -10; 0; 0; 0; 0; 0; 0];
    
    artificial_vars = [6, 7, 8];
    
    [x_opt, z_opt, feasible] = two_phase_method(A, b, c_original, ...
                                                 artificial_vars, true);
    
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

function problem3_twophase()
    fprintf('\nProblem 3: Max z = 3x1 + 2x2\n');
    fprintf('-----------------------------\n');
    
    A = [1, 1, 1, 0, 0;
         1, 3, 0, 1, 0;
         1, -1, 0, 0, 1];
    b = [2; 3; 1];
    c_original = [3; 2; 0; 0; 0];
    
    artificial_vars = [5];  % Only equality constraint needs artificial
    
    [x_opt, z_opt, feasible] = two_phase_method(A, b, c_original, ...
                                                 artificial_vars, true);
    
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

function [x_opt, z_opt, feasible] = two_phase_method(A, b, c_original, ...
                                                      artificial_vars, verbose)
    % Two-Phase Method Algorithm
    
    [m, n] = size(A);
    
    %% PHASE I
    if verbose
        fprintf('\n========================================\n');
        fprintf('PHASE I: Finding Initial BFS\n');
        fprintf('========================================\n\n');
    end
    
    % Create Phase I objective: Min sum of artificial variables
    % For maximization: Max (-sum of artificial variables)
    c_phase1 = zeros(n, 1);
    c_phase1(artificial_vars) = -1;  % Maximize negative sum (minimize sum)
    
    % Initial basis: artificial variables
    basis = artificial_vars;
    
    % Solve Phase I
    [x_phase1, z_phase1, ~] = simplex_solve(A, b, c_phase1, basis, verbose);
    
    % Check if Phase I found a feasible solution
    if abs(z_phase1) > 1e-6  % If min sum of artificial > 0
        if verbose
            fprintf('\nPhase I Result: z = %.4e\n', z_phase1);
            fprintf('Artificial variables have non-zero values!\n');
        end
        x_opt = x_phase1;
        z_opt = 0;
        feasible = false;
        return;
    end
    
    if verbose
        fprintf('\nPhase I completed successfully!\n');
        fprintf('Initial BFS found: ');
        fprintf('%.4f ', x_phase1);
        fprintf('\n');
    end
    
    % Get final basis from Phase I (excluding artificial variables)
    basis_phase1 = [];
    for i = 1:length(basis)
        if ~ismember(basis(i), artificial_vars)
            basis_phase1 = [basis_phase1, basis(i)];
        end
    end
    
    % If basis is incomplete, add non-artificial variables
    if length(basis_phase1) < m
        for j = 1:n
            if ~ismember(j, artificial_vars) && ~ismember(j, basis_phase1)
                basis_phase1 = [basis_phase1, j];
                if length(basis_phase1) == m
                    break;
                end
            end
        end
    end
    
    %% PHASE II
    if verbose
        fprintf('\n========================================\n');
        fprintf('PHASE II: Optimizing Original Objective\n');
        fprintf('========================================\n\n');
    end
    
    % Remove artificial variable columns from A
    n_original = n - length(artificial_vars);
    A_phase2 = A(:, 1:n_original);
    c_phase2 = c_original(1:n_original);
    
    % Update basis indices (remove artificial variables)
    basis_phase2 = basis_phase1(basis_phase1 <= n_original);
    
    % If basis is incomplete, complete it with slack/surplus variables
    while length(basis_phase2) < m
        for j = 1:n_original
            if ~ismember(j, basis_phase2)
                B_test = A_phase2(:, [basis_phase2, j]);
                if rank(B_test) == min(size(B_test))
                    basis_phase2 = [basis_phase2, j];
                    break;
                end
            end
        end
    end
    
    % Solve Phase II
    [x_phase2, z_phase2, ~] = simplex_solve(A_phase2, b, c_phase2, ...
                                             basis_phase2, verbose);
    
    % Extend solution to include artificial variables (set to 0)
    x_opt = zeros(n, 1);
    x_opt(1:n_original) = x_phase2;
    z_opt = z_phase2;
    feasible = true;
end

function [x_opt, z_opt, iterations] = simplex_solve(A, b, c, basis, verbose)
    % Simplex solver for both phases
    
    [m, n] = size(A);
    iteration = 0;
    max_iterations = 100;
    
    while iteration < max_iterations
        iteration = iteration + 1;
        
        B = A(:, basis);
        
        % Check if basis is singular
        if abs(det(B)) < 1e-10
            % Try to find a better basis
            for j = 1:n
                if ~ismember(j, basis)
                    test_basis = basis;
                    test_basis(1) = j;
                    B_test = A(:, test_basis);
                    if abs(det(B_test)) > 1e-10
                        basis = test_basis;
                        B = B_test;
                        break;
                    end
                end
            end
        end
        
        B_inv = inv(B);
        XB = B_inv * b;
        
        X = zeros(n, 1);
        X(basis) = XB;
        
        CB = c(basis);
        z = CB' * XB;
        
        if verbose && mod(iteration, 5) == 1
            fprintf('Iteration %d: z = %.4e\n', iteration, z);
        end
        
        % Calculate Zj - cj
        Zj_minus_cj = zeros(n, 1);
        for j = 1:n
            Aj = A(:, j);
            Zj_minus_cj(j) = CB' * (B_inv * Aj) - c(j);
        end
        
        [min_val, k] = min(Zj_minus_cj);
        
        if min_val >= -1e-10
            x_opt = X;
            z_opt = z;
            iterations = iteration;
            return;
        end
        
        alpha_k = B_inv * A(:, k);
        
        if all(alpha_k <= 1e-10)
            warning('Problem is unbounded');
            x_opt = X;
            z_opt = inf;
            iterations = iteration;
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
    iterations = iteration;
end
