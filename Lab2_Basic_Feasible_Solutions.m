% =========================================================================
% Thapar Institute of Engineering and Technology, Patiala
% School of Mathematics
% Optimization Methods (UMA-035)
% Lab Experiment - 2: Basic Feasible Solutions
% =========================================================================

function Lab2_Basic_Feasible_Solutions()
    clc;
    clear;
    
    fprintf('========================================\n');
    fprintf('Basic Feasible Solutions of LPP\n');
    fprintf('========================================\n\n');
    
    fprintf('Select Problem:\n');
    fprintf('1. Max z = x1 + 2x2, s.t. -x1+x2<=1, x1+x2<=2, x1,x2>=0\n');
    fprintf('2. Check basis (x1,x4) and (x3,x2)\n');
    fprintf('3. Max z = -x1+2x2-x3, s.t. x1<=4, x2<=4, -x1+x2<=6, -x1+2x3<=4\n');
    fprintf('4. Check for degenerate BFS\n');
    problem = input('Enter problem number (1-4): ');
    
    switch problem
        case 1
            problem1_all_bfs();
        case 2
            problem2_check_basis();
        case 3
            problem3_solve_by_bfs();
        case 4
            problem4_degenerate_bfs();
        otherwise
            error('Invalid problem selection');
    end
end

function problem1_all_bfs()
    % Convert to standard form: Max z = x1 + 2x2
    % -x1 + x2 + s1 = 1
    %  x1 + x2 + s2 = 2
    % x1, x2, s1, s2 >= 0
    
    fprintf('\nProblem 1: Finding all BFS\n');
    fprintf('---------------------------\n');
    
    C = [1; 2; 0; 0];
    A = [-1, 1, 1, 0;
          1, 1, 0, 1];
    b = [1; 2];
    
    [m, n] = size(A);
    fprintf('Number of constraints (m): %d\n', m);
    fprintf('Number of variables (n): %d\n', n);
    fprintf('Number of possible bases: C(%d,%d) = %d\n\n', n, m, nchoosek(n,m));
    
    % Find all combinations of m variables
    all_bases = nchoosek(1:n, m);
    num_bases = size(all_bases, 1);
    
    fprintf('Checking all possible bases:\n');
    fprintf('============================\n\n');
    
    bfs_count = 0;
    all_bfs = [];
    
    for i = 1:num_bases
        basis_indices = all_bases(i, :);
        fprintf('Basis %d: Variables {', i);
        fprintf('x%d ', basis_indices);
        fprintf('}\n');
        
        [is_bfs, X_full, type] = check_bfs(A, b, basis_indices, n);
        
        if is_bfs
            bfs_count = bfs_count + 1;
            all_bfs = [all_bfs; X_full'];
            z_value = C' * X_full;
            fprintf('  ✓ %s\n', type);
            fprintf('  Solution: ');
            for j = 1:n
                fprintf('x%d=%.4f ', j, X_full(j));
            end
            fprintf('\n  Objective value: z = %.4f\n\n', z_value);
        else
            fprintf('  ✗ %s\n\n', type);
        end
    end
    
    fprintf('============================\n');
    fprintf('Total BFS found: %d out of %d bases\n', bfs_count, num_bases);
    
    if bfs_count > 0
        obj_values = all_bfs * C;
        [max_val, max_idx] = max(obj_values);
        fprintf('\nOptimal Solution (Maximum):\n');
        fprintf('x1 = %.4f, x2 = %.4f\n', all_bfs(max_idx,1), all_bfs(max_idx,2));
        fprintf('Optimal Value: z = %.4f\n', max_val);
    end
end

function problem2_check_basis()
    % Max z = x1 + 2x2 - x3 + x4
    % x1 + x2 - x3 + 3x4 = 15
    % 5x1 + x2 + 4x3 + 15x4 = 12
    
    fprintf('\nProblem 2: Checking specific bases\n');
    fprintf('-----------------------------------\n');
    
    A = [1, 1, -1, 3;
         5, 1,  4, 15];
    b = [15; 12];
    
    % Check basis (x1, x4)
    fprintf('\n(1) Checking basis (x1, x4):\n');
    basis1 = [1, 4];
    [is_bfs1, X1, type1] = check_bfs(A, b, basis1, 4);
    fprintf('Result: %s\n', type1);
    if is_bfs1
        fprintf('Solution: x1=%.4f, x2=%.4f, x3=%.4f, x4=%.4f\n', X1);
    end
    
    % Check basis (x3, x2)
    fprintf('\n(2) Checking basis (x3, x2):\n');
    basis2 = [3, 2];
    [is_bfs2, X2, type2] = check_bfs(A, b, basis2, 4);
    fprintf('Result: %s\n', type2);
    if is_bfs2
        fprintf('Solution: x1=%.4f, x2=%.4f, x3=%.4f, x4=%.4f\n', X2);
    end
end

function problem3_solve_by_bfs()
    % Convert to standard form: Max z = -x1 + 2x2 - x3
    % x1 + s1 = 4
    % x2 + s2 = 4
    % -x1 + x2 + s3 = 6
    % -x1 + 2x3 + s4 = 4
    
    fprintf('\nProblem 3: Solve by finding all BFS\n');
    fprintf('------------------------------------\n');
    
    C = [-1; 2; -1; 0; 0; 0; 0];
    A = [1, 0, 0, 1, 0, 0, 0;
         0, 1, 0, 0, 1, 0, 0;
        -1, 1, 0, 0, 0, 1, 0;
        -1, 0, 2, 0, 0, 0, 1];
    b = [4; 4; 6; 4];
    
    [m, n] = size(A);
    all_bases = nchoosek(1:n, m);
    num_bases = size(all_bases, 1);
    
    fprintf('Total possible bases: %d\n\n', num_bases);
    
    bfs_count = 0;
    all_bfs = [];
    
    for i = 1:num_bases
        basis_indices = all_bases(i, :);
        [is_bfs, X_full, ~] = check_bfs(A, b, basis_indices, n);
        
        if is_bfs
            bfs_count = bfs_count + 1;
            all_bfs = [all_bfs; X_full'];
        end
    end
    
    fprintf('Total BFS found: %d\n', bfs_count);
    
    if bfs_count > 0
        obj_values = all_bfs * C;
        [max_val, max_idx] = max(obj_values);
        fprintf('\nOptimal Solution:\n');
        fprintf('x1 = %.4f, x2 = %.4f, x3 = %.4f\n', ...
                all_bfs(max_idx,1), all_bfs(max_idx,2), all_bfs(max_idx,3));
        fprintf('Optimal Value: z = %.4f\n', max_val);
    end
end

function problem4_degenerate_bfs()
    % Max z = x1 + x2 + x3
    % x1 + x2 + s1 = 1
    % -x2 + x3 + s2 = 0
    
    fprintf('\nProblem 4: Checking for degenerate BFS\n');
    fprintf('---------------------------------------\n');
    
    C = [1; 1; 1; 0; 0];
    A = [1, 1, 0, 1, 0;
         0, -1, 1, 0, 1];
    b = [1; 0];
    
    [m, n] = size(A);
    all_bases = nchoosek(1:n, m);
    num_bases = size(all_bases, 1);
    
    degenerate_bfs = [];
    degenerate_bases = {};
    
    for i = 1:num_bases
        basis_indices = all_bases(i, :);
        [is_bfs, X_full, type] = check_bfs(A, b, basis_indices, n);
        
        if is_bfs && contains(type, 'Degenerate')
            degenerate_bfs = [degenerate_bfs; X_full'];
            degenerate_bases{end+1} = basis_indices;
            
            fprintf('Degenerate BFS found with basis: ');
            fprintf('x%d ', basis_indices);
            fprintf('\n  Solution: ');
            for j = 1:n
                fprintf('x%d=%.4f ', j, X_full(j));
            end
            fprintf('\n\n');
        end
    end
    
    if ~isempty(degenerate_bfs)
        fprintf('Total degenerate BFS found: %d\n', size(degenerate_bfs,1));
        fprintf('All these correspond to the same point!\n');
    else
        fprintf('No degenerate BFS found.\n');
    end
end

function [is_bfs, X_full, type_str] = check_bfs(A, b, basis_indices, n)
    % Check if given basis produces a BFS
    % Returns: is_bfs (boolean), X_full (complete solution), type_str (description)
    
    [m, ~] = size(A);
    
    % Step 1: Construct basis matrix B
    B = A(:, basis_indices);
    
    % Step 2: Check if det(B) != 0
    if abs(det(B)) < 1e-10
        is_bfs = false;
        X_full = zeros(n, 1);
        type_str = 'Not a basic solution (det(B) = 0)';
        return;
    end
    
    % Find XB = B^(-1) * b
    XB = B \ b;
    
    % Create full solution vector
    X_full = zeros(n, 1);
    X_full(basis_indices) = XB;
    
    % Step 3: Check feasibility
    tolerance = 1e-6;
    
    if any(XB < -tolerance)
        is_bfs = false;
        type_str = 'Not a BFS (negative components)';
        return;
    end
    
    % Check for degeneracy
    if any(abs(XB) < tolerance)
        is_bfs = true;
        type_str = 'Degenerate BFS';
    else
        is_bfs = true;
        type_str = 'Non-degenerate BFS';
    end
end
