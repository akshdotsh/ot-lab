% =========================================================================
% Thapar Institute of Engineering and Technology, Patiala
% School of Mathematics
% Optimization Methods (UMA-035)
% Lab Experiment - 1: Graphical Method for Linear Programming Problems
% =========================================================================

function Lab1_Graphical_Method()
    clc;
    clear;
    
    fprintf('========================================\n');
    fprintf('Graphical Method for LPP\n');
    fprintf('========================================\n\n');
    
    % Problem Selection
    fprintf('Select Problem:\n');
    fprintf('1. Max/Min (3x1 + 2x2) s.t. 2x1+4x2<=8, 3x1+5x2>=15, x1,x2>=0\n');
    fprintf('2. Max/Min (3x1 + 2x2) s.t. 2x1+4x2>=8, 3x1+5x2>=15, x1,x2>=0\n');
    fprintf('3. Max/Min (3x1 + 2x2) s.t. 2x1+4x2<=8, 3x1+5x2<=15, x1,x2>=0\n');
    problem = input('Enter problem number (1/2/3): ');
    
    % Objective: Maximize or Minimize
    obj_type = input('Enter 1 for Maximize, 0 for Minimize: ');
    
    % Define problem parameters based on selection
    switch problem
        case 1
            C = [3; 2];  % Objective coefficients
            A = [2, 4; -3, -5];  % Constraint coefficients (converted to <= form)
            B = [8; -15];  % RHS (converted to <= form)
            constraint_signs = {'<=', '>='};
        case 2
            C = [3; 2];
            A = [-2, -4; -3, -5];
            B = [-8; -15];
            constraint_signs = {'>=', '>='};
        case 3
            C = [3; 2];
            A = [2, 4; 3, 5];
            B = [8; 15];
            constraint_signs = {'<=', '<='};
        otherwise
            error('Invalid problem selection');
    end
    
    % Call graphical method solver
    [x_opt, z_opt] = graphical_method(C, A, B, obj_type);
    
    % Display results
    fprintf('\n========================================\n');
    fprintf('RESULTS\n');
    fprintf('========================================\n');
    if ~isempty(x_opt)
        fprintf('Optimal Solution: x1 = %.4f, x2 = %.4f\n', x_opt(1), x_opt(2));
        fprintf('Optimal Objective Value: z = %.4f\n', z_opt);
    else
        fprintf('Problem is INFEASIBLE or UNBOUNDED\n');
    end
end

function [x_opt, z_opt] = graphical_method(C, A, B, obj_type)
    % Step 1: Data entered in arrays C, A, B
    
    % Step 2: Select range of x1 for plotting
    x1_range = linspace(0, 10, 1000);
    
    % Step 3 & 4: Plot constraints
    figure('Position', [100, 100, 800, 600]);
    hold on;
    grid on;
    xlabel('x_1', 'FontSize', 12);
    ylabel('x_2', 'FontSize', 12);
    title('Graphical Method for LPP', 'FontSize', 14);
    
    % Plot each constraint
    [m, ~] = size(A);
    colors = ['r', 'b', 'g', 'm', 'c'];
    
    for i = 1:m
        if A(i,2) ~= 0
            % x2 = (B(i) - A(i,1)*x1) / A(i,2)
            x2_i = (B(i) - A(i,1)*x1_range) / A(i,2);
            % Keep only non-negative values
            x2_i(x2_i < 0) = 0;
            plot(x1_range, x2_i, colors(i), 'LineWidth', 2, ...
                 'DisplayName', sprintf('Constraint %d', i));
        else
            % Vertical line: x1 = B(i)/A(i,1)
            x1_val = B(i)/A(i,1);
            plot([x1_val x1_val], [0 10], colors(i), 'LineWidth', 2, ...
                 'DisplayName', sprintf('Constraint %d', i));
        end
    end
    
    % Step 5-12: Find corner points
    solution = [];
    
    % Intersections of constraints
    for i = 1:m
        for j = i+1:m
            % Step 6-10: Extract rows and combine
            A1 = A(i,:);
            B1 = B(i);
            A2 = A(j,:);
            B2 = B(j);
            
            A3 = [A1; A2];
            B3 = [B1; B2];
            
            % Step 11: Solve system A3*X = B3
            if rank(A3) == 2
                X = A3 \ B3;
                solution = [solution; X'];
            end
        end
    end
    
    % Intersections with axes
    for i = 1:m
        % Intersection with x1-axis (x2=0)
        if A(i,1) ~= 0
            x1_val = B(i)/A(i,1);
            solution = [solution; x1_val, 0];
        end
        % Intersection with x2-axis (x1=0)
        if A(i,2) ~= 0
            x2_val = B(i)/A(i,2);
            solution = [solution; 0, x2_val];
        end
    end
    
    % Add origin
    solution = [solution; 0, 0];
    
    % Remove duplicates
    solution = unique(solution, 'rows');
    
    % Step 13-20: Filter feasible solutions
    feasible_solutions = [];
    
    for k = 1:size(solution,1)
        x1 = solution(k,1);
        x2 = solution(k,2);
        
        % Check non-negativity
        if x1 >= -1e-6 && x2 >= -1e-6
            % Check all constraints
            feasible = true;
            for i = 1:m
                if A(i,:)*[x1; x2] > B(i) + 1e-6
                    feasible = false;
                    break;
                end
            end
            
            if feasible
                feasible_solutions = [feasible_solutions; x1, x2];
            end
        end
    end
    
    % Plot feasible region
    if ~isempty(feasible_solutions)
        % Sort points for polygon
        if size(feasible_solutions,1) > 2
            k = convhull(feasible_solutions(:,1), feasible_solutions(:,2));
            fill(feasible_solutions(k,1), feasible_solutions(k,2), ...
                 'y', 'FaceAlpha', 0.3, 'DisplayName', 'Feasible Region');
        end
        
        % Plot corner points
        plot(feasible_solutions(:,1), feasible_solutions(:,2), ...
             'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'r', ...
             'DisplayName', 'Corner Points');
    end
    
    legend('Location', 'best');
    xlim([0, 10]);
    ylim([0, 10]);
    
    % Step 21-24: Find optimal solution
    if isempty(feasible_solutions)
        fprintf('No feasible solution found!\n');
        x_opt = [];
        z_opt = [];
        return;
    end
    
    OBJ = feasible_solutions * C;
    
    if obj_type == 1  % Maximize
        [z_opt, idx] = max(OBJ);
    else  % Minimize
        [z_opt, idx] = min(OBJ);
    end
    
    x_opt = feasible_solutions(idx, :)';
    
    % Highlight optimal point
    plot(x_opt(1), x_opt(2), 'p', 'MarkerSize', 20, ...
         'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k', ...
         'LineWidth', 2, 'DisplayName', 'Optimal Solution');
    
    % Plot objective function contours
    [X1, X2] = meshgrid(0:0.1:10, 0:0.1:10);
    Z = C(1)*X1 + C(2)*X2;
    contour(X1, X2, Z, 20, '--', 'LineWidth', 0.5);
    
    legend('Location', 'best');
    hold off;
end
