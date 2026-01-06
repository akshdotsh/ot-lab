% =========================================================================
% Thapar Institute of Engineering and Technology, Patiala
% School of Mathematics
% Optimization Methods (UMA-035)
% Lab Experiment - 7: Least Cost Method for Transportation Problem
% =========================================================================

function Lab7_Least_Cost_Method()
    clc;
    clear;
    
    fprintf('========================================\n');
    fprintf('Least Cost Method - Transportation Problem\n');
    fprintf('========================================\n\n');
    
    fprintf('Select Problem:\n');
    fprintf('1. 3 Sources × 4 Destinations\n');
    fprintf('2. 4 Sources × 5 Destinations\n');
    problem = input('Enter problem number (1-2): ');
    
    switch problem
        case 1
            problem1_least_cost();
        case 2
            problem2_least_cost();
        otherwise
            error('Invalid problem selection');
    end
end

function problem1_least_cost()
    % Problem 1: 3 Sources × 4 Destinations
    fprintf('\nProblem 1: Transportation Problem\n');
    fprintf('----------------------------------\n');
    
    % Cost matrix
    C = [2, 10, 4, 5;
         1, 6, 12, 8;
         3, 9, 5, 7];
    
    % Supply (availability)
    a = [12; 25; 20];
    
    % Demand (requirement)
    b = [25; 10; 15; 5];
    
    fprintf('\nCost Matrix:\n');
    disp(C);
    fprintf('Supply: ');
    fprintf('%d ', a);
    fprintf('\n');
    fprintf('Demand: ');
    fprintf('%d ', b);
    fprintf('\n');
    
    % Check if balanced
    total_supply = sum(a);
    total_demand = sum(b);
    fprintf('\nTotal Supply: %d\n', total_supply);
    fprintf('Total Demand: %d\n', total_demand);
    
    if total_supply ~= total_demand
        fprintf('WARNING: Problem is unbalanced!\n');
        if total_supply > total_demand
            % Add dummy destination
            C = [C, zeros(size(C,1), 1)];
            b = [b; total_supply - total_demand];
            fprintf('Added dummy destination with demand %d\n', ...
                    total_supply - total_demand);
        else
            % Add dummy source
            C = [C; zeros(1, size(C,2))];
            a = [a; total_demand - total_supply];
            fprintf('Added dummy source with supply %d\n', ...
                    total_demand - total_supply);
        end
    end
    
    [X, total_cost] = least_cost_method(C, a, b, true);
    
    fprintf('\n========================================\n');
    fprintf('FINAL SOLUTION\n');
    fprintf('========================================\n');
    fprintf('Allocation Matrix:\n');
    disp(X);
    fprintf('Total Transportation Cost: %.2f\n', total_cost);
end

function problem2_least_cost()
    % Problem 2: 4 Sources × 5 Destinations
    fprintf('\nProblem 2: Transportation Problem\n');
    fprintf('----------------------------------\n');
    
    C = [3, 11, 4, 14, 15;
         6, 16, 18, 2, 28;
         10, 13, 15, 19, 17;
         7, 12, 5, 8, 9];
    
    a = [15; 25; 10; 15];
    b = [20; 10; 15; 15; 5];
    
    fprintf('\nCost Matrix:\n');
    disp(C);
    fprintf('Supply: ');
    fprintf('%d ', a);
    fprintf('\n');
    fprintf('Demand: ');
    fprintf('%d ', b);
    fprintf('\n');
    
    total_supply = sum(a);
    total_demand = sum(b);
    fprintf('\nTotal Supply: %d\n', total_supply);
    fprintf('Total Demand: %d\n', total_demand);
    
    if total_supply ~= total_demand
        fprintf('WARNING: Problem is unbalanced!\n');
        if total_supply > total_demand
            C = [C, zeros(size(C,1), 1)];
            b = [b; total_supply - total_demand];
        else
            C = [C; zeros(1, size(C,2))];
            a = [a; total_demand - total_supply];
        end
    end
    
    [X, total_cost] = least_cost_method(C, a, b, true);
    
    fprintf('\n========================================\n');
    fprintf('FINAL SOLUTION\n');
    fprintf('========================================\n');
    fprintf('Allocation Matrix:\n');
    disp(X);
    fprintf('Total Transportation Cost: %.2f\n', total_cost);
end

function [X, total_cost] = least_cost_method(C, a, b, verbose)
    % Least Cost Method for Transportation Problem
    % Inputs:
    %   C: cost matrix (m × n)
    %   a: supply vector (m × 1)
    %   b: demand vector (n × 1)
    %   verbose: display steps (true/false)
    
    [m, n] = size(C);
    X = zeros(m, n);  % Allocation matrix
    
    % Make copies of supply and demand
    supply = a;
    demand = b;
    
    % Create a working cost matrix
    C_work = C;
    
    % Total allocations needed: m + n - 1
    total_allocations = m + n - 1;
    k = 0;  % Allocation counter
    
    if verbose
        fprintf('\n========================================\n');
        fprintf('Least Cost Method - Step by Step\n');
        fprintf('========================================\n\n');
    end
    
    while k < total_allocations
        k = k + 1;
        
        if verbose
            fprintf('Step %d:\n', k);
            fprintf('-------\n');
        end
        
        % Step 1: Find minimum cost cell
        [min_cost, linear_idx] = min(C_work(:));
        [p, q] = ind2sub([m, n], linear_idx);
        
        if verbose
            fprintf('Minimum cost: C(%d,%d) = %.2f\n', p, q, min_cost);
        end
        
        % Step 2: Allocate min(supply(p), demand(q))
        allocation = min(supply(p), demand(q));
        X(p, q) = allocation;
        
        if verbose
            fprintf('Allocation: x(%d,%d) = %.2f\n', p, q, allocation);
            fprintf('  min(supply(%d), demand(%d)) = min(%.2f, %.2f) = %.2f\n', ...
                    p, q, supply(p), demand(q), allocation);
        end
        
        % Step 3: Update supply and demand
        if supply(p) <= demand(q)
            % Supply exhausted first
            demand(q) = demand(q) - supply(p);
            supply(p) = 0;
            
            if verbose
                fprintf('  Supply at source %d exhausted\n', p);
                fprintf('  Remaining demand at destination %d: %.2f\n', q, demand(q));
            end
            
            % Mark entire row as unavailable (set to large number)
            C_work(p, :) = 1e10;
        else
            % Demand satisfied first
            supply(p) = supply(p) - demand(q);
            demand(q) = 0;
            
            if verbose
                fprintf('  Demand at destination %d satisfied\n', q);
                fprintf('  Remaining supply at source %d: %.2f\n', p, supply(p));
            end
            
            % Mark entire column as unavailable
            C_work(:, q) = 1e10;
        end
        
        if verbose
            fprintf('\n');
        end
    end
    
    % Calculate total cost
    total_cost = sum(sum(C .* X));
    
    if verbose
        fprintf('========================================\n');
        fprintf('All allocations completed!\n');
        fprintf('========================================\n\n');
    end
end
