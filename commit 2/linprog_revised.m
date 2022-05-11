function [x, fval, B_ind, exitflag] = linprog_revised(c, A, b, mute)
    if ~exist('mute','var')
      mute = false;
    end
    [m, n] = size(A);
    exitflag = 0;
    tol = 1e-12;

    %% Phase 1
    % Save original problem
    ORG = true(n, 1);

    % Construct LP_aux
    b_aux = abs(b);
    A_aux = A; A_aux(:, b < 0) = - A_aux(:, b < 0); A_aux = [A_aux, eye(m)];
    c_aux = [zeros(n, 1); ones(m, 1)];

    % Construct tableau
    B = [false(n, 1); true(m, 1)]; B_ind = find(B);
    B_inv = inv(A_aux(:, B_ind)); p = c_aux(B_ind)' * B_inv;
    c_bar = c_aux' - p * A_aux;
    x = zeros(m+n, 1);
    x(B_ind) = B_inv * b_aux;
    
    
    % Start iteration
    iter = 0;
    while sum(c_bar < -tol)
        ROWS = true(1, m);
        
        % Basis in
        [~, pivot_c] = find(c_bar < -tol, 1, 'first');
        
        % Stepsize
        u = B_inv * A_aux(:, pivot_c);
        if sum(u < tol) == m
            break
        end
        theta = x(B_ind) ./ u; theta(u < tol) = inf;
        [~, pivot_r] = min(theta);
        

        % Basis out
        B_ind(pivot_r) = pivot_c;
        B_inv(pivot_r, :) = B_inv(pivot_r, :) / u(pivot_r);
        ROWS = true(1, m); ROWS(pivot_r) = false;
        B_inv(ROWS, :) = B_inv(ROWS, :) - B_inv(pivot_r, :) .* u(ROWS);
        
        
        % Update c_bar
        p = c_aux(B_ind)' * B_inv;
        x = zeros(m+n, 1);
        x(B_ind) = B_inv * b_aux;
        c_bar = c_aux' - p * A_aux;
        iter = iter + 1;
    end
    fval = c_aux(B_ind)' * B_inv * b_aux;
    
    % Rule out infeasibility
    if abs(fval - 0) > tol
        fprintf('Linprog stopped because no point satisfies the constraints.\n');
        x = [];
        fval = [];
        exitflag = -2;
        return
    end
    
    % Degeneracy
    if any(B_ind > numel(ORG))
        BA = inv(A(:, B_ind)) * A;
        for row = find(B_ind > numel(ORG))'
            for col = 1:numel(ORG)
                if ~sum(B_ind == col) && (BA(row, col) ~= 0)
                    B([B_ind(row), col]) = ~B([B_ind(row), col]);
                    B_ind(row) = col;
                end
            end
        end
    end
    
    %% Phase 2
    B_inv = inv(A(:, B_ind)); p = c(B_ind)' * B_inv; x(B_ind) = B_inv * b;
    x = zeros(n, 1); x(B_ind) = B_inv * b;
    c_bar = c' - p * A;
    
    % Start iteration
    iter = 0;
    while sum(c_bar < -tol)
        ROWS = true(1, m);
        
        % Basis in
        [~, pivot_c] = find(c_bar < -tol, 1, 'first'); 
        
        % Stepsize
        u = B_inv * A(:, pivot_c);
        if sum(u < tol) == m
            break
        end
        theta = x(B_ind) ./ u; theta(u < tol) = inf;
        [~, pivot_r] = min(theta);
        

        % Basis out
        B_ind(pivot_r) = pivot_c;
        B_inv(pivot_r, :) = B_inv(pivot_r, :) / u(pivot_r);
        ROWS = true(1, m); ROWS(pivot_r) = false;
        B_inv(ROWS, :) = B_inv(ROWS, :) - B_inv(pivot_r, :) .* u(ROWS);
        
        % Update c_bar
        p = c(B_ind)' * B_inv;
        x = zeros(n, 1);
        x(B_ind) = B_inv * b;
        c_bar = c' - p * A;
        iter = iter + 1;
    end
    fprintf("Optimal solution found.\n");
    fval = c(B_ind)' * B_inv * b;
end