function [x, fval, exitflag] = linprog_dual_revised(c, A, b, init, mute)
    if ~exist('mute','var')
      mute = false;
    end
    [m, n] = size(A);
    exitflag = 0;
    tol = 1e-16;
    
    B_ind = sort(init);
    B_inv = inv(A(:, B_ind)); p = c(B_ind)' * B_inv; x = zeros(n, 1); x(B_ind) = B_inv * b;
    c_bar = c' - p * A;
    
    % Start iteration
    iter = 0;
    while sum(x < -tol)
        ROWS = true(1, m);
        
        % Basis in
        [pivot_r, ~] = find(x < -tol, 1, 'first'); 
        
        % Stepsize
        v = B_inv * A;
        v = v(pivot_r, :);
        if sum(v > -tol) - sum(v(B_ind) > -tol) == n - m
            fprintf('Linprog stopped because no point satisfies the constraints.\n');
            x = [];
            fval = [];
            exitflag = -2;
            return
        end
        theta = c_bar ./ -v;
        theta(v > -tol) = inf; theta(B_ind) = inf;
        [~, pivot_c] = min(theta);
        u = B_inv * A(:, pivot_c);

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