function [x, fval, exitflag] = linprog_tableau(c, A, b, mute)
    if ~exist('mute','var')
      mute = false;
    end
    [m, n] = size(A);
    exitflag = 0;
    eps = 1e-16;
    tol = 1e-12;

    %% Phase 1
    % Save original problem
    ORG = true(n, 1);

    % Construct LP_aux
    b_aux = abs(b);
    A_aux = A; A_aux(:, b < 0) = - A_aux(:, b < 0); A_aux = [A_aux, eye(m)];
    c_aux = [zeros(n, 1); ones(m, 1)];

    % Construct tableau
    B = [false(n, 1); true(m, 1)];
    B_ind = find(B);
    BA = A_aux(:, B) \ A_aux;
    Bb = A_aux(:, B) \ b_aux;
    c_bar = c_aux' - c_aux(B)' * BA;
    fval = - c_aux(B)' * Bb;

    % Start iteration
    iter = 0;
    while sum(c_bar < -tol)
        ROWS = true(1, m);

        % Basis out
        [~, pivot_c] = find(c_bar < -tol, 1, 'first');

        % Stepsize
        if sum((BA(:, pivot_c) < 0) | (abs(Bb .* BA(:, pivot_c)) < tol))  == m
            break
        end
        theta = Bb ./ BA(:, pivot_c);
        theta(BA(:, pivot_c) <= tol) = inf;
        [~, pivot_r] = min(theta);

        % Basis in
        ROWS(pivot_r) = ~ROWS(pivot_r);
        B([pivot_c, B_ind(pivot_r)]) = ~B([pivot_c, B_ind(pivot_r)]);
        B_ind(pivot_r) = pivot_c;

        % Swapping
        Bb(pivot_r) = Bb(pivot_r) / BA(pivot_r, pivot_c);
        Bb(ROWS) = Bb(ROWS) - BA(ROWS, pivot_c) * Bb(pivot_r);
        BA(pivot_r, :) = BA(pivot_r, :) / BA(pivot_r, pivot_c);
        BA(ROWS, :) = BA(ROWS, :) - BA(ROWS, pivot_c) * BA(pivot_r, :);
        fval = fval - c_bar(pivot_c) * Bb(pivot_r);
        c_bar = c_bar - c_bar(pivot_c) * BA(pivot_r, :);
        BA(abs(BA) < tol) = 0; Bb(Bb < 0) = 0; c(abs(c) < eps) = 0;
        
        iter = iter + 1;
    end
    
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

    % Construct tableau
    B = B(ORG);
    B_ind = find(B);
    BA = A(:, B) \ A;
    Bb = A(:, B) \ b;
    c_bar = c' - c(B)' * BA;
    fval = - c(B)' * Bb;

    % Start iteration
    while sum(c_bar < 0)
        ROWS = true(1, m);

        % Basis out
        [~, pivot_c] = find(c_bar < -tol, 1, 'first');

        % Stepsize
        if sum(BA(:, pivot_c) < 0) == m
            exitflag = -3;  % Unbounded
            x = zeros(n, 1);
            x(B_ind) = Bb;
            fval = -inf;
            fprintf('The problem is unbounded.\n');
            return
        end
        theta = Bb ./ BA(:, pivot_c);
        theta(BA(:, pivot_c) < 0) = inf;
        [~, pivot_r] = min(theta);

        % Basis in
        ROWS(pivot_r) = ~ROWS(pivot_r);
        B([pivot_c, B_ind(pivot_r)]) = ~B([pivot_c, B_ind(pivot_r)]);
        B_ind(pivot_r) = pivot_c;

        % Swapping
        Bb(pivot_r) = Bb(pivot_r) / BA(pivot_r, pivot_c);
        BA(pivot_r, :) = BA(pivot_r, :) / BA(pivot_r, pivot_c);
        Bb(ROWS) = Bb(ROWS) - BA(ROWS, pivot_c) * Bb(pivot_r);
        BA(ROWS, :) = BA(ROWS, :) - BA(ROWS, pivot_c) * BA(pivot_r, :);
        fval = fval - c_bar(pivot_c) * Bb(pivot_r);
        c_bar = c_bar - c_bar(pivot_c) * BA(pivot_r, :);
    end
    x = zeros(n, 1);
    x(B_ind) = Bb;
    fval = -fval;
    if ~mute
        fprintf('Optimal solution found.\n');
    end
end