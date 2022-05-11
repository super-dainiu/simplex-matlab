%% Test I
A =  [0 2 1 -1;...
      2 1 3 0;...
      1 0 0 1];
b =  [3;8;2];
c = [-2;4;-2;3];

fprintf('Solve the problem I with revised simplex method.\n')
[x, fval, B_ind] = linprog_revised(c, A, b)
[x, fval] = linprog_dual_revised(c, A, b, B_ind)   % this is trivially correct

% substitute b and convert the primal problem to an infeasible one
fprintf('Update b, which makes problem infeasible\n')
b = [-5;1;-2];
[x, fval] = linprog_dual_revised(c, A, b, B_ind)
[x, fval] = linprog(c, [], [], A, b, zeros(numel(c), 1))

% substitute b and convert the primal problem to a feasible one
fprintf('Update b, which makes problem feasible\n')
b =  [2;4;2];
[x, fval] = linprog_dual_revised(c, A, b, B_ind)
[x, fval] = linprog(c, [], [], A, b, zeros(numel(c), 1))