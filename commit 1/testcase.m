%% Test I (Degenerate)
A =  [1 0 2 1;...
      0 1 2 0;...
      -1 0 0 1];
b =  [3;2;-3];
c = [2;4;6;-5];

fprintf('Problem 1: Degenerate\n')
[x, fval] = linprog_revised(c, A, b)
[x, fval] = linprog_tableau(c, A, b)
[x, fval] = linprog(c, [], [], A, b, zeros(numel(c), 1))

%% Test II (Normal)
A =  [0 2 1 -1;...
      2 1 3 0;...
      1 0 0 1];
b =  [3;8;2];
c = [-2;4;-2;3];

fprintf('Problem 2: Normal\n')
[x, fval] = linprog_revised(c, A, b)
[x, fval] = linprog_tableau(c, A, b)
[x, fval] = linprog(c, [], [], A, b, zeros(numel(c), 1))

%% Test III (Infeasible)
A =  [1 3 2 1;...
      0 1 2 0;...
      -1 0 2 1];
b =  [3;2;-3];
c = [2;3;1;2];

fprintf('Problem 3: Infeasible\n')
[x, fval] = linprog_revised(c, A, b)
[x, fval] = linprog_tableau(c, A, b)
[x, fval] = linprog(c, [], [], A, b, zeros(numel(c), 1))