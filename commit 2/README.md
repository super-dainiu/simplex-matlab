## Commit 2

The following programs aim to solve really simple standard form linear optimization problems, with dual simplex method.
$$
\begin{aligned}
&\min \quad c'x\\
&\textrm{s.t.}\quad  Ax = b\\
&\quad\quad\quad x\ge 0
\end{aligned}
$$


- [linprog_revised.m](linprog_revised.m): Revised simplex method. (Same as commit 1)

  Usage: [x, fval, B_ind, exitflag] = linprog_revised(c, A, b, mute)

  x: optimal solution

  fval: optimal value

  B_ind: optimal basis

  exitflag: 0 for *Optimal solution found*, -2 for *Linprog stopped because no point satisfies the constraints*.

  mute(=false): whether to mute the warnings

- [linprog_dual_revised.m](linprog_dual_revised.m): Revised dual simplex method.  (You need to input a dual feasible solution, and it will converge to primal feasibility.)

​		Usage: [x, fval, exitflag] = linprog_dual_revised(c, A, b, B_ind, mute)

​		x: optimal solution

​		fval: optimal value

​		exitflag: 0 for *Optimal solution found*, -2 for *Linprog stopped because no point satisfies the constraints*.

​		mute(=false): whether to mute the warnings



