## Commit 1

The following programs aim to solve really simple standard form linear optimization problems.  

$$
\begin{aligned}
&\min \quad c'x\\
&\textrm{s.t.}\quad  Ax = b\\
&\quad\quad\quad x\ge 0
\end{aligned}
$$


- [linprog_revised.m](linprog_revised.m): Revised simplex method. 

  Usage: [x, fval, exitflag] = linprog_revised(c, A, b, mute)

  x: optimal solution

  fval: optimal value

  exitflag: 0 for *Optimal solution found*, -2 for *Linprog stopped because no point satisfies the constraints*.

  mute(=false): whether to mute the warnings

- [linprog_tableau.m](linprog_tableau.m): Tableau simplex method. 

​		Usage: [x, fval, exitflag] = linprog_tableau(c, A, b, mute)

​		x: optimal solution

​		fval: optimal value

​		exitflag: 0 for *Optimal solution found*, -2 for *Linprog stopped because no point satisfies the constraints*.

​		mute(=false): whether to mute the warnings



