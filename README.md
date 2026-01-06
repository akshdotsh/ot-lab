# List of experiments

https://sites.google.com/thapar.edu/meenakshirana/Current-Semester-2020/optimization-techniques_1

**1. Graphical method for LPP**
Formulate the LPP → plot each constraint as an equation → identify the feasible region → find corner (extreme) points → evaluate the objective function at each corner → the best value gives the optimal solution.

---

**2. Basic solutions & bounded LPP**
Convert constraints to equations → select (m) variables (where (m) = number of equations) and set remaining to zero → solve to get **basic solutions** → feasible basic solutions are corner points → evaluate objective function at each → choose optimal (bounded if maximum/minimum exists).

---

**3. Simplex method (≤ type constraints)**
Convert inequalities to equations using slack variables → write initial simplex table → identify entering variable (most negative (C_j - Z_j)) → identify leaving variable (minimum positive ratio) → pivot → repeat until all (C_j - Z_j \ge 0).

---

**4. Big M method (≥ type constraints)**
Convert constraints using surplus and artificial variables → add large penalty (M) to artificial variables in objective function → apply simplex method → eliminate artificial variables → final table gives optimal solution.

---

**5. Two-phase method**
Phase I: remove artificial variables by minimizing their sum → check feasibility → Phase II: optimize original objective function using feasible solution from Phase I.

---

**6. Dual simplex method**
Used when solution is infeasible but optimality condition holds → choose most negative RHS as leaving variable → select entering variable using minimum ratio → iterate until feasibility and optimality are achieved.

---

**7. Transportation problem (Least Cost Method)**
Identify the cell with minimum cost → allocate maximum possible supply/demand → adjust table → repeat until all supplies and demands are satisfied → gives initial basic feasible solution.

---

**8. Weighted sum method (Multi-objective LPP)**
Assign weights to each objective → convert multiple objectives into a single objective function → solve using standard LPP techniques → solution depends on chosen weights.

---

**9. Fibonacci search technique**
Used for 1D unconstrained optimization → interval is reduced using Fibonacci ratios → evaluate function at two points → eliminate non-optimal interval → repeat until desired accuracy → gives optimal point.
