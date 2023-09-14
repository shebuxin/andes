# Q&A

1. general architecture of AMS
2. check the detail of LDCOPF
   then find the source file (didn't find it)
3. let Jinning know the new model, and compare the efforts of hard coding and general coding
   let him know the way to modify the cvxpy problem
4. vector formulation in cvxpy
   (check the source code for problem formulation)

# Problem 20230904

1. cannot define lb and ub of cvxpy variables
   **cvx doesnot support this functionality**
2. cannot call value with ssp.LDOPF.om.mdl, can call with ssp.LDOPF.om
   **we can extract value of the decision variable using ssp.LDOPF.om (Note, the original cvx decision variable) after the new problem is solved
   Cannot call it suing ssp.LDOPF, because it is the routine variables values of the original model**
3. check the dimension of x and np.array([1, 2]) and x @ np.array([1, 2])
   **use * for coefficient multiply, and use @ for matrix multiply (incluiding 1dim matrix multipy 1 dim matrix)**
4. How to prepare a new case
   and load it using ams
   **keep topology list and add Gcost data**

# Problem 20230906

1. for main_vis_v0: the dispatch results are different from those after modification
   the results are wired, PV_10 and PV_11 are 0
   check the formulated obj
   check the coefficients (may be not matching the idx)
2. cannot call gurobi as solver
   what is the default solver?

# Tips

check new obj formulation and compare it with the original ams LDOPF

# Problem summary

**Problem summary**

1. variable idx issue: c2 and c1 doesn't match gen variable if the total number of gen is large than 10
2. cannot call gurobi as a solver due to c++ issue

**Some expected features**

* be able to freely add new decision variables, new constraints, and modify the existing constraints
* be able to define modify objective function, or defined customized function based on existing obj
* be able to include new parameters and easily call them




# NN optimization debug note

hup and hdown will imapct the optimization results

hup = 1000 (hdown = -1000) will result in infeasibility

Too large inertia cost coefficient will also results in infeasibility
