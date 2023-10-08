# 20230920

V6 ams case is good, but the RoCoF is small when M/D is set as low value

Hence, scale load to 3.5 times based on V6 file, and name it as V7

# 20230919

#### Note

**Base_DC_OPF (123 system DCOPF) :**

    engenvalue calculation after togeler (SG 2 trip) is not correct right after the new calculation

    for TDS, let it happen at 2 seconds to show the stable and sastifed constraints before and after the gen trip

**Test dynamics of 123 system based on DCOPF results :**

load change resutls is different from gen trip results

set all MD to 1 will result in stability issue after gen trip but not for load change

# 20230914 Note

After transferring to gurobi, the NN linearization could be solved.

The problem is the NN prediction accuracy. Both for RoCof, nadir, and eigenvalue

    M/D will be solved to lower boundary, then check the NN prediction, they satisfy the constriants

#### TODO

So need to 1) **generate more typical unstable scenarios**, and also 2) **check the RoCof at lower boundary (check generated data)**

whether to include pg as data input

    generate two group of data: one with pg; and the other without pg

Also, remember to change the upper boundary of VSG, and lower boundary of slack bus, when the NN linearization constraints take effect, pg will be solved to boundary too.

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
