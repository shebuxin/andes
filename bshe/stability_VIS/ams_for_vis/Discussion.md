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
2. cannot call value with ssp.LDOPF.om.mdl, can call with ssp.LDOPF.om
3. check the dimension of x and np.array([1, 2]) and x @ np.array([1, 2])


# Tips

check new obj formulation and compare it with the original ams LDOPF
