"""
@author: zguo19
"""

import gurobipy as gp
from gurobipy import GRB
import sys

# proj3 Decide how much to order each week to minimize the total cost

week, demand, minInventory, holdCostUnit = gp.multidict({
    1 : [10, 1, 1],
    2 : [10, 1, 1],
    3 : [10, 1, 1],
    4 : [0, 0, 1],
    5 : [0, 0, 1],
    6 : [15, 1.5,1],
    7 : [20, 2, 1],
    8 : [20, 2, 1],
    9 : [0, 0, 1],
    10 :[10, 1, 1]})
# Model
m = gp.Model("IE522_zguo19_proj3")

# Create decision variables
x = m.addVars(week,lb= 0,  name="x")  #order amount each week
y = m.addVars(week,lb= 0,  name="y")  #Inventory held at end of week
z = m.addVars(week,4,vtype=GRB.BINARY, name="z")   #auxiliary variable

M= sum(demand)
p={i:2*z[i,1]+z[i,2]+0.5*z[i,3] for i in week} #price in unit
# The objective is to minimize the costs
m.setObjective(x.prod(p)+y.prod(holdCostUnit), GRB.MINIMIZE)
                          
# Add constraints
m.addConstrs(y[i]>=minInventory[i] for i in week)
m.addConstr(y[1]==20+x[1]-demand[1])
m.addConstrs(y[i+1]==y[i]+x[i+1]-demand[i+1] for i in range(1,len(week)))
m.addConstrs(z.sum(i,'*')==1 for i in week)
m.addConstrs(x[i]<=9.99*z[i,1]+14.99*z[i,2]+M*z[i,3] for i in week)
m.addConstrs(x[i]>=5*z[i,1]+10*z[i,2]+15*z[i,3] for i in week)

             
m.optimize()
status=m.status
if status == GRB.UNBOUNDED:
    print('The model cannot be solved because it is unbounded')
    sys.exit()
if status == GRB.OPTIMAL:
    print('\nCost: %g' % m.ObjVal)
    print('\norder:')
    for v in m.getVars():
        print('%s %g' % (v.varName, v.x))  
    sys.exit()
if status !=GRB.INF_OR_UNBD and status != GRB.INFEASIBLE:
    print('Optimization was stopped with status %d'  % status)
    sys.exit()

# do IIS
print('The model is infeasible; computing IIS')
m.computeIIS()
if m.IISMinimal:
    print('IIS is minimal\n')
else:
    print('IIS is not minimal\n')
print('\nThe following constraint(s) cannot be satisfied:')
for c in m.getConstrs():
    if c.IISConstr:
        print('%s' % c.ConstrName)
