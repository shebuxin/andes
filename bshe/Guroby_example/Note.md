# Q&A

**class: optimization method**

Q1: the functinality of GRB

A1: 调用条件/参数状态   check more on gurobi documentation

---

Q2: gurobi multi dict 和python自己的dict有啥区别？

python dict 不能批量初始化

可以从pandas或者python的dict转过来吗

A2: 经过测试： index 可以是string，但不能是未定义的变量

如果idex是数字，可以直接调用 sum(demand)；否则需要 用 sum(demand[idx] for idx in demand)

对于批量定义dict, 例如：note: python的dict

```
p = {i:2z[i,1]+z[i,2]+0.5z[i,3] for i in week} #price in unit
```

---

Q3: z那里是一种什么定义？

A3: 多维变量，在week 10维上扩展成了 10*4，然后idx，变成(i, j)，有点像matlab

用list来add var，生产的dict var

call variable: z[i,1]     i 可以是string

---

Q4: x.prod(p) 是啥意思？prod函数

A4: product 内积

---

Q5: 展示gurobi求解结果相关

m. status 一般有哪几种

GRB.UNBOUNDED

GRB.OPTIMAL

computeIIS( )

m. RC: reduced cost

m. SAObjUp: the range of obj, 这个范围内一直这个约束bounded

m. SAObLow

m. getRow( ).getValue( ) :   代入变量计算这个约束的值

c in m.getConstrs( )          c.Pi, c.RHS, c. SARHSUp, c.SARHSLow: **Pi shadow price**; RHS right hand side; SARHSUp/Low  range of b

---

Q6: addMvar (6,) 和 week的效果相同？然后可以用(6, 4), 就是6*4的变量？

A6: 新建矩阵变量

name.[n, m]

---

Q7: Product@x 和 A@x>b的含义

@是表示矩阵乘法吗？对数据格式有什么要求？A可以是含变量的list吗？

A7: 必须是array，声明变量的时候用addMvar，gurobi dict不能这么写

---

Q8: %s %g and %d

A8: used for formatting strings.

%s acts a placeholder for a string

%d acts as a placeholder for a number.

Refer to the below code for reference:

> name ='giacomo'
>
> number = 4.3
>
> print('%s %s %d %f %g' % (name, number, number, number, number))

The output:

> giacomo 4.3 4 4.300000 4.3
