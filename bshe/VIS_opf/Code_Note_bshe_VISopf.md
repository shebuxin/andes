## Probelm formulation (without analytical P_vsg)

### Variable in objective function

$pg$: power of generator

$pru$: up researve

$prd$: down researve

$H_{vsg}$: virtual inertia

$D_{vsg}$: virutal damping      ***Use andes to plot Eig***

### Some intern variable

$$
M_{sys} = \frac {\sum (M_gP_{bg} + M_{vsg}P_{bvsg})}
{\sum P_{bg}+ \sum P_{bvsg}} \quad

D_{sys} = \frac {\sum (D_gP_{bg} + D_{vsg}P_{bvsg})}
{\sum P_{bg}+ \sum P_{ bvsg}}
$$

$$
R_{gsys} = \sum \frac {K_g P_{bg}}{R_g \sum P_{bg}} \quad

F_{gsys} = \sum \frac {K_g P_{bg} F_g}{R_g \sum P_{bg}} \quad

F_g = \frac {pg}{P_{bg}}
$$

Note:
$P_{bg}$ and $P_{bvsg}$ is stored in 'PV' andes file;

$M_g$ and $D_g$ is stored in 'GENROU' andes file;

$M_{vsg}$ and $D_{vsg}$ is stored in 'REGCV1' or 'REGCV2' file;

$R_g$ is soted in 'TGOV1N' andes file;

# Coding Note

### **Andes_var (from andes, padas dataframe):**

1) **gen**: ['idx', 'u', 'name', 'Sn', 'Vn', 'bus', 'p0', 'pmax', 'pmin', 'v0', 'ctrl', 'ramp5', 'ramp10', 'ramp30' ]
   function 'build_group_table (ssa, 'StaticGen', ...) uniforms PV and slack generators and forms group 'StaticGen'

   'ctrl==1' means the gen is controllable; 'ctrl==0' means the gen is uncontrollable;
2) **bus**: ['idx', 'u', 'name', 'Vn', 'vmax', 'vmin', 'v0', 'a0', 'area', 'zone', 'owner']
3) **load**: ['idx', 'u', 'name', 'bus', 'Vn', 'p0', 'q0', 'vmax', 'vmin', 'owner']
4) **line**: ['idx', 'u', 'name', 'bus1', 'bus2', 'Sn', 'fn', 'Vn1', 'Vn2', 'trans', 'tap', 'phi', 'rate_a', 'rate_b', 'rate_c']       *normalized in class function 'from andes'*
5) **gsf:** two dimentional GSF matrix
   row: transfer coefficients of each bus to a line
   column: bus power transfer coefficients to each line
6) **cost**: ['idx', 'c2', 'c1', 'c0', 'cr']

   **---TO DO---**
7) add two column to **gen** [type, M_vsg, D_vsg]

the above dataframe is stored in a dictionary with new_name = var_name+'dict', i.g., gen → gendict

### **Coding_var (python dict var, to fomulate opf):**

#### **gorubi_var**

1. **pg**: **output of generators**
   call a specific generator output: **pg [ gendict.keys ]**, the keys are 'idex', the index of **gen** (generators)

#### constant py dict var

1. **costdict: cost of generators**
   call a specific generator cost: **condict [ gendict.keys ]**, the keys are 'idex' of **gen** (generators)

generator output:

| Constraints                                     | idx | formulation_var  | Coding_var                                                                                | Check |
| ----------------------------------------------- | --- | ---------------- | ----------------------------------------------------------------------------------------- | ----- |
| objective function                              |     | $G_i$, $c_i$ | pg [ gendict.keys ]<br />costdict [ gendict.keys ] ['c0/c1/c2']                           |       |
| Power balance                                   | 01  |                  | pg [ gendict.keys ]<br />ptotal = self.load.p0.sum( )                                     |       |
| Line limit                                      | 02  |                  | linedict[ idx ] [ 'sup' ]  %% 'sup' is load contribution to line, define it in from andes |       |
| SG limit (type I)                               |     |                  | defined in '_build_var'                                                                   |       |
| SG ramping rate                                 |     |                  |                                                                                           |       |
| VSG generation limit (type II)                  |     |                  | defined in '_build_var'                                                                   |       |
| VSG power reserve: (upreserve and down reserve) |     |                  |                                                                                           |       |
| RoCof                                           |     |                  |                                                                                           |       |
| Nadir                                           |     |                  |                                                                                           |       |

RoCo constraint:
$$
RoCof= \Delta P_e / M_{sys}
$$

$$
M_{sys} * RoCof_{lim} \ge \Delta P_e 
$$
frequency nadir constraint:
$$
- f_{nadir-lim} \le f_{nadir-pred} \le f_{nadir-lim}
$$

---

## Q&A

**Q1:** What is $D_t$ in 'TGOV1N' andes?

$D_t$ is damping of governor;

GGEROW also has a damping parameter, which is used in frequency response model.

**Q2:** Where to find $K_g$?

Andes has no $K_g$, assume to be 1

$R_g$ maps power to frequency, droop coefficient

$K_g/R_g$ maps frequecy to power

**Q3:** How to merge parameters in different andes files to a single table?

staticGen = PV + Slack

use function merge

```
res = pd.merge(left=stg.rename(columns={'idx':'gen'}), right=vsg, on='gen', how='left')
```

on shows key name to merge
rename the header with same key

**Check pandas "merge" function documentation:**

https://realpython.com/pandas-merge-join-and-concat/

**Fill pandas nan with " "**

dataframe.fillna('empty', inplace=true)   %把merge表格时产生的nan用给定字符串填充，这里用了‘enpty'，防止andes调用时候报错

inplace = false, 则原表不变，返回一个修改过的表

**Q4:** Should add intern variables as indepent variable?

可以增加中间变量，但是中间变量不需要声明成gurobipy的变量

What if there variable are not included in constriants? i.e. the neural output of MLP, and some binary variables

需要声明成gurobipy变量
