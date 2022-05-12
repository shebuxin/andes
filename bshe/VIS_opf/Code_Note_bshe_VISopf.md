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
