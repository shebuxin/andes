## Probelm formulation (without analytical Pvsg)

### variable in objective function

$pg$: power of generator

$pru$: up researve

$prd$: down researve

$H_{vsg}$: virtual inertia

$D_{vsg}$: virutal damping

### some intern variable

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

Q1: What is $D_t$ in 'TGOV1N' andes?

Q2: Where to find $K_g$?

Q3: How to merge parameters in different andes files to a single table?

Q4: Should add intern variables as indepent variable? 

What if there variable are not included in constriants? i.e. the neural output of MLP, and some binary variables
