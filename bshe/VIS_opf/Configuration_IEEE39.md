## Case overview

**Standard case**: 39 Bus, 10 generator (connect 9 PV bus and 1 slack bus)

**Modified vis case**: replace 4 SG by turning off existing GENROU (u=0) and adding REGCV1

Existing generator: Bus 30-39

replace SG connected to ***bus 30, 35, 37, and 38***

## Typical parameters

Tips:

1. increase cost of slack bus, since it will compensate for the power mismatch in dynamic simulation

dpe = 0.006 p0

andes base value: $S_{base}$ = 100 MVA

ss.config.mva = 1000

$$
M_{sys} \in [6, 12] \\

D_{sys} \in [0.5, 3] \\

R_{sys} \in [15, 25] \\

F_{sys} \in [0, 20]
$$

TODO:

1. 确认andes的base 然后框定39case各个参数的优化空间

   andes base: 100

   ieee 39 base: 1000
2. 在优化空间内，重新生成训练数据，并训练网络。把神经网络的训练预测结果和matlab仿真结果对比验证
3. 设计ieee 39 的cost，然后验证下dcopf的结果
4. 顺一遍 visopf fnadir的code，包括formulation + code
5. debug 优化结果

   和飞哥确认，gorubypy 能不能看哪些约束bounde，

## Case results

### Base case: DC opf

cost dic: the cost of slack bus gen is higher than the other gen because it need more room for TDS simulation

{'PV_1': {'c2': 0.0, 'c1': 1.0, 'c0': 0.0, 'cr': 0.0, 'cru': 0.0, 'crd': 0.0},
 'PV_2': {'c2': 0.0, 'c1': 1.0, 'c0': 0.0, 'cr': 0.0, 'cru': 0.0, 'crd': 0.0},
 'PV_3': {'c2': 0.0, 'c1': 1.0, 'c0': 0.0, 'cr': 0.0, 'cru': 0.0, 'crd': 0.0},
 'PV_4': {'c2': 0.0, 'c1': 1.0, 'c0': 0.0, 'cr': 0.0, 'cru': 0.0, 'crd': 0.0},
 'PV_5': {'c2': 0.0, 'c1': 1.0, 'c0': 0.0, 'cr': 0.0, 'cru': 0.0, 'crd': 0.0},
 'PV_6': {'c2': 0.0, 'c1': 1.0, 'c0': 0.0, 'cr': 0.0, 'cru': 0.0, 'crd': 0.0},
 'PV_7': {'c2': 0.0, 'c1': 1.0, 'c0': 0.0, 'cr': 0.0, 'cru': 0.0, 'crd': 0.0},
 'PV_8': {'c2': 0.0, 'c1': 1.0, 'c0': 0.0, 'cr': 0.0, 'cru': 0.0, 'crd': 0.0},
 'PV_9': {'c2': 0.0, 'c1': 1.0, 'c0': 0.0, 'cr': 0.0, 'cru': 0.0, 'crd': 0.0},
 'Slack_10': {'c2': 0.0,
  'c1': **1.5**,
  'c0': 0.0,
  'cr': 0.0,
  'cru': 0.0,
  'crd': 0.0}}

Optimal objective  60.064

(        gen         pg  pru  prd
 0      PV_1  12.515789  0.0  0.0
 1      PV_2   5.895637  0.0  0.0
 2      PV_3   8.000000  0.0  0.0
 3      PV_4   5.851096  0.0  0.0
 4      PV_5   7.000000  0.0  0.0
 5      PV_6   4.800000  0.0  0.0
 6      PV_7   3.501478  0.0  0.0
 7      PV_8   7.000000  0.0  0.0
 8      PV_9   1.000000  0.0  0.0
 9  Slack_10   3.000000  0.0  0.0,
     gen  Mvsg  Dvsg
 0  PV_1   0.0   0.0
 1  PV_6   0.0   0.0
 2  PV_8   0.0   0.0
 3  PV_9   0.0   0.0)

TODO: Dynamic index under disturbances

### Base case: DC opf + RoCof

Dpe max = 




# Q&A

1. EV2: M of GENROW_is set as 100

   change to 10
2. the way to replace GENROW with REGCV1

   除了这个删除 GENROU 还有啥要删除的，excitor/IEEEST，IEEEX1似乎和GENROU是match的

   采用u控制离线 在线

   var up/down limit=0

   obj里面乘上u

   计算合成M,D之前 把MD乘上u
3. 负荷的合理性（负载率选多少？）    怎么同时使用在14和39系统里边的

   IEEE39: 1~1.4*P0_base   可行负荷范围和static gen pmax有关

   IEEE14: 0.8-1.6*P0_base
4. set TGGOV1 Tn as generator capacity, otherwise, system base is used
5. cost：gencost分别是1和1.5，reserve cost是0

   check matpower data by searching **'github ieee39 matpower'**
