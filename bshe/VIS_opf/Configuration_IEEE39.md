## Case overview

**Standard case**: 39 Bus, 10 generator (connect 9 PV bus and 1 slack bus)

**Modified vis case**: replace 4 SG by turning off existing GENROU (u=0) and adding REGCV1

Existing generator: Bus 30-39

replace SG connected to ***bus 30, 35, 37, and 38***

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
