# 整体思路

根据下一个interval (5mins, 300s)的 preidcted average load curve 做 economic dispatch, 然后改setting point

在接下来的300s内不断有 load change 发生，其中有个时刻可以设置 5%的大波动，然后圈出来观察 primary frequency response

# My code framework

There are two attrbutes

    - andes attributes: ssa

    - visopf attribute: ssv

ssv is updated based on new ssa, ssv.from_andes(ssa)

'andes.to_list' outputs the latest data after calling ssa.PV.p0 = p0_new

---

for time in range(end_time):

    - load update t_inv = 1s

    modify andes TDS attribute, ssa.PQ

    set a large load deviation in load file around 100s after Dispatch update

    - AGC signal update t_inv = 4s

    modify andes TDS attribute, ssa.TurbineGov.set

    - Dispatch update t_inv = 300s

    get latest state from andes TDS attribute, including ssa.PQ

    then, update visopf attribute using ssv.from_andes(ssa), and new pmax

    new setting point for andes attribute ssa.PV.p0

# Code framework of jwang co-sim39_ev

1. load andes case

   - set output config as manual
   - set constant PQ load
   - turn on 'numba' to acclerate andes
2. load synthetic for 1 hour

   - real time load ratio, interval = 1s
   - synthetic laod raito for dispatch (average load each 5 mins)
3. build instance

   - pandapower instance: converted from andes instance
   - build RTED instance → opf class we wrote

   a key action: set p_pre `<br>`
4. some preparation

   - link talbe: table linked pandapower and opf attribute
   - function get_pre: Get the active power (TurbineGov/DG) after TDS
   - Calc SFR requirements
5. set TDS parameters

   - TDS loop para
   - AGC table
   - ACE var
   - ini load
6. TSD loop

TODO:

add to dc

add M, D

add set_ppre

---

# Code Q&A

Q1: the functionality of to_dcopf( )? returns a gurobi attribute, and can be directly solved?

A1: a gurobi class, can call build and get_res( )

---

Q2: What is make_link_table(andes_attribute)?

A2: Map pandapower to andes

---

Q3: setpoints for DG?

for example: agc for DG

PVD1 (ideal curent source), EV, EV2

---

Q4: In cell[20], Reserve some capacity to avoid TDS crush?
ssp.gen.max_p_mw = ssp.gen.max_p_mw

A4: reserve power for TDS

---

Q5: Does andes continue run or run from the beginning?
ssa.TDS.config.tf = end_time

A5: 往后跑

---

Q6: The functionality of 'smoth setting point'?

A6: limit the AGC setting point single sharp increase

---

Q7: 哪些量是不允许改变的？load, gov, PV,

load scale 和 接入toggler 有区别吗？

A7: REGC.M, D可以

潮流分布不同

single load change 更critical

---

Q8: opf 的结果给 governor or PV 的 p0

PV的p0 是算power flow用的

A8: ofp 给governor

---

Q9: go through AGC code

A9: andes build in AGC block

---

Q10: Plan A is better?
