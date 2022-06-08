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

---

# Code Q&A 1

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

# Code Q&A 2

Q1: Co-sim 里的Benchmark结果似乎不一样

A1: 多解造成的，cost需要设置的不一样

---

Q2: relationship between gurobi and pandapower

A2: gorubi solves the reserve capacity using dcopf, and change the power limit in pandapower

pandapower AC opf package include:

- dcopf
- acopf → line limit, ramping limit

---

Q3: link table doesn't include vsg gen

Q3: solved.  ref pdpower.py

---

Q4: function *get_pe*? Why not get the pe of GENROW?

get pe or Pref2 in REGCV1?

Q5: get_pre is use to limit ramping limit

ignore it if not consier ramping limit

---

Q5: return of get load function

    - dpd_u: regu up reauirement, for sfr reserve

    - dpd_d: regu down requirement, for sfr reserve

    - load_exp: load forecasted value

vsg reserve can match dpd_u and dpd_d

ref def_sfr in jwang's code

---

Q6: what is step change interval?

intv_step = 100 #  step change; smooth the setpoitns

---

Q7: What is runopp_map ?
A7: pandapwoer uses dataframe index, i.g. 0, 1, 2 ....
andes uses idx to call variable

runopp_map forces pandapower to calll variable using andes idx

---

Q8: what is "bu" and "bd" in rted get_res( ) and used in AGC?

A8: agc parcipation factor power

... set more room for sfr

把一个gen的cost调整最大，用来承接 AGC mismatch

---

Q9: How to set agc signal to vsg gen?

A9:

---

Q10: 如果相邻两次的dispatch结果相差较大，影响TDS吗？

比如：governor的ref改的过大，这在vsg capacity变化时会比较明显

A 10: To test

smooth in 100s

5% total load

PV profile, check fangxing's paper

**set large disturbance around 150s in ED to avoid the smooth progress**

---

Q11: save final TDS results to file?

A11: In progress ...


Q12: difference of pandapower sn_mva and andes sn_mva?



---
