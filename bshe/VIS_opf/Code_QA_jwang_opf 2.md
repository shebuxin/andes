## Q&A dcopf

__Phase one Q&A__

Q0: How to check the self variable? `<br>`

class_obj_name.var_name

Q1: How to get the transformed data from andes?

dic data is used for optimization formulation

pandas data is transformed from andes

update dic using function __update_dic__ after modifying any dataframe data

q1: the meaning of "_default_cost" `<br>`

$c_2, c_1, c_0$ is for type I normal generators; `<br>`

$c_r$ is for secondary reserve generators;

q2: the meaning of "setattr(self, mdl+'dict', mdl_df.T.to_dict( ) )" `<br>`

To add self variable to the class

q3: super( ).__init__(name) `<br>`

get all the self parameters from father class

q1: why three _build_ function return mdl since they just update the self.mdl? `<br>`

To better shows the process of data passing

__Phase two Q&A__

Q0: How to set the defalt login in path in server?

google ssh defalt path

Q1: What does hasattr( ) mean? has attribute ''cost''? why include ''self and cost'' as input?

```python
        if not hasattr(self, 'cost'):
            self._default_cost()
```

Q2: Why model dict can be updated when there is no model as input? model is a global var?

default model is the complete gurobypy model → update dic with pandasd dataframe

model can be a __list__ to represent a sepcific parameter

Q3: function 'build_group_table'

a function in andes: 把andes一个group下的器件整合起来

andes group: 一个类型的器件，比如 Static Gen group 包含 PV generator, slack generator

reset_index is a pandas funciton to re-sort the data index

Q4: seems there is no column 'ctrl' at the begining？can add a new column directly?

Add a column to panda dataframe

Q5: actural value used in Andes line? And why there are rate_a, rate_b, and rate_c?

denefition from matpower

rate_a: thermal limit, used in dcopf

rate_b: 。。。rate_c。。。  有个是故障载流，比如短时过载，再查下

Q6: gsf matrix → two dimential matrix

row: row means line, 各个bus在这个line上的转移系数

column: a bus的功率按比例分配到各个line

Q7: meaning of sup function uncontrolled gen? why formulate sup like this?

sup is a kind of load aggration, aggrate load (__negative number__) to each line

set it as a line attribute

Q8: why use p0 as the limit of uncontrollable gen?

set both the min and max as p0, then it is uncontrollable

Q9: GEN (gendic.key) is var name in grubipy? What about name='pq'? what does 'obj=0' mean?

gorubypy.addVars: pg is dic key of the vars, i.e. pg['PV_3']

gorubypy.addVar: set 'name' attribute as the name, the dic key 'pg' is not required when calling the var, i.e. PV_3

obj=0, coefficients in obj, can set it seperately

Q10: how to check mdl is solved?

'X' is an attribute fed by solved gorubypi model

'X' value is the optimizaitn solution

Q11: Why doesnot output the sedond logger.info  res_cost?

info 显示在网页端的白色输出

warning 是红色输出

info 在vscode jupyter 默认不显示

想通一个问题：

M D 的都可以写成等式约束，相当于引入中间变量，然后再列写其他等式、不等式约束

NPCC: robust to Gen trip

VEC: two weak system connection, easy osscialliation

## Q&A rted

**link:** *https://github.com/jinningwang/andes/blob/sfr/jwang/notes2/jams.py#L334*

Q1: **def**to_dcopf(**self**):   calculate previous setting point

Why this function can initiallize a system? (No input system)

Run 'from_andes' before run solve optimization


Q2: input system of rted? No ''from_andes''

data comes from system data 'from_andes'

Manually define **self**.**gen**[**'p_pre'**]


Q3: opf results pru, prd, bu?

pru: up power reserve (dispatch results)

bu: up regulation paticipant factor up reserve

bu = pru / pru.sum()

| gen | pg   | pru   | prd | bu  | bd  |
| --- | ---- | ----- | --- | --- | --- |
| 0   | PV_2 | 0.100 | 0.0 | 0.0 | 0.0 |

line 449, 451, why equality constraints?

line 619, 621

1) step one: run dispatch

2. step two: run TDS, keep generating AGC

3. step three: AGC is dispatched to gen accordin go bu and bd

4. step four: run dispatch (every 5 minutes)


Q4: what is sfr du and dd? 

sfr requirements for up and down regulation, not necessarily met

estimated by operator, i.e. ISO


How to formulate the pru and prd constraints? 

rted: hard contraints for sfr regulation


Q5(rted2): gen_idx for type II is manually defined? 

mannually defined


Q6: **for** (**new_key**, **new_value**) **in**gendict.**items**()

dict is stored in a dict


Q7: soft sfr?

remove sfr du and dd constraints

add violation cost to obj


## Q&A  Co-sim ieee14

**Link:** *https://github.com/jinningwang/andes/blob/sfr/jwang/notes2/co-sim_ieee14.ipynb*

Q1: Where to find synthetic load?

From PJM: *https://dataminer2.pjm.com/list*

---

用同样的load profile数据

比较economic dispatch 加入不同的constraints 会带来怎样不同的优化结果

敏感度分析：1）cost 对 load 的灵敏度 → LMP : locational marginal price

2) cost 对 其他参数的灵敏度
   比如：cost 对 (sfr, inertia support) reserve 的大小的灵敏度，看看有没有critical change，类似于CLL ciritical load level
3) 而reserve的大小又是由 inertia 的参数或者sfr的参数决定的，分析对这些参数的灵敏度

***关键点：ED加进新的cosntraints，就会产生critical state，可以进行灵敏度分析***

develop analytical 的框架?


ED里面0，1变量不可导：0、1变量找不到对偶，不可导

---


Q2: gen_cost matirx (in [7]) seems different from that in dcopf

gen_cost is matpower format，use the last three column


目前 二次 一次都有的结果需要validate

只用二次或者只用一次结果正确


Q3: dataframe? set uncontrollable gen before initialize opf class?

```
# set EV generator as uncontrollable
ssp.gen.controllable.iloc[4] = False
```

iloc: index location, 就是row


Q4: Some functions imported from andes

```
add_gencost, build_group_table, make_GSF
to_pandapower, make_link_table, runopp_map
```

make line table

runopp_map → run with acopf with pandapower and map to andes


Q5: How to calculate AGC signal? Use PI?

PI regu in TDS

accumulated AGC, dispatch AGC results every every 4s


Q6: Tool to run acopf? data source?

pandapower, update load file before running opf


Q7: Smooth setpoint?

smoote setpoint (calculated by economic dispatch) for TDS


Q8: use andes build in plot or matplot?  how to set different colors when break through the limit ?
