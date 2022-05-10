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

link: https://github.com/jinningwang/andes/blob/sfr/jwang/notes2/jams.py#L334

Q1: **def**to_dcopf(**self**):   calculate previous setting point

Why this function can initiallize a system? (No input system)

Q2: input system of rted? No ''from_andes''

Q3: what is sfr du and dd?

line 449, 451, why equality constraints?

line 619, 621

Q4: How to formulate the pru and prd constraints? Line 454 should be '+'?

Q5(rted2): gen_idx for type II is manually defined? where 

Q6: **for** (**new_key**, **new_value**) **in**gendict.**items**()

what is dic.item?

Q7: soft sfr?
