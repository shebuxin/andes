## Q&A

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

q2: the meaning of "setattr(self, mdl+'dict', mdl_df.T.to_dict())" `<br>`

To add self variable to the class

q3: super().__init__(name) `<br>`

get all the self parameters from father class

q1: why three _build_ function return mdl since they just update the self.mdl? `<br>`

To better shows the process of data passing

__Phase two Q&A__

Q1: What does hasattr( ) mean? has attribute ''cost''? why include ''self and cost'' as input?

```python
        if not hasattr(self, 'cost'):
            self._default_cost()
```

Q2: Why model dict can be updated when there is no model as input? model is a global var?

Q3: function 'build_group_table'

Q4: seems there is no column 'ctrl' at the begining？can add a new column directly?

Q5: actural value used in Andes line? And why there are rate_a, rate_b, and rate_c?

Q6: gsf matrix → two dimential matrix

row:

column:

Q7: meaning of sup function      uncontrolled gen? why formulate sup like this?

Q8: why use p0 as the limit of uncontrollable gen?

Q9: GEN (gendic.key) is var name in grubipy? What about name='pq'? what does 'obj=0' mean?

Q10: how to check mdl is solved?

Q11: Why doesnot output the sedond logger.info  res_cost?



| col1 | col2 | col3 |
| ---- | ---- | ---- |
|      |      |      |
|      |      |      |


想通一个问题：

M D 的都可以写成等式约束，相当于引入中间变量，然后再列写其他等式、不等式约束
