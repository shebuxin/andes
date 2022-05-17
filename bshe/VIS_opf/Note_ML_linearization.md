## Machine learning assisted linearization

remember to normalize and denormalize

### frequency response

1. Input and input normalization
   input variable: $M_{sys}$,  $D_{sys}$, $R_{sys}$, $F_{sys}$
   Output: nadir
   input normalization
2. fedfordward and linearlize ReLU

   for each unit of multilayer perseption

   $$
   z' = w*x + b \\
   z  = max \{ z', 0 \}
   $$

   introduce binary variable $a \in \{0, 1\} $, then **both $a$  and $z$ are decision variables**, then

   $$
   z' = w*x + b \\
   z \le z' - h_{down}(1-a) \\
   z \ge z' \\
   z \le h_{up}a \\
   z \ge 0
   $$

   Eliminate $z'$ since x contains gurobi var, then

   $$
   z \le w*x + b  - h_{down}(1-a) \\
   z \ge w*x + b  \\
   z \le h_{up}a \\
   z \ge 0
   $$

where $a \in \{0, 1\} $, $a$ and $z$ are decision variable.

3. output

$$
z_{output} = \sum z * w + b \\

f'_{nadir} = z_{output} * f_{std} + f_{mean}
$$

4. denormlization
$$
f'_{nadir} = \Delta P_e f'_{nadir}
$$

### vsg power response

Compared with frequency prediction, power predictin have two more vairable

1. input variable:
   $M_{sys}$,  $D_{sys}$, $R_{sys}$, $F_{sys}$, $H_{vsg}$, $D_{vsg}$

the rest is similar to those above
