# ams_for_vis

# dynamic_for_vis

### Batch_data_gen_NN_train_v2:

input:   [M, D]

output:  [rocof_max, fnadir, dtheta_max, eig_max]

nn_diam: [16, 128, 4]

### Batch_data_gen_NN_train_v2b:

input:   [M, D]

output:  [rocof_max, fnadir, dtheta_max, eig_max]

nn_diam: [16, 64, 4]

### Batch_data_gen_NN_train_v2c:

divide M and D by 10 based on v2b

### Batch_data_gen_NN_train_v3:

input:   [M, D]

output:  [rocof_max]

nn_diam: [16, 128, 1]
