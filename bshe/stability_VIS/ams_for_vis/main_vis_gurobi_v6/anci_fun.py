import torch
import numpy as np

def extract_nn_para(nn_name):
    """
    Extract the parameters of a fully connected neural network
    using torch
    
    input
    ---------
        nn_name: string
            the name of the neural network
                  name

    return
    ---------
        nn_weight: list
            the weight of the neural network
        nn_bias: list
            the bias of the neural network
        nn_dim: list
            the dimension of the neural network
    """

    # load the neural network
    nn_para = torch.load(nn_name)

    # extract the weight and bias
    nn_weights = []
    nn_bias = []
    nn_dim = []
    items_list = list(nn_para.items())
    
    for index, (key, value) in enumerate(items_list):
        if index % 2 == 0:
            nn_weights.append(value.numpy().astype(np.float64))
        else:
            nn_bias.append(value.numpy().astype(np.float64))
            # get the dimension of the layer
            dim = nn_bias[-1].shape[0]
            nn_dim.append(dim)
    
    # get the input dimension
    input_dim = nn_weights[0].shape[1]
    # add the input dimension to the nn_dim
    nn_dim.insert(0, input_dim)

    return nn_weights, nn_bias, nn_dim