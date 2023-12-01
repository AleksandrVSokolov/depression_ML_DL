import tensorflow as tf


#################################################################################
##########################    Section    ########################################
#################################################################################
# Small classifiers with different parameters

def DNNclassfier(input_data, scaler=1, inner_activations="relu", outer_activation="sigmoid"):

    hidden_layer1 = tf.keras.layers.Dense(units=100*scaler, activation=inner_activations, 
                                          kernel_regularizer="l2", name="hidden_layer1")(input_data)
    hidden_layer2 = tf.keras.layers.Dense(units=64*scaler, activation=inner_activations, 
                                          kernel_regularizer="l2", name="hidden_layer2")(hidden_layer1)
    hidden_layer3 = tf.keras.layers.Dense(units=32*scaler, activation=inner_activations, 
                                          kernel_regularizer="l2", name="hidden_layer3")(hidden_layer2)
    hidden_layer4 = tf.keras.layers.Dense(units=16*scaler, activation=inner_activations, 
                                          kernel_regularizer="l2", name="hidden_layer4")(hidden_layer3)
    output_layer = tf.keras.layers.Dense(units=1, activation=outer_activation, 
                                         kernel_regularizer="l2", name="output")(hidden_layer4)

    return output_layer


def DNNclassfiernoreg(input_data, scaler=1, inner_activations="relu", outer_activation="sigmoid"):

    hidden_layer1 = tf.keras.layers.Dense(units=100*scaler, activation=inner_activations, name="dense_layer1")(input_data)
    hidden_layer2 = tf.keras.layers.Dense(units=64*scaler, activation=inner_activations, name="dense_layer2")(hidden_layer1)
    hidden_layer3 = tf.keras.layers.Dense(units=32*scaler, activation=inner_activations, name="dense_layer3")(hidden_layer2)
    hidden_layer4 = tf.keras.layers.Dense(units=16*scaler, activation=inner_activations, name="dense_layer4")(hidden_layer3)
    output_layer = tf.keras.layers.Dense(units=1, activation=outer_activation, name="output")(hidden_layer4)

    return output_layer


def DNNclassfierL1L2(input_data, scaler=1, inner_activations="relu", outer_activation="sigmoid"):

    hidden_layer1 = tf.keras.layers.Dense(units=100*scaler, activation=inner_activations, 
                                        kernel_regularizer="l1_l2", name="dense_layer1")(input_data)
    hidden_layer2 = tf.keras.layers.Dense(units=64*scaler, activation=inner_activations, 
                                        kernel_regularizer="l1_l2", name="dense_layer2")(hidden_layer1)
    hidden_layer3 = tf.keras.layers.Dense(units=32*scaler, activation=inner_activations, 
                                        kernel_regularizer="l1_l2", name="dense_layer3")(hidden_layer2)
    hidden_layer4 = tf.keras.layers.Dense(units=16*scaler, activation=inner_activations, 
                                        kernel_regularizer="l1_l2", name="dense_layer4")(hidden_layer3)
    output_layer = tf.keras.layers.Dense(units=1, activation=outer_activation, 
                                         kernel_regularizer="l1_l2", name="output")(hidden_layer4)

    return output_layer


def DNNclassfierL1L2drop(input_data, scaler=1, inner_activations="relu", outer_activation="sigmoid"):

    hidden_layer1 = tf.keras.layers.Dense(units=100*scaler, activation=inner_activations, 
                                          kernel_regularizer="l1_l2", name="dense_layer1")(input_data)
    hidden_layer2 = tf.keras.layers.Dense(units=64*scaler, activation=inner_activations, 
                                          kernel_regularizer="l1_l2", name="dense_layer2")(hidden_layer1)

    dropout_1 = tf.keras.layers.Dropout(0.15, name="drop_layer1")(hidden_layer2)

    hidden_layer3 = tf.keras.layers.Dense(units=32*scaler, activation=inner_activations, 
                                        kernel_regularizer="l1_l2", name="dense_layer3")(dropout_1)
    hidden_layer4 = tf.keras.layers.Dense(units=16*scaler, activation=inner_activations, 
                                        kernel_regularizer="l1_l2", name="dense_layer4")(hidden_layer3)
    output_layer = tf.keras.layers.Dense(units=1, activation=outer_activation, 
                                        kernel_regularizer="l1_l2", name="output")(hidden_layer4)

    return output_layer


def DNNclassfierL1L2manydrops(input_data, scaler=1, inner_activations="relu", outer_activation="sigmoid"):

    x = tf.keras.layers.Dense(units=100*scaler, activation=inner_activations,
                              kernel_regularizer="l1_l2", name="dense_layer1")(input_data)
    x = tf.keras.layers.Dropout(0.1, name="drop_layer1")(x)

    x = tf.keras.layers.Dense(units=64*scaler, activation=inner_activations,
                              kernel_regularizer="l1_l2", name="dense_layer2")(x)
    x = tf.keras.layers.Dropout(0.1, name="drop_layer2")(x)

    x = tf.keras.layers.Dense(units=32*scaler, activation=inner_activations,
                              kernel_regularizer="l1_l2", name="dense_layer3")(x)
    x = tf.keras.layers.Dropout(0.1, name="drop_layer3")(x)

    x = tf.keras.layers.Dense(units=16*scaler, activation=inner_activations,
                              kernel_regularizer="l1_l2", name="dense_layer4")(x)
    x = tf.keras.layers.Dropout(0.1, name="drop_layer4")(x)

    output_layer = tf.keras.layers.Dense(units=1, activation=outer_activation,
                                         kernel_regularizer="l1_l2", name="output")(x)

    return output_layer


def DNNclassfierL1L2manyBigdrops(input_data, scaler=1, inner_activations="relu", outer_activation="sigmoid"):

    x = tf.keras.layers.Dense(units=100*scaler, activation=inner_activations,
                              kernel_regularizer="l1_l2", name="dense_layer1")(input_data)
    x = tf.keras.layers.Dropout(0.3, name="drop_layer1")(x)

    x = tf.keras.layers.Dense(units=64*scaler, activation=inner_activations,
                              kernel_regularizer="l1_l2", name="dense_layer2")(x)
    x = tf.keras.layers.Dropout(0.3, name="drop_layer2")(x)

    x = tf.keras.layers.Dense(units=32*scaler, activation=inner_activations,
                              kernel_regularizer="l1_l2", name="dense_layer3")(x)
    x = tf.keras.layers.Dropout(0.3, name="drop_layer3")(x)

    x = tf.keras.layers.Dense(units=16*scaler, activation=inner_activations,
                              kernel_regularizer="l1_l2", name="dense_layer4")(x)
    x = tf.keras.layers.Dropout(0.3, name="drop_layer4")(x)

    output_layer = tf.keras.layers.Dense(units=1, activation=outer_activation,
                                         kernel_regularizer="l1_l2", name="output")(x)

    return output_layer


def DNNclassfierOnlyDrops(input_data, scaler=1, inner_activations="relu", outer_activation="sigmoid"):

    x = tf.keras.layers.Dense(units=100*scaler, activation=inner_activations, name="dense_layer1")(input_data)
    x = tf.keras.layers.Dropout(0.1, name="drop_layer1")(x)

    x = tf.keras.layers.Dense(units=64*scaler, activation=inner_activations, name="dense_layer2")(x)
    x = tf.keras.layers.Dropout(0.1, name="drop_layer2")(x)

    x = tf.keras.layers.Dense(units=32*scaler, activation=inner_activations, name="dense_layer3")(x)
    x = tf.keras.layers.Dropout(0.1, name="drop_layer3")(x)

    x = tf.keras.layers.Dense(units=16*scaler, activation=inner_activations, name="dense_layer4")(x)
    x = tf.keras.layers.Dropout(0.1, name="drop_layer4")(x)

    output_layer = tf.keras.layers.Dense(units=1, activation=outer_activation, name="output")(x)

    return output_layer


def DNNclassfierSimple(input_data, scaler=1, inner_activations="relu", outer_activation="sigmoid"):

    x = tf.keras.layers.Dense(units=100*scaler, activation=inner_activations, name="dense_layer1")(input_data)
    x = tf.keras.layers.Dropout(0.15, name="drop_layer1")(x)

    x = tf.keras.layers.Dense(units=32*scaler, activation=inner_activations, name="dense_layer2")(x)

    output_layer = tf.keras.layers.Dense(units=1, activation=outer_activation, name="output")(x)

    return output_layer


def DNNclassfierSimpleBatch(input_data, scaler=1, inner_activations="relu", outer_activation="sigmoid"):

    x = tf.keras.layers.Dense(units=100*scaler, activation=inner_activations, name="dense_layer1")(input_data)
    x = tf.keras.layers.Dropout(0.1, name="drop_layer1")(x)

    x = tf.keras.layers.Dense(units=64*scaler, name="dense_layer2")(x)
    x = tf.keras.layers.BatchNormalization(name="BN1")(x)
    x = tf.keras.layers.Activation(inner_activations)(x)
    x = tf.keras.layers.Dropout(0.1, name="drop_layer2")(x)

    x = tf.keras.layers.Dense(units=16*scaler, name="dense_layer3")(x)
    x = tf.keras.layers.BatchNormalization(name="BN2")(x)
    x = tf.keras.layers.Activation(inner_activations)(x)

    output_layer = tf.keras.layers.Dense(units=1, activation=outer_activation, name="output")(x)

    return output_layer


#################################################################################
##########################    Section    ########################################
#################################################################################
# Inner classifier for joint models
def innerDNNClassifier(input_dim1D,
                       BatchNorm,
                       inner_activations="relu",
                       outer_activation="sigmoid",
                       DropRate=0,
                       **kwargs):
    
    if BatchNorm:

        def dense_block(x,
                        nodes,
                        inner_activations=inner_activations,
                        DropRate=DropRate,
                        **kwargs):
            x = tf.keras.layers.Dense(units=nodes, **kwargs)(x)
            x = tf.keras.layers.BatchNormalization()(x)
            x = tf.keras.layers.Activation(inner_activations)(x)
            x = tf.keras.layers.Dropout(DropRate)(x)
            return x
    else: 

        def dense_block(x,
                        nodes,
                        inner_activations=inner_activations,
                        DropRate=DropRate,
                        **kwargs):
            x = tf.keras.layers.Dense(units=nodes, activation=inner_activations, **kwargs)(x)
            x = tf.keras.layers.Dropout(DropRate)(x)
            return x
        
    x = tf.keras.layers.Dropout(DropRate)(input_dim1D)
    x = dense_block(x, nodes=32, **kwargs)
    x = dense_block(x, nodes=32, **kwargs)
    x = dense_block(x, nodes=16, **kwargs)
    x = tf.keras.layers.Dense(units=1, activation=outer_activation, name="out", **kwargs)(x)

    return x