import numpy as np
import matplotlib.pyplot as plt
import tvm
import tensorflow as tf
from tvm import relay
from tvm.contrib import graph_executor

# This import may vary depending on your environment:
# - If using full TensorFlow: `import tensorflow as tf`
# - If using tflite-runtime: `from tflite_runtime import interpreter as tflite_interpreter`
# - If installed via pip install tflite, you can do:
import tflite  # For reading .tflite files

#####################################
# 1. Load the TFLite model from file
#####################################
tflite_model_file = "tf_sine_model.tflite" 
with open(tflite_model_file, "rb") as f:
    tflite_model_buf = f.read()

# Parse the TFLite model with the FlatBuffers-based TFLite module
tflite_model = tflite.Model.GetRootAsModel(tflite_model_buf, 0)



# Load the TFLite interpreter to obtain input/output details
interpreter = tf.lite.Interpreter(model_path=tflite_model_file)
interpreter.allocate_tensors()
input_details = interpreter.get_input_details()
for i in range(len(input_details)):
    print("== Input details for: " + tflite_model_file + " ==")
    input_name=input_details[i]['name']
    print("name:", input_name)
    print("shape:", input_details[i]['shape'])
    print("type:", input_details[i]['dtype'])
    
#output_details = interpreter.get_output_details()
#for i in range(len(output_details)):
#    print("==== output_details details for: " + tflite_model_file + " ====")
#    output_name=output_details[i]['name']
#    print("i,name:", i, output_name)
#    print("shape:", output_details[i]['shape'])
#    print("type:", output_details[i]['dtype'])




#####################################
# 2. Convert TFLite model to Relay IR
#####################################
# Change these based on your TFLite model's expected input name/shape/dtype
input_tensor_name = "serving_default_keras_tensor_38:0"
input_shape = (1, 1)
input_dtype = "float32"

print("Converting TFLite model to Relay IR...")
mod, params = relay.frontend.from_tflite(
    tflite_model,
    shape_dict={input_tensor_name: input_shape},
    dtype_dict={input_tensor_name: input_dtype}
)


print(mod)

##################################
# 3. Compile the Relay IR with TVM
##################################
target = "llvm"  # You can also use "llvm -mcpu=core-avx2" or other targets
with tvm.transform.PassContext(opt_level=3):
    lib = relay.build(mod, target=target, params=params)


###################################################
# 4. Create a runtime executor and run in a loop 0-2π
###################################################
# Create a device (CPU) and a GraphModule
dev = tvm.device(target, 0)
module = graph_executor.GraphModule(lib["default"](dev))

# Generate 100 points between 0 and 2π
x_values = np.linspace(0, 2 * np.pi, 100, dtype=np.float32)
y_values = []

# Inference loop
for x in x_values:
    # Model expects shape [1,1], so wrap x in that shape
    input_data = np.array([[x]], dtype=np.float32)
    module.set_input(input_tensor_name, input_data)
    module.run()

    # Get output data, assume it's shape [1,1]
    output_data = module.get_output(0).asnumpy()
    y_values.append(output_data.squeeze())  # Convert array to scalar

# Compute the actual sine for comparison
sin_values = np.sin(x_values)

    
###############################
# 5. Plot the model predictions
###############################
plt.figure(figsize=(8, 5))
plt.plot(x_values, y_values, label="TFLite Model Output (compiled by TVM)")
plt.plot(x_values, sin_values, label="sin(x)", color="red")
plt.title("TVM-Compiled TFLite Model Inference")
plt.xlabel("x")
plt.ylabel("model output")
plt.legend()
plt.grid(True)
plt.show()

