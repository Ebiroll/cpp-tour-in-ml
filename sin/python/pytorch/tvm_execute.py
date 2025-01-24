import tvm
from tvm import relay
from tvm.contrib import graph_executor
import onnx
import matplotlib.pyplot as plt
import numpy as np

# Load the ONNX model
onnx_model_path = "sin_model.onnx"  # Replace with your ONNX model file
onnx_model = onnx.load(onnx_model_path)

# Convert ONNX model to Relay IR
shape_dict = {"input": (1, 1)}  # Replace with the input shape of your model
mod, params = relay.frontend.from_onnx(onnx_model, shape_dict)

# Print the TIR (Tensor Intermediate Representation)
print("Relay IR:")
print(mod)

# Build the model with LLVM target
target = "llvm"
with tvm.transform.PassContext(opt_level=3):
    lib = relay.build(mod, target=target, params=params)


# Export the compiled model
compiled_model_path = "compiled_model.tar"
lib.export_library(compiled_model_path)
print(f"Compiled model saved to {compiled_model_path}")




# Create a runtime executor
dev = tvm.device("llvm", 0)
module = graph_executor.GraphModule(lib["default"](dev))

# Prepare data for plotting
x_values = np.linspace(0, 2 * np.pi, 100, dtype=np.float32)  # 100 points from 0 to 2π
y_values = []

# Loop over x_values, run inference for each, collect output
for x in x_values:
    # Convert scalar to a batch of shape (1,1)
    input_data = np.array([[x]], dtype="float32")
    module.set_input("input", input_data)
    module.run()
    # Retrieve output from the model
    output = module.get_output(0).asnumpy()
    y_values.append(output.item())  # .item() to get a Python float

# Plot results
plt.figure(figsize=(8, 5))
plt.plot(x_values, y_values, label="Model Output")
plt.title("TVM Model Output from x=0 to x=2π")
plt.xlabel("x")
plt.ylabel("model output")
plt.legend()
plt.grid(True)
plt.show()