
Steps to build tensorflow and this example

git submodule update --init --recursive

# Run in codespace
I have set up the dockefile so that it should run in github codespaces

If you are not, you might want to edit the Dockerfile
// Add a user with the specified UID and username
   RUN useradd --no-log-init -m -u ${USER_ID} -s /bin/bash ${USER_NAME} && \
       echo "${USER_NAME} ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers

This is the output you are looking for. 
```
________                               _______________
___  __/__________________________________  ____/__  /________      __
__  /  _  _ \_  __ \_  ___/  __ \_  ___/_  /_   __  /_  __ \_ | /| / /
_  /   /  __/  / / /(__  )/ /_/ /  /   _  __/   _  / / /_/ /_ |/ |/ /
/_/    \___//_/ /_//____/ \____//_/    /_/      /_/  \____/____/|__/


You are running this container as user with ID 1000 and group 1000,
which should map to the ID and group for your user on the Docker host. Great!
```
If you are having problems, try
export USER_ID=$(id -u)


# Another option
I never tried this, but this like another reasonable path forward.

https://github.com/FloopCZ/tensorflow_cc


# Steps to build

This is not necessary if you cloned the submodule
Clone the TensorFlow Repository:

    git clone https://github.com/tensorflow/tensorflow/



Install dependencies:
    sudo apt-get install build-essential cmake git curl zip unzip autoconf automake libtool  gcc g++
    sudo apt install libboost-all-dev

Bazel is used to build tensorflow

    wget https://github.com/bazelbuild/bazelisk/releases/download/v1.25.0/bazelisk-linux-amd64
    mv this to i.e. $HOME/bin/bazel if it is in your path


Configure the Build Environment:

    Run the ./configure script to set up the build environment, specifying appropriate options for your system.


Build library
    Executed the Bazel build command to compile the TensorFlow C++ library:

bazel build --config=opt //tensorflow:libtensorflow_cc.so



mkdir build
cd build
cmake ..
make
./tf_sine_example




Some external libraries when building tensorflow, 

ls bazel-bin/external/
boringssl             com_google_protobuf        curl               fft2d        FXdiv        icu            llvm-project           nsync       png          snappy      upb
com_github_grpc_grpc  com_googlesource_code_re2  double_conversion  flatbuffers  gif          jsoncpp_git    local_config_cuda      onednn      pthreadpool  stablehlo   XNNPACK
com_google_absl       cpuinfo                    farmhash_archive   FP16         highwayhash  libjpeg_turbo  local_config_tensorrt  org_sqlite  ruy          tf_runtime  zlib

Try builing headers,

bazel build //tensorflow:install_headers


# Compiling with conda

As I had conda installed, I tried but never got it to work.


# Never to early to give up,

Never got it to work so I switched to using a Dockerfile

docker build --build-arg USER_ID=$(id -u) --build-arg USER_NAME=$(whoami) -t tensorflow_cpp_env .devcontainer/


'''
xhost +local:docker
docker run -it \
    -e DISPLAY=$DISPLAY \
    -v /tmp/.X11-unix:/tmp/.X11-unix \
    -v $(pwd):/workspace \
    --user $(id -u):$(id -g) \
    tensorflow_cpp_env
'''

# Cleaning up the container
   docker rmi tensorflow_cpp_env

# bare metal

Inspired by this article, 

https://www.codeproject.com/Articles/1237026/Simple-MLP-Backpropagation-Artificial-Neural-Netwo

sin/bare_metal.cpp

It uses simple toy-example gradient descent


```
     	Input Layer      	Hidden Layer (N neurons)        	Output Layer
     	┌─────────┐      	 ┌────────────────────────┐      	┌──────────┐
     	│     	  │      	│                    	  │      	│      	   │
     	│         │  W[0]	│  c[0] + W[0]*x          │  V[0]	│      	   │
     	│     	  ├──────────► (sigmoid activation)   ├──────────►         │
     	│     	  │      	│                    	  │      	│      	   │
     	│     	  │  W[1]	│  c[1] + W[1]*x     	  │  V[1]	│  Σ(V[i]  │
     	│	x	  ├──────────► (sigmoid activation)   ├──────────► *h[i])+ │
     	│     	  │      	│                    	  │         │	b 	   │
     	│     	  │   ...	│      	...       	      │   ...   │      	   │
     	│     	  │      	│                    	  │      	│      	   │
     	│     	  │  W[N-1] │  c[N-1] + W[N-1]*x 	  │  V[N-1] │      	   │
     	│     	  ├──────────► (sigmoid activation)   ├──────────►         │
     	└─────────┘      	└────────────────────────┘      	└──────────┘
```

# label_image

Example from tensorflow codebase, it has its own Readme

# eigen 

eigen/simple_nn.cpp

Simple backpropagation network using Eigen

https://medium.com/@chizy7/implementing-a-simple-neural-network-in-c-with-eigen-f555a664f7b8

Implementing a neural network leraning Xor with  eigen

eigen/main.cpp
eigen/NeuralNetwork.cpp

# skynet
Testing with a n-body simulation

# httpserver
Senders example from net29
https://github.com/bemanproject/net

Senders For Network Operations

# Deep learning examples
Here are many examples of deep learning models
https://github.com/rasbt/deeplearning-models