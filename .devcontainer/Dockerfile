# Use the official TensorFlow 2.17 image as the base
FROM tensorflow/tensorflow:2.17.0

# Set build arguments for user ID and username
ARG USER_ID=1000
ARG USER_NAME=defaultuser

# Install necessary packages
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    wget \
    gnuplot \
    libboost-all-dev \
    libx11-dev \
    sudo \
    && rm -rf /var/lib/apt/lists/*

# Add a user with the specified UID and username
RUN useradd --no-log-init -m -u ${USER_ID} -s /bin/bash ${USER_NAME} && \
    echo "${USER_NAME} ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers

# Switch to the new user
USER ${USER_NAME}

# Set the working directory
WORKDIR /home/${USER_NAME}

# Clone the gnuplot-iostream repository
RUN git clone https://github.com/dstahlke/gnuplot-iostream.git

# Build gnuplot-iostream examples (optional)
WORKDIR /home/${USER_NAME}/gnuplot-iostream
#RUN make

# Copy your C++ source code into the container
#COPY . .

# Compile your C++ example
#RUN g++ -std=c++17 -o my_example main.cpp -lboost_iostreams

# Set the entry point to a bash prompt
ENTRYPOINT ["/bin/bash"]
