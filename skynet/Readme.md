# Newer cmake is required

     you can get it here
    wget https://github.com/Kitware/CMake/releases/download/v3.31.5/cmake-3.31.5-linux-x86_64.tar.gz


Then buiild with 
     ../cmake-3.31.5-linux-x86_64/bin/cmake ..


The nbody porgram opens a websocket on localhost:1234 that the webpage html/index.html
Will connect to. It calculates this,

## N-Body simulation

The N-body gravitational problem begins with a set of
N particles with initial positions P<sub>i</sub> = (x<sub>i</sub> , y<sub>i</sub> , z<sub>i</sub>) in three dimensional space, where i ∈ {1, 2, 3, ..., N }, moving with initial velocities V<sub>i</sub> = (v<sub>i</sub><sup>x</sup>, v<sub>i</sub><sup>y</sup>, v<sub>i</sub><sup>z</sup>). Each of these particles are under the influence of a gravitational force F<sub>i</sub> = (f<sub>i</sub><sup>x</sup>, f<sub>i</sub><sup>y</sup>, f<sub>i</sub><sup>z</sup>) due to their masses m i , according to Newton’s inverse-square law of gravitation. The component of the force in the x direction is given by:

<p style="text-align: center;">f<sub>i</sub><sup>x</sup> = GΣ<sub>j≠i</sub> ((m<sub>i</sub>m<sub>j</sub>) / (d<sup>2</sup>(i,j))) ((x<sub>i</sub>x<sub>j</sub>) / (d(i,j)))</p>

where G is the gravitational constant and d is the distance between particles i and j. Similar components of force apply in the y and z directions. The acceleration on a particle, A<sub>i</sub> = (a<sub>i</sub><sup>x</sup>, a<sub>i</sub><sup>y</sup>, a<sub>i</sub><sup>z</sup>), as a result of the force experienced in the x direction is given by:

<p style="text-align: center;">a<sub>i</sub><sup>x</sup> = f<sub>i</sub><sup>x</sup> / m<sub>i</sub></p>

This results in a new velocity v<sub>i</sub><sup>x'</sup> after time δt:

<p style="text-align: center;">v<sub>i</sub><sup>x'</sup> = v<sub>i</sub><sup>x</sup> + a<sub>i</sub><sup>x</sup> δt</p>

The N-body simulation then continues for a predetermined time period T , integrating with δt discrete time steps. In order to achieve accurate results, the time intervals between each integration step must be short.


https://github.com/dileban/nbody-simulation

![alt text](image.png)

The Idea was to use P2300 https://github.com/NVIDIA/stdexec to perform the calculation,
but this is not done yet.