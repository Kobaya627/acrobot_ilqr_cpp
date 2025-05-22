# acrobot_ilqr_cpp
Acrobot swing up task, using Iterative Linear Quadratic Regulator(ilqr) for offline trajectory optimization and Time-varing Linear Quadratic Regulator(tvlqr) for online trajectory tracking

![swing_up_mujoco](./assets/mujoco_swing_up.gif)

## Download & Build

```
git clone https://github.com/Kobaya627/acrobot_ilqr_cpp.git
cd acrobot_ilqr_cpp
mkdir build && cd build
cmake ..
make
```

2 executable files(**acrobot_ilqr_node** and **acrobot_mujoco_node**) will be generated under current directory

## Getting Started

### Offline Swing Up Trajectory Optimization with ilqr
```
./acrobot_ilqr_node
```
Files containing optimized swing up trajectory will be created in `../data/acrobot_ilqr/`

### Online Trajectory Tracking with tvlqr 

```
./acrobot_mujoco_node
```
