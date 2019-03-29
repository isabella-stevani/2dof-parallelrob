# Robust Control of a 2-DOF Parallel Mechanism

## Overview

This project refers to the master thesis entitled "Robust Control of a 2-DOF Parallel Mechanism Combining Feedback Linearization and H-infinity Design" developed by Isabella Stevani with the supervision of Professor Diego Colón, both from the Laboratory of Automation and Control of the University of São Paulo - LAC/USP. More details about the control and modeling methodologies used in the codes presented here could be found at the master thesis qualification text available soon at Research Gate.

The purpose of this project is to combine feedback linearization and H-infinity design to achieve a robust control design with simple tuning and applicable to a wide range of nonlinear systems, specially robotic mechanisms with parallel architecture. The purposed robust control design must stabilize the system even in the presence of uncertain parameters in the mechanism model.

#### -- Project Status: Active

## Code

The code was developed in MATLAB 2015a, so early versions might be incompatible.

### Instructions

- Clone the repository:
```sh
$ git clone https://github.com/isabella-stevani/2dof-parallelrob.git
```

- In MATLAB, make the cloned repository the current folder:
```sh
>> cd ..\2dof-parallelrob
```

- For continuous time simulation:
```sh
>> run('main_c.m')
```

- For discrete time simulation:
```sh
>> run('main_d.m')
```

## Project Details

### Partners
* André Garnier Coutinho, GMSIE/USP

### Methods Used
* Parallel mechanism modeling method developed by André Garnier Coutinho
* Feedback linearization
* H-infinity design
* 4th order Runge-Kutta

### Contributing Members

**Team Lead: [Isabella Stevani](https://github.com/isabella-stevani)**

## Contact
* Feel free to contact team leads with any questions or if you are interested in contributing!
