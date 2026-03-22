## Overview
The project models the Universal Robots UR5, a 6-degree-of-freedom (6-DOF) industrial manipulator, using both the traditional Denavit-Hartenberg (DH) convention and the modern Screw Theory (Product of Exponentials) framework.

### Key Comparative Modules
Forward Kinematics: Validates end-effector position and orientation across cyclic trajectories.

Inverse Kinematics: Employs numerical optimization (fmincon) to solve for joint angles using both models.

Jacobian Analysis: Compares the Geometric Jacobian (DH) and Body Jacobian (Screw Theory) through singular value decomposition.

Static Dynamics: Calculates gravity-induced torques at each joint based on simplified mass properties.

Robot Visualization: Includes a real-time animation of the UR5 performing 3D motion.

## Requirements
MATLAB (Tested on R2021b or later).

Optimization Toolbox: Required for the Inverse Kinematics solver (fmincon).

## How to Use
Download the .m script

Run the script in the MATLAB environment.

Outputs:

Figure 1: Displays the 3D trajectory, error analysis, IK comparisons, and gravity torques.

Figure 2: An interactive animation of the UR5 robot model in motion.

Command Window: Prints a results summary including position/orientation errors and Jacobian condition numbers.

## Core Conclusions

The simulation demonstrates that both methods yield identical kinematic results within numerical precision. However, the script highlights that Screw Theory offers several advantages:

Geometric Intuition: Easier to visualize axes in the space frame.

Mathematical Elegance: Simpler Jacobian derivation and a natural extension to robot dynamics.

Efficiency: Global representation that avoids the tedious link-by-link frame assignments required by DH.

## Author
Sarthak Bikram Panta
