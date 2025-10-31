# XY spin glass model on a 3D cubic lattice
In this project, I studied the XY spin glass model on a 3D cubic lattice using advanced Monte Carlo simulation techniques. The system consists of planar spins with continuous angular degrees of freedom, interacting via random Gaussian couplings, leading to frustration and disorder typical of spin glasses. The numerical implementation combined Metropolis updates, over-relaxation moves, and parallel tempering (replica exchange) to efficiently sample the rugged energy landscape and equilibrate at low temperatures. The main focus was on investigating the thermodynamic behavior and searching for signs of a spin glass transition in three dimensions, such as the behavior of the spin overlap distribution. The simulations were implemented from scratch in C++ and analyzed using Python, with attention to finite-size scaling and thermalization diagnostics.

# How the code works:
The file montecarlo.hpp contains the core of the simulation, which relies on the "randnumgen.hpp" library to generate pseudo-random numbers. The "main.cpp" file calls these components to execute the program. Finally, an example configuration file named "config.txt" is included, from which the program reads the system control parameters.

# General structure of the simulation
| **Class**       | **Description** |
|------------------|-----------------|
| `sys` | Represents a single physical system at a fixed temperature. It stores the spin configuration, the coupling matrix, and provides the main update routines (Metropolis and over-relaxation). Each `sys` object can evolve independently and output its energy and configuration. |
| `replica` | Manages a set of systems at different temperatures and performs **Parallel Tempering** exchanges between them. This allows the system to overcome energy barriers and reach equilibrium faster. Each replica evolves all its temperature points in parallel using multithreading. |
| `replicas2` | Contains two independent replicas that share the same disorder realization (same $J_{ij}$). It is mainly used to compute overlap observables and correlation lengths related to spin-glass ordering. |
| `montecarlo` | Handles multiple disorder samples. Each sample corresponds to a different random realization of the couplings, and the class averages physical quantities over all samples to produce quenched averages. |

For further details on the architecture, you can check the "XY spin glass model.pdf", where each method and class is explained in detail.
---
If you have any questions about the implementation, would like further clarifications, or have suggestions for improvement, feel free to reach out — I’ll be happy to discuss and help.

**Mattia Liguoro**    
mattialiguoro17@gmail.com
