🌀 Computational-Mechanics-Simulation
Advanced Process Engineering Portfolio | Summer Semester 2025 Developed at the Chair of Mechanical Process Engineering, Otto-von-Guericke University Magdeburg.

📋 Project Summary
This repository houses a suite of numerical and stochastic models designed to analyze complex mechanical operations. The focus is on predicting particulate behavior and fluid-solid interactions through high-fidelity simulations. These tools were developed to bridge the gap between theoretical process physics and practical computational implementation.

🛠 Simulation Modules
1. Deterministic Trajectory Analysis (Winnowing)
Methods: Euler & Runge–Kutta Integration This module maps the precise flight paths of particles within aerodynamic separators. By comparing first-order and higher-order integration schemes, the simulation evaluates numerical stability and the accuracy of drag-force calculations in a winnowing environment.

2. Stochastic Motion Modeling (Monte Carlo)
Method: Probabilistic Sampling Recognizing that real-world particles are never identical, this implementation uses Monte Carlo techniques to inject randomness into particle properties (mass, size, and initial velocity). The result is a statistical distribution of trajectories, providing a more realistic view of process efficiency.

3. Particulate Interaction Dynamics (DEM)
Framework: Discrete Element Method (DEM) A comparative study of collision physics through two distinct lenses:

Hard-Sphere Model: Analyzing instantaneous momentum exchange during ideal elastic impacts.

Soft-Sphere Model: Implementing realistic contact mechanics to simulate material deformation and energy dissipation during particle-to-particle collisions.

🚀 Technical Stack
Numerical Methods: ODE Integration (RK4), Stochastic Sampling.

Physics Engines: Collision detection and Contact Force modeling.

Environment: [MATLAB]

🏛 Academic Disclaimer
This repository is curated as a digital portfolio for coursework completed at Otto-von-Guericke University. It serves as an independent implementation of core process engineering principles and is intended for educational and review purposes.
