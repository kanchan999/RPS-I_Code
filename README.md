# Regenerative Population Strategy-I (RPS-I)

This repository contains Python implementations of RPS-I for multiple metaheuristic algorithms, illustrating how RPS-I dynamically addresses structural bias during optimization. 

# The code includes:
6 files (RPS_I_GA, RPS_I_PSO, RPS_I_DE, RPS_I_GWO, RPS_I_WOA, RPS_I_HHO) for comparison of strauctual bias with the base version through Generalised Signature Test.
6 files (RPS_I_GA_POPULATION_PLOT, RPS_I_PSO_POPULATION_PLOT, RPS_I_DE_POPULATION_PLOT, RPS_I_GWO_POPULATION_PLOT, RPS_I_WOA, RPS_I_HHO_POPULATION_PLOT) for visuliazation of popualtion.
Tension/Compression Spring Design Problem (continuous convex benchmark).
Pressure Vessel Design Problem (engineering design benchmark).

# How to Run

Install Dependencies
Python 3.7+
NumPy
Matplotlib (for plotting convergence curves)
Clone or Download the repository.
Navigate to the repository folder, then run either problem’s script (e.g. for the spring design):
python compression_spring_design_problem.py
or

python Pressure_vessel_design_problem.py
Each script will:

Initialize the population
Run all 12 algorithms (standard + RPS-I variants)
Print best solutions and constraint feasibility
Optionally plot a multi-line convergence chart
Check Output
The console/log will show each algorithm’s best solution, penalized fitness, objective value, and whether constraints are satisfied.
A popup window with the convergence plot may appear if plotting is enabled.

# RPS-I Overview

Regenerative Population Strategy-I is a dynamic approach for mitigating structural bias. At each generation:

Measure:
Population diversity (α)
Improvement rate (β)
Compute γ
Reinitialize N indvidulas based on \(N = \psi(\gamma) = \lfloor (1-\gamma) \cdot (P-1) \rfloor\) 
This helps the algorithm avoid premature convergence and maintain better exploration/exploitation trade-offs.

# Citation

If you use or reference this code in your publications, please cite the paper:

Kanchan Rajwar et al., “Regenerative Population Strategy-I: A Dynamic Methodology to Mitigate Structural Bias in Metaheuristic Algorithms.”

# License

This code is provided for academic and research purposes. 

# Contact

For questions or collaboration, feel free to contact:

Author: Kanchan Rajwar
Email: kanchanrajwar1519@gmail..com
