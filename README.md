# Regenerative Population Strategy-I (RPS-I)

This repository contains Python implementations of RPS-I for multiple metaheuristic algorithms, illustrating how RPS-I dynamically addresses structural bias during optimization. 

# The code includes:
- 6 files (e.g. RPS_I_GA) for comparison of strauctual bias with the base version through Generalised Signature Test.
- 6 files (e.g. RPS_I_GA_POPULATION_PLOT) for visuliazation of popualtion.
- Tension/Compression Spring Design Problem (continuous convex benchmark).
- Pressure Vessel Design Problem (engineering design benchmark).

# How to Run

- Install Dependencies: Python 3.7+, NumPy, Matplotlib (for plotting convergence curves)
- Clone or Download the repository.
- Navigate to the repository folder, then run either problem’s script (e.g. for the spring design):

# RPS-I Overview

Regenerative Population Strategy-I is a dynamic approach for mitigating structural bias. At each generation:

Measure:
- Population diversity (α)
- Improvement rate (β)
- Compute γ
- Reinitialize N indvidulas based on Eq. (14)
This helps the algorithm avoid premature convergence and maintain better exploration/exploitation trade-offs.

# Citation

If you use or reference this code in your publications, please cite the paper:

Kanchan Rajwar et al., “Regenerative Population Strategy-I: A Dynamic Methodology to Mitigate Structural Bias in Metaheuristic Algorithms.”

# License

This code is provided for academic and research purposes. 

# Contact

For questions or collaboration, feel free to contact:

- Author: Kanchan Rajwar
- Email: kanchanrajwar1519@gmail..com
