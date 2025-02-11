# Digital Twin Calibration using Point Clouds

# Introduction

This study presents a new approach for the automated system identification and model calibration of Digital Twins (DTs) for prefabricated spaceframe structures subjected to static displacement using point clouds. The methodology includes: (i) extracting prefabricated structural elements from point clouds using artificial intelligence (AI)-based robust model fitting; (ii) estimating nodal coordinates, positional uncertainties, and relative displacements from temporal point clouds; (iii) modeling the unknowns in the governing equations as decision variables; (iv) recovering the decision variables through AI-based minimization of the difference in relative nodal displacements between the DT prediction and the real-world measurements; and (v) accounting for the uncertainties in the measured nodal coordinates through relaxation strategies used in robust optimization. 

The following image shows the case study used for validation of the methodology. The controlled static vertical displacement before and after the removal of the temporary bracing (support) is estimated to calibrate the predicted displacements in the finite element model (FEM).


<img width="600" alt="image" src="https://github.com/user-attachments/assets/37a6dfa2-682f-4399-a8bb-f913b498b17a" />


# Algorithm and Optimization Schematics

The code employs a metaheuristic AI-based algorithm to recover the joint fixity factors of vertical elements in the longitudinal direction as well as the loading coefficients to account for load case uncertainty. The AI-based algorithm aims to find the decision variables that minimize the difference between the measured (from laser scanner) and predicted (from FEM simulation) of the vertical displacements for all nodes, while incorporating the uncertainty in the estimated real-world nodal coordinates. The following image shows the process of employing the AI-based optimization to maximize agreement between the DT and reality (the image on the left shows the real vertical displacement with 10x magnification for clarity).


<img width="710" alt="image" src="https://github.com/user-attachments/assets/067db487-ed1c-4a7e-89ab-6ff2da7ff2f7" />


# Results of the Study

The code provides the user with the opportunity to perform the optimization using three metaheuristic algorithms, namely, Particle Swarm (PS), Genetic Algorithm (GA), and Simulated Annealing (SA). The code is implemented in Matlab and requires the Parallel Computing Toolbox as well as the Optimization Toolbox. The code also includes a script to graph the results of the metaheuristic algorithms, which run for a pre-specified number of simulations (here, the default is set to 50). The following image is a sample generated from the "General Plotting" script after the results of the "Main Code" are obtained.


<img width="566" alt="image" src="https://github.com/user-attachments/assets/054aa2b8-3b25-4bae-ad54-5762d806dc53" />


# Explanation on Code and Provided Files


As part of this Repository, three folders are provided:

1- Inputs, which includes the required input data, such as real and FEM displacements.

2- Dependent Functions, which includes the new functions called in the main script.

3- Demonstration, which includes the main script along with the code for the general plotting.

To successfully run the code, all items must be accessible to Matlab at the time of execution.


# Citation
The study is expected to be published in the proceedings of the ISPRS 6th Joint International Conference on Deformation Monitoring (JICDM), where the mathematical formulations are explained in more detail. You may cite the study using the following information:

Maalek, R. 2025. Digital Twin Model Calibration of Structural Systems under Visible Deformation: Recovering the Joint Fixity of a Steel Cantilever Spaceframe Structure Experiencing Static Displacement using Laser Scanning, Proceedings of the ISPRS 6th Joint International Conference on Deformation Monitoring (JICDM), Karlsruhe, Germany, April 2025.

# Acknowledgements
The author wishes to acknowledge the generous endowment provided by GOLDBECK GmbH to the Karlsruhe Institute of Technology (KIT) for the establishment of the Professorship in Digital Engineering and Construction at the Institute of Technology and Management in Construction (TMB).
