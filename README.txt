This repository contains code to fully reproduce the simulation results and application illustration in Knowlton, R., Parast, L. (2024) ``Efficient Testing Using Surrogate Information," which is under review. For questions or comments, please contact R. Knowlton at rknowlton@utexas.edu.

The code has been written using R version 4.3.3 (platform: x86_64-w64-mingw32/x64, 64-bit) with packages hetsurr_1.0, matrixStats_1.3.0, quantreg_5.97, Rsurrogate_3.2, ggplot2_3.5.1, and etsi_1.0. Note that the ``etsi'' package is available on Github at https://github.com/rebeccaknowlton/etsi (Commit 93e07c7).

Additionally, to reproduce the AIDS examples, you will need to request the data directly from the AIDS Clinical Trial Group. Specifically, you may request the ACTG 320 study and the ACTG 193A study at https://actgnetwork.org/submit-a-proposal.

The following steps should be completed in order.

To reproduce Table 1 in the paper:
Open the file etsi_sims_110424.R and input the desired setting on line 5. (Table 1 includes settings 1-3, and setting 4-5 are included in Table A1 in Appendix C.) The data from Study A (which is fixed for every iteration) and the simulation output (1000 iterations) will be saved in the folder sim_output. Then, open etsi_sims_readin.R and enter the corresponding setting number in line 5, and run the code. The results for Table 1 for each setting will be saved in the folder titled results.

To reproduce Table 2: 
Open the file etsi_design_sims_110424.R and input the desired setting on line 6. (Note that this script relies on functions from etsi_sims_110424.R, so these should still be loaded in the current working session from creating Table 1.) The results for Table 2 for each setting will be saved in the folder titled results.

To reproduce Figure 2 and Table 3: 
Run the code in aids_example.R. You must have the ACTG data available in a suitable working directory. Figure 2 and Table 3 will be saved in the results folder. (The data from Study A will also be saved in the folder aids_output.)

To reproduce Figure 3:
Run the code in aids_design.R. Figure 3 will be saved in the results folder.

Note that the file etsi_sims_truth_and_assumptions.R calculates the truth and checks necessary assumptions for the simulation settings. It is included for completeness, but not required to produce any tables or figures.

