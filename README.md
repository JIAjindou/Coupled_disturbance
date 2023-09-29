# Coupled_disturbance
### Paper 

Estimating Coupled Disturbance via Variable Separation

### Author

Jindou Jia, Yuhang Liu, Kexin Guo, Xiang Yu, Lihua Xie and Lei Guo

### Contant email

jdjia@buaa.edu.cn, kxguo@buaa.edu.cn

### Updates

- September 29, 2023- First release.

### Description

This shared folder contains some supplementary codes and data to support the paper ' Estimating Coupled Disturbance via Variable Separation'. 

1. Simulation.slx

   The simulation project of Section V.B, running on Matlab R2019a.

2. Main_training.m

   The offline training program of Section V.B, running on Matlab R2019a. It will call sub-programs: B_X_fun.m, Cheby_poly.m, xi_fun.m.

3. Main_meshplot.m

   Plot the true and learning surface diagrams of the chosen three nonlinears functions in Section V.A. It will call sub-programs: B_X_fun.m, Cheby_poly.m, xi_fun.m.

4. Main_heatmap.m

   Learning errors of the learning algorithm under different noise variances ${\sigma}_x^2$ and parameters $p$ on the test dataset in Section V.A, i.e., the program of Fig. 2 B.  It will call sub-programs: B_X_fun.m, Cheby_poly.m, xi_fun.m.
