# Thoracic-Electrical Impedance Tomography-Imaging-Algorithms
This MATLAB package include code and data for the execution of several difference imaging algorithms applied on dynamic thoracic imaging. It also includes quantitative metrics (Figures of Merit-FoM) and comparison between state-of-the-art approaches and the recently proposed Method-of-Moment based Sparse Bayesian Learning one. 

Requires the EIDORS package libary (version 3.9 or later) to be executed. https://eidors3d.sourceforge.net/download.shtml
The in-vivo data used is included in the EIDORS library online. 

Testbench.m script file: 
Performs EIT image reconstruction based on dynamic (5 breathing states/variant boundaries and admittances) 3D simulated thoracic structures and in-vivo data. Algorithms used: Iterative Tikhnonov Regularization, Total Variation Regularization, Absolute EIT imaging (difference of reconstructed conductivities), GREIT, Non-Linear Difference EIT imaging, Linear EIT imaging with movement prior and Method-of-Moment radial basis function based EIT imaging. 

BSBL_Method_of_Moment.m script file:
Performs EIT image reconstruction based on the above structures and the recently proposed Sparse-Bayesian Learning EIT combined with a Method-of-Moment radial basis function approach.

FOM_testbench.m script file:
Extracts the quantitative metrics of the imaging results based on the above reconstruction approaches. Including GREIT metrics, Pearson Correlation Coefficient and RMSE.

Please Cite:

---Dimas, C., Alimisis, V., Uzunoglu, N., & Sotiriadis, P. P. (2021). A point-matching method of moment with sparse Bayesian learning applied and evaluated in dynamic lung electrical impedance tomography. Bioengineering, 8(12), 191.

---Dimas, C., Uzunoglu, N., & Sotiriadis, P. P. (2021). An efficient Point-Matching Method-of-Moments for 2D and 3D Electrical Impedance Tomography Using Radial Basis functions. IEEE Transactions on Biomedical Engineering, 69(2), 783-794.

---Liu, S., Jia, J., Zhang, Y. D., & Yang, Y. (2018). Image reconstruction in electrical impedance tomography based on structure-aware sparse Bayesian learning. IEEE transactions on medical imaging, 37(9), 2090-2102.





