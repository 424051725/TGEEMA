Matlab code for "Longitudinal Tensor Data Learning via Optimal Model Averaging"

DATE:  Mar 2024

Tested on MATLAB 2017a and 2024a
-------------------------------------------------------------------------------

The following MATALB script files are used to estimate the TGEEMA model on the datasets studied in the main paper (in Section 4)

The main function for simulation is 'simulation_demo.m'

To load the related toolbox, click "Set Path" in Matlab and choose "add with subfolders" for the folder "TensorGEE".

* simulation_demo.m
Demo script, generate 2-D synthetic dataset and estimate the TGEEMA model. This script will generate 2-D data with a sample size of n=200 and m=4, then fit the CP tensor GEE model of rank 1-5, and calculate the optimal model averaging weight. The code will eventually output the restored image of each candidate model and the restored image of the TGEEMA method. Finally, the rmse of each model on the testing data and the RMSE of parameter prediction will be recorded in the "2d_simulation" folder.

-------------------------------------------------------------------------------

Main functions:

* data_2d.m
Generate 2-D simulation data

* model_averaging.m
Calculate CP tensor GEE model from rank-1 to rank-maxrank and different woeking correlation matrixes. Finally calculate optimal model averaging weight

*obj_func.m
The square loss function for optimization

*simulation_func.m
The main function part, including calculate the optimal weight, the estimator of each candidate model, and save the final results and recovered images

*split_data.m
Split the ehole samples into train and test data set

*tensor_gee_ma.m
The adaptive tensor GEE estimation function

*test_simulation.m
Calculate RMSE of y and B on the testing data
-------------------------------------------------------------------------------
