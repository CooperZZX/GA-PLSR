# GA-PLSR
There are mainly three '.m' file you need to modify to run the code: **demo_general_train_V4_2.m Beta_integrate_stop_GA_intime_V3_Oxide.m and Beta_R2_Setsplit_inpro_plot_V4.m**

There is one data file: **train_dat.mat**, that helps you run this demo.

## main function
**demo_general_train_V4_2.m**

This is the main function that trains the GA-PLSR model. Just run the code after setting proper file directories in the very first part of the code.

In this demo, you may also want to change the parameters in the "model parameters" part. Make sure you have understood these parameters before changing them. You can also just change the values and see what will happen.

## postprocessing
**Beta_integrate_stop_GA_intime_V3_Oxide.m**
This function processes the data produced by the main function. It integrates all the solutions trained by GA-PLSR, selectes the top 50, and shows the results.

**Beta_R2_Setsplit_inpro_plot_V4.m**
This function visualizes the final solution of the GA-PLSR model, and evaluates the uncertainties.

# Description
## Training of a GA-PLSR model

<img src="https://github.com/CooperZZX/GA-PLSR/blob/main/images/GA-PLSR%20architecture.jpg" width="45%" height="45%">

**Architecture of the GA-PLSR model**

Here, m, n, and b are the number of chromosomes (band combinations) before selection, the number of chromosomes after selection, and the number of spectral bands, respectively. Each chromosome carries binary-coded genes that toggle the spectral bands on or off. The relationship between m and n is defined as n = m×gap, where gap (generation gap) represents the proportion to generate the next generation. The chromosomes exchange genes during the cross step and acquire new genes during the mutation step. The RMSEP is calculated using the validation set. The SD is the standard deviation of the predicted oxide contents from CE-4 D10 spectra and represents the effects of viewing geometry on the model. The fitness is obtained by adding of RMSEP and SD, effectively representing the overall quality of a chromosome.

<img src="https://github.com/CooperZZX/GA-PLSR/blob/main/images/Flowchart.jpg" width="60%" height="60%">

**The flowchart of the code**

(a) The preprocessing of training data. The train_dat.mat is the result of this step.

(b) The training process of GA-PLSR. Model parameters of every epoch are stored together for the next selection. The code was run 120 times in 6 different generational differences with 20 solutions each time to avoid local optima. STD is the standard deviation of each oxide’s contents estimated from CE-4 D10 spectra. RMSEP is calculated from the validation set. K is the efficient between RMSEP and STD. 

(c) Eliminate overfitting and select the top 50 models with the smallest RMSEP from 2400 (120×20) solutions for each oxide. R2(n) is the R-square of the model on validation set at the nth epoch. THR is threshold.

<img src="https://github.com/CooperZZX/GA-PLSR/blob/main/images/Train_process.jpg" width="60%" height="60%">

**The training process of GA-PLSR model (take FeO for an example)**

(a) Split of training dataset and validation dataset of reflectance spectra. About 2/3 spectra are selected as the training dataset and the other 1/3 spectra become the validation dataset. Both of them are evenly selected from the reference data based on the oxide contents. 

(b) Split of training dataset and validation dataset of the 1st derived spectra. The vertical line marks the stop time of the training process when the rate of fitness decrease slowed. 

(c) Changing of root mean square error of prediction (RMSEP), standard deviation (STD), and fitness during iterations. 

(d) Changing of R-square (R2) during iterations.
