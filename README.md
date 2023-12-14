# GA-PLSR
There are mainly three '.m' file you need to modify to run the code: **demo_general_train_V4_2.m Beta_integrate_stop_GA_intime_V3_Oxide.m and Beta_R2_Setsplit_inpro_plot_V4.m**

There is one data file: **train_dat.mat**, that helps you run this demo.

## main function
**demo_general_train_V4_2.m**

This is the main function that trains the GA-PLSR model. Just run the code after setting propers file directories in the very first part of the code.

In this demo, you may also want to change the parameters in the "model parameters" part. Make sure you have understood these parameters before changing them. You can also just change the parameter value and see what will happen.

## postprocessing
**Beta_integrate_stop_GA_intime_V3_Oxide.m**
This function processes the data produced by the main function. It integrates all the solutions trained by GA-PLSR, selectes the top 50, and shows the results.

**Beta_R2_Setsplit_inpro_plot_V4.m**
This function visualizes the final solution of the GA-PLSR model, and evaluates the uncertainties.
