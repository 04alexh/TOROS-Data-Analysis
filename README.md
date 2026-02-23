# TOROS-Data-Analysis
This function will process the data from many TOROS fields and create a variety of graphical and analytical products.

## What code is doing
After star detections from all desired nights have been aligned, this program will read through whatever folder contains the photometry data of all stars found in the frame and do a variety of analysis: \
*1.* Generate light curves for EVERY detected star present among all fields. The program will do a median-based systematic removal per night for each star, as well as do offset corrections if there are multiple nights of data as well. Once created, the user can also save all light curve csvs to disk if desired. \
*2.* Option to generate light curve plots for stars if desired as well. \
*3.* Option to generate a phased light curve for a specific star if a period of variability and StarID is known. \
*4.* Generate J Stetson values for all stars in field and a histogram of J values. \
*5.* Generate an error analysis of all stars, plotting RMS vs. magnitude with stars flagged as variables from either Δ90 or J values marked on plot. \
*6.* Option to run tests to determine the most optimal number of observations and the most optimal time interval of observation for a known star.

##Use
This section will explain the function parameters. (**function_parameter** [data type]: Explanation of parameter.) \
\
\
**night_folders** [str]: Folder containing field data seperated into folders per night. Set to NULL by default. \
\
**preload_data** [bool]: If true, the program will not read new data from night folders and will instead read in already loaded data (essentially, this skips the first three subfunctions which load in the data from night_folders and then create light curves and perform systematic removal). Set to F by default. \
\
**preloaded_data** [str]: If preload_data is true, this is the data that gets preloaded into the program. Set to "" by default. \
\
**target** [int]: The star ID of the star you want to do any type of anaylsis to. If this ID does not exist in all observations, an error will be raised. Set to 0 by default. \
\
**write_lcs** [bool]: If true, the program will write the csvs of the light curves to disk. Set to F by default. \
\
**lccsv_folder** [str]: If write_lcs is true, this is the directory where LC csvs will be saved. Set to "" by default. \
\
**run_Var_Analysis** [bool]: If run_Var_Analysis is true, the program will plot and return brightness variability analysis. \
\
**run_RMS_Analysis** [bool]: If run_RMS_Analysis is true, the program will plot and return RMS data about the dataset. \
\
**plot_lcs** [bool]: If true, the program will write light curve plots to disk. Set to F by default. \
\
**lc_folder** [str]: If plot_lcs is true, this is the directory where LC plots will be saved. Set to "" by default. \
\
**field_name** [str]: The name of the field being observed. Only used for aesthetic purposes when labeling plots. Set to "Observed Field" by default. \
\
**plot_lc** [bool]: If true, the program will plot the light curve and systematic of the target. Set to F by default \
\
**fit_data** [bool]: If true, program will perform a simple sinusoid fitting to the target light curve. Set to F by default. \
\
**plot_LSP** [bool]: If true, program will perform a Lomb-Scargle periodogram on the target light curve. Set to F by default. \
\
**plc_use_fitted_P** [bool]: If true, program will use the period from the fitting function to the PLC. Set to F by default. \
\
**plc_use_LSP_P** [bool]: If true, program will use period(s) from the LSP to the PLC. Set to F by default. \
\
**plc_LSP_per_num** [int]: The number of periods you want to be plotted on the phase curve from the LSP. Set to 1 by default. \
\
**plot_plcs** [bool]: If true, the program will plot a phased light curve of a specified star. Set to F by default. \
\
**plc_period** [float]: The period of the star (in days) you want to plot the phased light curve of. Set to 0.0 by default. \
\
**run_obs_pred** [bool]: If true, the program will run the TOROS Observation Number optimization program. Set to F by default. \
\
**run_time_pred** [bool]: If true, the program will run the TOROS Observation Time Interval optimization program. Set to F by default. \
\
**run_LSST_test** [bool]: If true, program will run observation and interval test for the target star. Set to F by default. \
\
**LSST_target_name** [str]: The actual name of the star you want to run the optimization programs on (purely an aethetic thing for plot titles). Set to "" by default. \
\
**LSST_N** [int]: The number of iterations per sample size for the Observation Number program. Set to 0 by default. \
\
**LSST_sigma** [int]: The significance level above which a star is flagged as a variable for the Optimization Number program. Set to 0 by default. \
\
**LSST_M** [int]: The number of iterations for the Observation Time Interval program. Set to 0 by default. \
\
