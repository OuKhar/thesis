14.4.2020, Oskari Jakonen
---------------------------------------------------------------------------------
# instructions 
on how to run economic models with dynare on Oskari Jakonen Master's Thesis 
---------------------------------------------------------------------------------
* the models run with dynare (https://www.dynare.org/) on matlab (https://se.mathworks.com/?s_tid=gn_logo)
* the latest versions of the scripts were updated were
	- dynare: 4.6.1
	- matlab: R2020a
---------------------------------------------------------------------------------
# run the models
---------------------------------------------------------------------------------
* after obtaining matlab, and dynare and integrating dynare path to matlab, run the models here with:
* (default run argument is "oskari")
"run gk"  
- runs the base model with six-agent economy and financial intermediary being the banking sector
- performs basic checks
 "run mazelis"
- runs the model with shadow banking added to financial sector
- performs basic checks
"run oskari"
- runs the model with shadow banking and with altered monetary policy and certain parameters
- performs basic checks
"run build"
- runs a build with shadow banking models (oskari and mazelis) and draws IRF's used in the thesis
- for the build to have compatible models, the shocks-section of the .mod files have to identical! 
---------------------------------------------------------------------------------
# more information
---------------------------------------------------------------------------------
* for more information on the specific models, check: 
- "Modelling Shadow Banking in China, a DSGE Approach", chapters 4 and 5
---------------------------------------------------------------------------------
# modify
---------------------------------------------------------------------------------
* to modify the models (shocks or observed variables, modify the .mod files in /dynare_files
    * N.B.! The shocks and variables in the build run are defined in the run-script!
---------------------------------------------------------------------------------
