## Reproduction of statistical analyses:

Raphael Eisenhofer, Iñaki Odriozola and Antton Alberdi

The statistical analyses of the manuscript "Impact of microbial genome completeness on functional metagenomics" can be reproduced by consecutively running the scripts titled Step1 through Step4.

- **"Step1_Load_and_Prepare_Data.R"**

Loads the raw datasets in the /data/mag_data/DRAM/ folder of this repository, and combines and reformats them to use as imput for modelling in the Step2 script.

- **"Step2_Fit_Fullness-Completeness_GLMMs"**

Loads the datasets outputted by Step1 script and fits binomial GLMs to estimate the fullness-completeness association for each module and bacterial phylum. Then, it outputs the parameter estimates of the models as well as average module fullness predictions for different completeness values.

- **"Step3_Explore_GLMM_Parameters_and_Predictions"**

Loads the datasets outputted in Step2 and computes summary statistics of the binomial GLMs associating module fullness and MAG completeness. This script also generates Figure 2A in the main manuscript.

- **"Step4_Predictors_Fullness-Completness_strength_LMM"**

Loads the data outputed in Step2, reformats it to fit the LMM associating the slopes estimated in the binomial GLMs with the phyla of the MAGs, domains and number of steps of the modules. Then it makes inferences of the marginal influences of MAG Phylum, module domain and number of steps on the strength of the fullness-completeness associations. This script generates the Figures 2B, 2C and 2D in the main manuscript. 
