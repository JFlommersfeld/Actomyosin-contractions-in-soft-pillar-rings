# Actomyosin contractions in soft pillar rings

This is the code used to solve an analytical model that describes the contraction dynamics of actomyosin networks that are attached to a ring of ten soft elastic pillars. The models are defined by a set of coupled differnetial equations that are solved numerically.

**Reference:** https://www.researchsquare.com/article/rs-863696/v1

**Author:** Johannes Flommersfeld

**Contact:** J.Flommersfeld@physik.uni-muenchen.de


Developed in Python 3.7. Dependencies:

    - NumPy, SciPy
    - Optional: Matplotlib (for plotting)

### Contents:

**calculate_contraction_dynamics.py:** The front end of the code that can be used to solve the contraction dynamics of different models.

**plotting.py:** Contains various plotting routines that might be useful to visualize the results of the calculate_contraction_dynamics method.

**models.py:** Contains the two different models for the contraction dynamics implemented here. The 'density model' only considers a possible density dependence of the actomyosin contractility, while the 'full model' also accounts for myosin filament binding dynamics.

**work_and_power.py:** Contains different methods to calculate the dissipated and transmitted power and work from the results of the calculate_contraction_dynamics method.  

**auxiliaries.py:** Contains different auxiliary methods that can be used to account for the descreteness of experimental data and thus allow for a fairer comparison between theory and experiments.

**parameter_loader.py:** Reads in the parameter files.

**full_model_demo.py:** an instructive example of how to use and analyse the 'full model'. 

**density_model_demo.py:** an instructive example of how to use and analyse the 'density model'.
 
