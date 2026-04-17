# ARS_INPUT

This folder includes simulation inputs for the submitted paper "Bayesian optimization of accelerator-radiator scheme for efficient X-ray generation."

template_input_ars.py -- template for Smilei (4.8) simulation, 256 -- 384 CPUs were used to conduct these simulations under 24 h 
 
para_coords.txt -- parametric coordinates acquired from BO, names are in the first row, can be directly copied and pasted to the input script

Chi_tes -- folder containing 2 simulation inputs, carried out to study effects of lower $\chi$ on the the conversion effeciency
        -- they set to save checkpoint after 46 h of computational time, 2-3 reruns are required
