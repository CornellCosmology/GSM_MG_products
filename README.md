# GSM_MG_products 

This repository contains codes for the calculation of the anisotropic redshift-space correlation function for biased tracers in MG theories using the Gaussian Streaming Model (GSM), as published in

G.Valogiannis, R. Bean and A. Aviles (2020), https://arxiv.org/abs/1909.05261 , JCAP01(2020)055 

The real-space ingredients that enter the GSM are evaulated using the Convolution Lagrangian Perturbation Theory (CLPT) approach for MG, as it was presented in 

G.Valogiannis and R. Bean (2019), https://arxiv.org/abs/1901.03763

the codes for which have been separately released publicly in https://github.com/CornellCosmology/bias_MG_LPT_products . 

The codes are the MG expansions of the original GR versions developed in https://github.com/martinjameswhite/CLEFT_GSM , for 
calculating the various ingredients of the Gaussian Streaming Model (GSM) using EFT corrections to CLPT.
Note, however, that the MG version presented here Only produces predictions for biased tracers using a local-in-density Lagrangian bias, with an EFT-inspired correction only to the velocity dispersion.
The Kernels and subroutines that refer to the tidal bias terms and the rest of the EFT corrections exist as part of the original code and have NOT been adapted for MG, even though the some of the corresponding expressions are derived and presented in https://arxiv.org/abs/1909.05261 . 

Input description:

For each calculation, the user needs to provide 4 files: 
1) The MG linear power spectrum for the particular cosmology and redshift
2) The MG linear growth rate for the same model
3) The k-dependent Q_n(k) and R_n(k) functions in MG. For more detailed instruction of how to obtain predictions for these, users are referred to the python modules released in https://github.com/CornellCosmology/bias_MG_LPT_products
4) The k-dependent Q^f_n(k), R^f_n(k),Q^ff_n(k), R^ff_n(k) functions in MG, as defined in Appendix A2 of G.Valogiannis, R. Bean and A. Aviles (2020). These can be evaluated through straightforward expansions of the Mathematica and Python codes released previously in  https://github.com/CornellCosmology/bias_MG_LPT_products. Files with these functions are provided here, and a more detailed discussion of how to produce them, together with the corresponding codes, is under construction.

The code reads (2)-(4) in the subroutine setupQR of lsm.cpp. 

How to compile and run:

After compiling with 'make' in the main directory, the executable 'lesm' is produced, which can be called as:
./lesm <Pk-file> <ff> <b1> <b2> <bs2> <Aeft> <Aeftv> <s2FoG> <ngrav> 
where
<Pk-file> is MG linear power spectrum from 1)
<ff> is the GR linear growth rate for the given z
<b1> is the linear Lagrangian bias
<b2> is the 2nd order Lagrangian bias
<bs2> is the tidal bias -> currently not used set = 0
<Aeft>, <Aeftv> are EFT corrections -> currently not used set = 0
<s2FoG> is the EFT-insipred Fingers-Of-God offset corrections to the velocity dispersion, and 
<ngrav> the gravity model, =0 for GR and =1 for MG.

Output description:

The code produces a file named 'xi0GSM_<Pk-file>', which contains the monopole, \xi_0, quadrupole, \xi_2 and hexadecapole \xi_4 of the RSD correlation function at a given separation s, in the column following format:

s     s^2*\xi_0  s^2*\xi_2   s^2*\xi_4

As an example, for the Fr6 model at z=0.5, which is one of the results shown in G.Valogiannis, R. Bean and A. Aviles (2020), 
we can run 
./lesm plin_Fr6z05wmap9.txt 0.77 0.23 -0.49 0 0 0 -15.5 1

which reproduces the results of the paper in the file 'xi0GSM_plin_Fr6z05wmap9.txt'.

The input files for this evaluation are provided, with the following names
1) plin_Fr6z05wmap9.txt
2) fgrowth_Fr6z05wmap9.txt
3) plin_Fr6z05wmap9_Q.txt, plin_Fr6z05wmap9_R.txt
4) plin_Fr6z05wmap9_QDer.txt, plin_Fr6z05wmap9_RDer.txt

The inputs and outputs for the rest of the models in G.Valogiannis, R. Bean and A. Aviles (2020) are also provided.





