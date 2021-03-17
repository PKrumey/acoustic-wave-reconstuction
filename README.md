# acoustic-wave-reconstuction
Matlab optimization algorithm to reconstruct acoustic waves from time-resolved X-Ray diffraction data

# The Folder structure is as followed:
  
reconstruction: holds main programs for the optimization algorithm  
-constants: data class that contaigns all the fixed properties of the strain pulses and the substrate  
-DNA: data class that contains all data about the bipolar pulses  
-main: function to start the optimization algorithm  
-population: class where the reconstruction takes place  
-firstGeneration: calculates the first generation    
-calcCRC: function that calculates the rocking corves to a bipolar pulse  
-calcFitness: function that calculates the fitness between the measured and the calculated rocking curve after filtering
  
support_code: holds the subroutines for the reconstruction   
-CosCoefs: function that calculates the cosine coefficients out of a bipolar pulse  
-SinCoefs: function that calculates the sine coefficients out of a bipolar pulse  
-fSeries: function that calculates the fourier series out of sine and cosine coefficients  
-DataHash: function that gives every DNAObject a unique ID  
-logNorm: function that calculates the fitness between the measured and the calculated rocking curve  
-norm2unp: function that norms the rocking curves to a unpumped rocking curve  
-Thomsen: calulate a thomsen pulse with random parameters  
   
XRD calculation: holds the data and programs to calculate and convolute the rocking curves  
-asymmetric_lorentz_convolution_v5_GaAs111: convolution between the measured rocking curves, an asymmetric lorentzian function and some constant vlaues  
-d_DynXDiffEquSubsV11: contains the differential equation used for calculating the rocking curves  
-DXRD_ode23_C_mex.mexw: pre compilated function that solves the differential equation for calculating the rocking curves  
-RC_Ti_GaAs: File that contains the averared unpumped rocking curves for calculating the instrument function  
-substrate_convolution_v5_GaAs111sm: function that calculates the convolution of rocking curves using a asymmetric loretzian function and the averaged unpumped rocking curves  
-XRD_inputsV8_m_s: structure calculating the parameters used in the differential equation out of the constants of the substrate  
   
workspace: holds the workspace with the measured or the synthetic data in the following format   
-time:  vector contaigning the time points of the rocking curves  
-theta: vector contaigning the angles of the rocking curves  
-MRC: array containing the rocking curves in the form [theta,times]. first column has to be a unpumped rocking curve  
(-pulse: vector that contains a synthetic bipolare pulse (optional))   


# Start Instructions:
run the "start"-script to to call the main function with the parameters: (workspace, generation size, filter area for fitness calculation, precalculated first generation (optional).

# Constraints
  
number of fourier coefficients: the number of fourier coefficients included in the reconstruction can be found in firstGeneration: line 52, line 54   
random mean strain: the random mean strain for the bipolar pulses of the first generation can be set in firstGeneration: line 60   
fourier coefficient area: the area around the fourier coefficients that is searched for optimization can be found in population: line 132
Reflectivity and bipolar pulse length: the reflectivity and the bipolar pulse length can be set in constants: 



