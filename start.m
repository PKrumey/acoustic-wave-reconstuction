%START routine to starts optimization algorithm
%   This is the start script for running a optimization algorithm to retrieve the
%   strain pulse out of the measured rocking curves for different time
%   delays. Here the main function of the algorithm is called
% ########################################################################
% Philipp Krumey
% 14-03-2021, University of Duisburg-Essen
% *************************************************************************

%add paths from folder
addpath(genpath('reconstruction'));
addpath(genpath('support_code'));
addpath(genpath('workspace'));
addpath(genpath('XRD calculation'));

%start main routine with (rocking curves, generation size, filter, precalculated first generation)
main('workspaceTiGaAs111Full.mat',200,[-0.6 -0.1 0.1 0.6],'none');





