%START routine to starts optimization algorithm
%   This is the start script for running a optimization algorithm to retrieve the
%   strain pulse out of the measured rocking curves for different time
%   delays.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%add paths from folder
addpath(genpath('code'));
addpath(genpath('support_code'));
addpath(genpath('workspace'));

main('workspaceTiGaAs111Full.mat',200,[-0.6 -0.1 0.1 0.6],'none');





