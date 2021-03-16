function fitness = logNorm(MRC, CRC)
% ************************************************************************
%logNorm: calculates the fitness between two rocking curves by taking
%         the sum over the squared deveation of the logarithmic rocking
%         curves
% ########################################################################
% INPUTS
% ########################################################################
% MRC:          The measured rocking curve [angles,times]
% CRC:          The calculated rocking curvs [angles,times]
%#########################################################################
% OUTPUT
% ########################################################################
% fitness:      a single value > 0 expressing the fiteness between MRC
%               and CRC
% ########################################################################
% Philipp Krumey
% 14-03-2021
% University of Duisburg-Essen
% ************************************************************************

        %calculate logarithm of MRC, and CRC
        normMRC = log(MRC+1e-4);
        normCRC = log(CRC+1e-4);
        
        %determine the number of angles and time points for normalization
        [thetas,times] = size(CRC);
        
        % fitness is sum of normed deviation squares
        fitness = sum(sum((normMRC - normCRC).^2)/thetas)/times;
        
end

