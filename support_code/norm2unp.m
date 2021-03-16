function RCnorm = norm2unp(RC)
% ************************************************************************
%%Norm2unp This function normalizes the RC matrix to the maximum of the unpumped RC,
% which is supposed to be in the first column.
% ########################################################################
% INPUTS
% ########################################################################
% RC:   rocking curves to be normalized
%#########################################################################
% OUTPUT
% ########################################################################
% RCnorm: normalized rocking curves
% ########################################################################
% Philipp Krumey
% 14-03-2021
% University of Duisburg-Essen
% ************************************************************************

%calculate the maximum of the unpumped rocking curve
maxi = max(RC);
maxOFF = maxi(1);

%norm all rocking curves to the maximum of the unpumped rocking curve
RCnorm = RC / maxOFF;
end

