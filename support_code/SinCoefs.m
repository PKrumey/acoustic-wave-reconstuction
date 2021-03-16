function Coefs = SinCoefs(f,maxk,halfL)
% ************************************************************************
%SINCOEFS Calculates a number of sine fourier coefficients for the 
%         function f
% ########################################################################
% INPUTS
% ########################################################################
% f:      a function of the length 2*halfL+1
% maxk:   the number of sine coefficients to be calculated plus 1:
%         0:maxk (sine coefficient for k=0 is always 0)
% halfL:  half the length of the function f
%#########################################################################
% OUTPUT
% ########################################################################
% Coefs:  calculated sine coefficients
% ########################################################################
% Philipp Krumey
% 14-03-2021
% University of Duisburg-Essen
% ************************************************************************
x = -halfL:halfL;
k = 0:maxk;
Coefs = trapz(x,f'.*sin(pi*k.*x'/halfL)/halfL);
end

