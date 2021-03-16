function Coefs = CosCoefs(f,maxk,halfL)
% ************************************************************************
%SINCOEFS Calculates a number of cosine fourier coefficients for the 
%         function f
% ########################################################################
% INPUTS
% ########################################################################
% f:      a function of the length 2*halfL+1
% maxk:   the number of cosine coefficients to be calculated plus 1:
%         0:maxk
% halfL:  half the length of the function f
%#########################################################################
% OUTPUT
% ########################################################################
% Coefs:  calculated cosine coefficients
% ########################################################################
% Philipp Krumey
% 14-03-2021
% University of Duisburg-Essen
% ************************************************************************
x = -halfL:halfL;                               
k = 0:maxk;
%calculate the cosine coefficients according to the equation
%Coefs_k = 1/halfL*integrate_{-halfL}^{halfL} f*cos(pi*k*x/halfL) dx
Coefs = trapz(x,f'.*cos(pi*k.*x'/halfL)/halfL);
end

