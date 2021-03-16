function Func = fSeries(sinCoef,cosCoef,halfL)
% ************************************************************************
% FSERIES Calculates the fourier series from the sine coefficients sinCoef
%         and the cosine coefficients cosCoef
% ########################################################################
% INPUTS
% ########################################################################
% sinCoef:   sine coefficients
% cosCoef:   cosine coefficients
% halfL:     half the length of the function Func
%#########################################################################
% OUTPUT
% ########################################################################
% Func:      the resulting fourier series
% ########################################################################
% Philipp Krumey
% 14-03-2021
% University of Duisburg-Essen
% ************************************************************************
x = -halfL:halfL;
k1 = 0:(length(sinCoef)-1);       %number of sine coefficients
k2 = 1:(length(cosCoef)-1);       %number of cosine coefficients
%calculate the fourier series
Func = sum(sinCoef'.*sin(pi*x.*k1'/halfL))+sum(cosCoef(2:end)'.*cos(pi*x.*k2'/halfL))+0.5*cosCoef(1);
end

