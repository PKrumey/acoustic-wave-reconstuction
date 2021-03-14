function Func = fSeries(sinCoef,cosCoef,halfL)
% FSERIES Calculates the Fourier Series out of sinCoef and cosCoef
x = -halfL:halfL;
k1 = 0:(length(sinCoef)-1);
k2 = 1:(length(cosCoef)-1);
Func = sum(sinCoef'.*sin(pi*x.*k1'/halfL))+sum(cosCoef(2:end)'.*cos(pi*x.*k2'/halfL))+0.5*cosCoef(1);
end

