function Coefs = SinCoefs(f,maxk,halfL)
%SINCOEFS Calculates maxk Sin Fourier Coefficients of f
x = -halfL:halfL;
k = 0:maxk;
Coefs = trapz(x,f'.*sin(pi*k.*x'/halfL)/halfL);
end

