function Coefs = CosCoefs(f,maxk,halfL)
%SINCOEFS Calculates maxk Cos Fourier Coefficients of f
x = -halfL:halfL;
k = 0:maxk;
Coefs = trapz(x,f'.*cos(pi*k.*x'/halfL)/halfL);
end

