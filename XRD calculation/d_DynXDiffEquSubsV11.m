function dX_Subs = d_DynXDiffEquSubsV11(A, X, gamma_zero, gamma_H, b, g, k, y_a, y_c, BraggAngle, A_T, S_z)
% ************************************************************************
% this function provides the equation used in differential equation for
% Dynamic X-Ray diffraction for the substrate
% ########################################################################
% INPUTS
% ########################################################################
% A:    the reduced spatial coordinate of the substrate starting from the
%       interface
% X:    the complex scattering amplitude (this is the solution of the
%       equation)
% gamma_zero, gamma_H, b, g, k, y_a, y_c, BraggAngle, A_T:
%       coefficients in the DXRD (Larsson Paper)
% S_z:  strain Vector for the given time delay; size(S_z) = 2000 1
%#########################################################################
% OUTPUT
% ########################################################################
% dX_Subs: provides differential equation for DXRD
% ########################################################################
% Mohammadmahdi Afshari; edited by Fabian Brinks
% 02-02-2018, University of Duisburg-Essen
% % changelog 21.11.2019:
%       - included gamma_zero and gamma_H as new inputs to remove as much
%       calculations as possible from the function
% *************************************************************************

% calculate depth from A
depth = A_T.*A*sqrt(abs(gamma_zero.*gamma_H)); 

% choose right value of S_z for given 
if (floor(depth*10^9) + 1) > length(S_z)
    SB = S_z(length(S_z));
else
    SB = S_z(floor(depth*10^9)+1);
end

% calculate alphaH and y
alphaH = 4*sin(BraggAngle*pi/180)*( sin(BraggAngle*pi/180) + gamma_zero - gamma_H.*SB );
y = y_a * (b.*alphaH)./sqrt(abs(b)) +  y_c*((b-1)./sqrt(abs(b)));

% calculate the differential equation
dX_Subs = [k*(X(1).^2 - X(2).^2 +1) + 2*X(2).*( X(1) - y)-2*g*X(1); ...
    -(X(1).^2 - X(2).^2 +1) + 2*X(1).*(X(2)*k + y) - 2*g*X(2)];
end
