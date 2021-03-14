function pulse = ExpGauss3(ll)

lambda = (5+1005*rand)*10^-9;
v = 4730;
D=1.2*rand*v*lambda;

A =(-ll:ll)*10^-9;

% Define the vector tprime of the integral to be calculated numerical
tprimemin = .0001e-12; 
tprimemax = 500e-12;
dt = .01e-12;
tprime = tprimemin:dt:tprimemax; Nbrt = numel(tprime);

F1 = zeros(length(A),1);
F2 = zeros(length(A),1);
F = zeros(length(A),1);

% Define the functions phi and x, both for the cases +/- and C/D. The
% functions are defined in Philipp's evaluation of the derivative of sigma
% (photos 1 - 4). The subscripts +/- are defined there.
% The superscripts C and D belong to the two different delta-functions,
% which define where the derivative is evaluated. For C, the functions are
% evaluated at (A+vt), for D they are evaluated at (-A-vt). Those
% evaluation points as well as the integration borders are calculated in
% photos 5 - 6 for the equations in (14) and (15), respectively.
% The evaluation of the phi and x functions for the different cases is done
% in photo 7; In the following script they are just considered as functions
% of A and t
phiminC = @(varA, vart) 1/lambda*((D/lambda-v)*vart-varA);
phipluC = @(varA, vart) 1/lambda*((D/lambda+v)*vart+varA);
xminC = @(varA, vart) -(varA+(v-2*D/lambda)*vart)/sqrt(4*D*vart);
xpluC = @(varA, vart) -(varA+(v+2*D/lambda)*vart)/sqrt(4*D*vart);

phiminD = @(varA, vart) phipluC(varA, vart);
phipluD = @(varA, vart) phiminC(varA, vart);
xminD = @(varA, vart) (varA+(v+2*D/lambda)*vart)/sqrt(4*D*vart);
xpluD = @(varA, vart) (varA+(v-2*D/lambda)*vart)/sqrt(4*D*vart);

% Define the two integrals C and D, which are the derivative of sigma
% evaluated at (A+vt) and (-A-vt), respectively.
IntC = @(varA, vart) D/(2*lambda^2)*(exp(phipluC(varA, vart))*(1+erf(xpluC(varA,vart)))...
                +exp(phiminC(varA,vart))*(1-erf(xminC(varA,vart))))...
                -sqrt(D/(pi*lambda^2*vart))*exp(-(varA+v*vart)^2/(4*D*vart));

IntD = @(varA, vart) D/(2*lambda^2)*(exp(phipluD(varA, vart))*(1+erf(xpluD(varA,vart)))...
                +exp(phiminD(varA,vart))*(1-erf(xminD(varA,vart))))...
                -sqrt(D/(pi*lambda^2*vart))*exp(-(varA+v*vart)^2/(4*D*vart));                        

%% Calculation of the integral           
parfor a = 1:numel(A)  % do the numerical integration for every value of A; Note: This only holds for a single time point. Change the for loop to calculate more points in time
    F_C = 0;    % Set the integrals to 0
    F_D = 0;
    if A(a) > 0 % equation (14)
        for tt = 1:Nbrt
            F_C = F_C + IntC(A(a), tprime(tt)); % integration is done by a simple sum over the values of tprime and multiplication with dt in 3 lines
        end
        F1(a) = -1/2*exp(-(A(a)/lambda));   % first part of F
        F2(a) = -1/2*F_C*dt;                % second part of F
    elseif A(a) < 0 % equation (15)
        border = -A(a)/v;   % calculate the upper/lower border for the two integrals 
        borderidxA = find(tprime > border, 1, 'first'); % find the indices of the border
        borderidxB = find(tprime < border, 1, 'last');
        for tt = borderidxA:Nbrt    % calculate first integral
                F_C = F_C + IntC(A(a), tprime(tt));
        end
        for tt = 1:borderidxB       % calculate second integral
                F_D = F_D + IntD(A(a), tprime(tt)); 
        end
        F1(a) = 1/2*exp(A(a)/lambda);   % first part of F
        F2(a) = -1/2*(F_C - F_D)*dt;    % second part of F
    else % case A = 0
        F1(a) = 0;
        F2(a) = 0;
    end
    F(a) = F1(a) + F2(a); % sum up full F
end
x=-ll:ll;
pulse=F/trapz(x,abs(F));

end
