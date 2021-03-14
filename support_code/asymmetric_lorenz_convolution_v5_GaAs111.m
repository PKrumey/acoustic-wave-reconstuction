% this function builds a convolution function containing the measured
% rocking curve, an asymmetric lorentzian function and some constant values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   input:
% CRCa:             matrix with calculated rocking curves; rows
%                   correspond to angles, columns to time steps
% MRCa:             measured rocking curve together with angles
% Sigma*:           left and right lorentzian parameter
% Offset*:          offsets of the left and right lorentzian curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   output:
% convCRC:          matrix containing the convoluted calculated rocking
%                   curves
% ConvFun:          normalized convolution function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% original author: Mohammadmedhi Afshari
% v4, 02.07.2019, Fabian Brinks
% changelog:    - refraction correction is done in the RC instead of the
%                 theta vector; theta is not conveyed via convCRC anymore
%               - removed MRCintp and convCRCmTheta as outputs
%               - simplified variable names and removed some useless staff
% v3, 30.04.2019, Fabian Brinks
% changelog:    - introduced aa_left and aa_right as parameters for the
%                 amplitude of the half-lorentzians
% v2, 07.11.2018, Fabian Brinks
% changelog:    - fixed a bug, that MRC connection to theta_conv was wrong
%               - added comments and tidied up code
               

function [convCRC, MRCintp, ConvFun] = asymmetric_lorenz_convolution_v5_GaAs111(CRCa,MRCa,SigmaLeft,SigmaRight, OffsetLeft, OffsetRight)
    %% prepare data
    % refraction correction for theta
    cTheta      = CRCa(:,1);            % extract theta
    refrCorr    = 0.007174;             % correction value (XOP)
    ThetaCorr   = cTheta - refrCorr;    % refraction correction
    CRC         = CRCa(:,2:end);        % extract rocking curves

    % extract measured rocking curve, angle and normalize it
    mTheta      = MRCa(:,1);
    MRC         = MRCa(:,2);

    %  define theta for convolution
    dTheta      = 0.00001;
    convTheta   = -max(abs(ThetaCorr)):dTheta:max(abs(ThetaCorr));
    convTheta   = convTheta';
    NumAngles   = length(convTheta);
    [~, idx0]   = min(abs(convTheta));      %find position of theta = 0

    % interpolate the measured rc to all convTheta, normalize it to the 
    % integral being 1 and shift it such that the maximum is at theta = 0
    MRCintp     = interp1(mTheta, MRC, convTheta, 'linear', 0);
    MRCintp     = MRCintp/(sum(MRCintp)*median(diff(convTheta)));     % normalize rc
    [~, idxmax] = max(MRCintp);                                     % find the max index
    MRCintp     = circshift(MRCintp, (idx0 - idxmax));              % shift MRC 

    % interpolate calculated rocking curves to convTheta
    NbrTimeSteps    = size(CRC,2);
    CRCintp         = zeros(NumAngles, NbrTimeSteps);
    
    for i = 1:NbrTimeSteps
        % use of linear interpolation method is crucial here
        CRCintp(:,i)      = interp1(ThetaCorr, CRC(:,i) ,convTheta, 'linear', 0);
    end

    %% define convolution function
    % calculate the lorentzian
    LorLeft    =   @(x) OffsetLeft + 0.95./(1+4*(x/SigmaLeft).^2);        
    LorRight   =   @(x) OffsetRight + 0.85./(1+4*(x/SigmaRight).^2);

    ThetaConvLeft   = convTheta(1:idx0);          % define angel variables
    ThetaConvRight  = convTheta(idx0+1:end);

    L = [LorLeft(ThetaConvLeft); LorRight(ThetaConvRight)];

    % multiply with MRC to get Convolution function
    ConvFun     =   L.*MRCintp;
    ConvFun     =   ConvFun/(sum(ConvFun)*dTheta);

    %% take convolution
    convCRC      =     zeros(length(ThetaCorr),NbrTimeSteps);

    for i = 1:NbrTimeSteps
         temp_g             = dTheta*conv(CRCintp(:,i), ConvFun, 'same');           % convolution
         temp_g             = circshift(temp_g, - int32(refrCorr/dTheta));          % shift matrix for refraction correction
         temp_g             = interp1(convTheta, temp_g, ThetaCorr, 'linear');         % This interpolates the convCRC on the original simulation theta
         convCRC(:,i)       = temp_g;
    end
    
    % interpolate convolution function to cTheta for handing it to main
    % function
    
    ConvFun     = interp1(convTheta, ConvFun, cTheta, 'linear');
    MRCintp     = interp1(convTheta, MRCintp, cTheta, 'linear');
end