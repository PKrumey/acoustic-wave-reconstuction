% this a code for convolution of RCs using MRC multiplied by a asymmetric 
% lorentzian function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               input:
% ConvolutionFile:  matrix containing the angle in the first column and the
%                   calculated rocking curves in the other columns; rows
%                   correspond to angles, columns to time steps
% ResultFolder:     path for saving the obtained results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               output:
% convCRC:          matrix containing the convoluted calculated rocking
%                   curves
% MRCintp:          measured rocking curve interpolated on cTheta
% ConvFun:          function used for convolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% original author:  Mohammadmedhi Afshari
% v4, 02.07.2019, Fabian Brinks
% changelog:        - removed saving functions: This has to be done in 
%                     main program now  
%                   - removed convCRCthetaM as output
%                   - removed Pearson function for clarity
%                   - simplified variable names
% v3, 30.04.2019, Fabian Brinks
% changelog:        - introduced aa_left and aa_right as amplitudes of the two 
%                     half-lorentzians       
% v2, 07.11.2018, Fabian Brinks
% changelog:        - added comments and tidied up code


function [convCRC, ConvFun] = substrate_convolution_v5_GaAs111sm(ConvolutionFile,varargin)

    %% extract CRC out of convolution file and load measured rocking curve
    % Matrices containing angles and rocking curves
      MRCa = importdata('RC_Ti_GaAs111sm2.txt');
%    load('RC_Ti_GaAs111sm_new.mat');
    CRCa = ConvolutionFile; 
    
%     % lorentzian
%     SigmaLeft = 0.015;  SigmaRight = 0.04;
%     OffsetLeft = 0.875; OffsetRight = 0.95;	
    
%     % lorentzian
%     SigmaLeft = 0.02;  SigmaRight = 0.05;
%     OffsetLeft = 0.75; OffsetRight = 0.8125;

    % lorentzian
    SigmaLeft = 0.06;  SigmaRight = 0.06;
    OffsetLeft = 0.4; OffsetRight = 0.4;

    %% convolution
    [convCRC, MRCintp, ConvFun] = asymmetric_lorenz_convolution_v5_GaAs111(CRCa,MRCa,SigmaLeft,SigmaRight, OffsetLeft, OffsetRight);

    % plot the results if varargin = 'plot'
    if nargin > 1 && strcmp(varargin{1}, 'plot')
        figure; subplot(211)
        plot(CRCa(:,1),convCRC(:,1)/max(convCRC(:,1))); hold on;    % unpumped rocking curve
        plot(CRCa(:,1), MRCintp/max(MRCintp), 'r-');                % measured rocking curve
        plot(CRCa(:,1),ConvFun/max(ConvFun) ,'g-');                 % convolution function
        ax1 = gca; ax1.YLim = [1e-3 1]; %ax1.XLim = [-0.4 0.4]; 
        legend('convoluted CRC','Measured unpumped RC', 'convolution function');
        subplot(212);
        plot(CRCa(:,1),convCRC(:,1)/max(convCRC(:,1))); hold on;    % unpumped rocking curve
        plot(CRCa(:,1), MRCintp/max(MRCintp), 'r-');                % measured rocking curve
        plot(CRCa(:,1),ConvFun/max(ConvFun) ,'g-');                 % convolution function
        ax2 = gca; ax2.YScale = 'log'; ax2.YLim = [1e-4 1];% ax2.XLim = [-0.4 0.4]; 
        hold off;
    end
end