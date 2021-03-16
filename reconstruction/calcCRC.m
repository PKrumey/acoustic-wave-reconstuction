function convCRC = calcCRC(pulse,theta,time)
% ************************************************************************
% calcCRC:  This function calculates the convoluted rocking curve           
% ########################################################################
% INPUTS
% ########################################################################
% pulse:    Vector containing the bipolare strain pulse
% theta:    Vector containing the angles for convCRC
% time:     Vector containing the time points for convCRC
%#########################################################################
% OUTPUT
% ########################################################################
% convCRC:  Array containing the calculated convoluted rocking curves
%           [theta,time]
% ########################################################################
% Philipp Krumey
% 15-03-2021, University of Duisburg-Essen
% *************************************************************************
    % calculate rocking curve depending on polarization flag in constants
    if constants.Pol_flag == 1 || constants.Pol_flag == 0
        % calulate rocking curves for polarized x-rays
        CRC = calcCRCPol(pulse,theta,time,constants.Pol_flag);
    else
        % calulate rocking curves for un or partly polarized x-rays
        CRCspol = calcCRCPol(pulse,theta,time,0);
        CRCppol = calcCRCPol(pulse,theta,time,1);
        CRC = constants.Pol_flag*CRCspol+(1-constants.Pol_flag)*CRCppol;
    end

    % convolute the calculated rocking curve
    ConvolutionFile =   [theta CRC];
    [convCRC, ~] = substrate_convolution_v5_GaAs111sm(ConvolutionFile);
    convCRC(convCRC == 0) = 1e-9;        
end

function CRC = calcCRCPol(pulse,theta,time,Pol_flag)
    % ************************************************************************
    %% calcCRC: This function calculates the rocking curve with s or p
    %               Polarisation
    % ########################################################################
    % INPUTS
    % ########################################################################
    % pulse:    Vector containing the bipolare strain pulse
    % theta:    Vector containing the angles for convCRC
    % time:     Vector containing the time points for convCRC
    % Pol_flag  This flag determins the polarization of the X-Ray 
    %           [0 means s pol, 1 means p pol
    %#########################################################################
    % OUTPUT
    % ########################################################################
    % CRC:  Array containing the calculated rocking curves
    %       [theta,time]    
    % *************************************************************************
    
    % calculate parameters for DXRD
    subs = XRD_inputsV8_m_s(constants.F0_Re, constants.F0_Im, ...
        constants.Fpsi_Re, constants.Fpsi_Im, constants.DWF, Pol_flag, ...
        constants.Theta_B, constants.lambda, constants.a, constants.R_elctron);

    % build strain array
    strain = buildStrain(pulse,time);

    % prepare parfor variables; note: cells are used to provide 'sliced
    % variables' for parallelization
    NbrTimeSteps    = length(time);             %number of times
    N               = length(theta);            % number of angles
    RCcell          =   cell(N,1);	            % cell for RCs
    Scell           =   cell(1, NbrTimeSteps);	% put strain in a cell structure for sliced variable use in parfor

    for i = 1:NbrTimeSteps                      % find the appropriate column of S and put it into the cells
        Scell{i} = strain(i, :);
    end

    parfor j = 1:N                              % loop over angles
        RCcell{j}  =   zeros(1,NbrTimeSteps);   % cell that holds rocking curves

        for tt = 1:NbrTimeSteps                 % loop over times
            S_z        =   Scell{tt}';          % pick correct strain 
            ang        =   theta(j);            % pick correct angle             

            % solve the DXRD equation                                                                 
            [A, X] = DXRD_ode23_C_mex(ang, subs.k, subs.g_c, subs.y_a, subs.y_c, subs.BraggAngle, subs.A_T, constants.dSubs, S_z)
            RCcell{j}(tt)	=   (X(length(A),1).^2 + X(length(A),2).^2)';	%calculation of the X-ray reflectivity
        end
    end

    CRC  = zeros(N, NbrTimeSteps);              % put the RCs into one matrix
    for i = 1:N  
        CRC(i,:)     = RCcell{i};
    end
end

function strain = buildStrain(pulse,time)
    % ************************************************************************
    %% buildStrain: This function takes a bipolar pulse and creates a strain
    %                 matrix for the given time vector
    % ########################################################################
    % INPUTS
    % ########################################################################
    % pulse:    Vector containing the bipolare strain pulse
    % time:     Vector containing the time points for strain
    %#########################################################################
    % OUTPUT
    % ########################################################################
    % strain:   Array containing the strain pulsetrains for all times 
    %           [time,pulsetrain]    
    % *************************************************************************

    % create pulse train from pulses and with reflectivity R
    fullPulseTrain = [constants.R^7*pulse constants.R^6*pulse constants.R^5*pulse constants.R^4*pulse constants.R^3*pulse constants.R^2*pulse constants.R*pulse pulse];
 
    % initialize strain array
    S       =   zeros(length(time), 2000);  

    % build strain array out of pulse train
    for i = 1:length(time)
        z_t = round(time(i)*1e-12*constants.svSubs/1e-9);
        if length(fullPulseTrain) > z_t
            S(i, 1:z_t) =   fullPulseTrain(end - z_t + 1: end);
        else
            S(i, z_t - length(fullPulseTrain) + 1 : z_t) = fullPulseTrain;
        end  
    end 
    strain = S;         
end


