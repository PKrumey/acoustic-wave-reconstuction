classdef constants
    % ************************************************************************
    %CONSTANTS Structure that holds the TRXD constants for GaAs 111
    % ########################################################################
    % Philipp Krumey
    % 14-03-2021
    % University of Duisburg-Essen
    % ************************************************************************  
    
    properties (Constant = true)
        %% substrate properties
        dSubs       =   15e-6;                  % substrate thickness (m)
        svSubs      =   5400;                   % Speed of sound
        F0_Re       =   254.958412168;          % Structural factors in forward direction (source XOP) real part
        F0_Im       =   20.2621596;             % Structural factors in forward direction (source XOP) imaginary part
        Fpsi_Re     =   155.09282806143935;     % Structural factor for (111) reflection (source XOP)  real part
        Fpsi_Im     =   14.379811817126567;     % Structural factor for (111) reflection (source XOP)  imaginary part
        DWF         =   1;                      % Debye_waller factor
        Theta_B     =   24.898521445843546;     % Bragg angle for the substrate (�)
        lambda      =   2.748486e-10;           % X-ray wave length (m)
        a           =   5.65325e-10;            % Unit cell constant (m)
        R_elctron   =   2.817940325*1e-15;      % classical electron radius (m)
        Pol_flag    =   0.5;                    % This flag determins the polarization of the X-Ray [0 means s pol, 1 means p pol, 0.5 unpol]
        
        %% fixed properties of strain pulse
        R           =   0.3;                    % acoustic reflectivity of the metal/GaAs interface
        ll          =   88;                     % length of single comp and exp in GaAs (nm)
        maxStrain   =   4.5e-3;                 % maximum strain amplitude
    end  
end

