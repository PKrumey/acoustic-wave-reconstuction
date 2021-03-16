    function Sample = XRD_inputsV8_m_s(varargin)
% ************************************************************************
% this function provides a structure contains parameter used in XRD
% simulation
% ########################################################################
% INPUTS
% ########################################################################
% varargin{1} = F0_Re Structural factors in forward direction (source XOP) real part
% varargin{2} = F0_Im Structural factors in forward direction (source XOP) im part
% varargin{3} = Fpsi_Re Structural factors in an desired direction (source XOP) real part
% varargin{4} = Fpsi_Im Structural factors in an desired direction (source XOP) im part
% varargin{5} = DWF % Debye_waller factor in an desired direction
% varargin{6} = flag chooses the polarization [0 means s pol and 1 means p pol
% varargin{7} = Bragg angle (°)
% varargin{8} = lambda the wave length of the x-ray (m)
% varargin{9} = a unit cell constant (m)
% varargin{10} = R_electron  electron radius (m)
%#########################################################################
% OUTPUT
% ########################################################################
% Sample.KPol : polarization factor
% Sample.Vc : unit cell volum
% Sample.Thickness_to_Acomplex : the coefficient used to define reduced spatial coordinate
% Sample.g : a parameter in larson paper
% Sample.k : a parameter in larson paper
% Sample.Angle_to_y : the coefficient used to define reduced angular coordinate
% ########################################################################
% Mohammadmahdi Afshari
% 12-04-2017
% University of Duisburg-Essen
% ************************************************************************
switch nargin
    case 10
        if varargin{6} == 0
            Sample.KPol = 1; % s pol
        else
            Sample.KPol = cos(2*varargin{7}*pi/180);   % P polarization
        end
        
        Sample.Vc  = varargin{9}*varargin{9}*varargin{9};
        Sample.A_T = Sample.Vc./(Sample.KPol * varargin{10} * varargin{5} * varargin{3}*varargin{8});
        Sample.g_c = 0.5*(varargin{2}/(Sample.KPol*varargin{5}*varargin{3}));
        Sample.k = varargin{4}/varargin{3};
        Sample.y_a = (pi*Sample.A_T)/(2*varargin{8}); % this is for symetric case not for asymetric diffraction
        Sample.y_c = 0.5*(varargin{1}/(Sample.KPol*varargin{5}*varargin{3}));
        Sample.BraggAngle = varargin{7};
        Sample.pol = varargin{6};
    otherwise
        error(' Invalid number of input arguments. It has to contain 10 inputs')
end