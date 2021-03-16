classdef DNA < handle
    %DNA: Class of CRCs
    %   This is a Data Class consisting of the calculated bipolar Pulses
    %   and their Fourier Coefficients. It also holds the calculated
    %   Rocking Curves and the Fitness as a quantity to determine the
    %   difference to the input Data
    % ########################################################################
    % Philipp Krumey
    % 14-03-2021, University of Duisburg-Essen
    % *************************************************************************
    
    properties     
        ID      % gives a unique ID to every DNA Object
        pulse   % shape of the pulse of a DNA object
        sinCoef % Uneven Fourier coefficients
        cosCoef % Even Fourier ceofficients
        CRC     % CRC matrix of a DNA object
        fitness % Fitness Value of the DNA object
    end
    
    methods
        function obj = DNA()
            %%%% DNA: Construction method assign class properties
        end
       
        function createID(obj)
            %%%% createID: creates unique ID            
            obj.ID = DataHash(obj.pulse);
        end
        
        function s = saveobj(obj)
            %%%% saveobj: outputs the DNA parameters as a struct for saving 
            s.ID = obj.ID;
            s.pulse = obj.pulse;
            s.sinCoef = obj.sinCoef;
            s.cosCoef = obj.cosCoef;
            s.CRC = obj.CRC;
            s.fitness = obj.fitness;
        end
    end
    
    methods(Static)
        function obj = loadobj(a)
            %%%% loadobj: creates DNA Object from input parameter a
            if isstruct(a)                
                newObj = DNA(); 
                newObj.pulse = a.pulse;
                newObj.ID = a.ID;                
                newObj.sinCoef = a.sinCoef;
                newObj.cosCoef = a.cosCoef;
                newObj.CRC = a.CRC;
                newObj.fitness = a.fitness;
                obj = newObj;
            else
                obj = a;
            end
        end
    end
end
