classdef DNA < handle
    %DNA: Class of CRCs
    %   This is a Data Class consisting of the calculated bipolar Pulses
    %   and their Fourier Coefficients. Also it holds the calculated
    %   Rocking Curves and the Fitness as a quantity to determine the
    %   difference to the input Data
    
    properties     
        ID      % gives a unique ID to everypulse shape
        pulse   % shape of the pulse of a DNA object
        sinCoef % Uneven Fourier coefficients
        cosCoef % Even Fourier ceofficients
        CRC     % CRC matrix of a DNA object
        fitness % Fitness Value of the DNA object
    end
    
    methods
        function obj = DNA()
            %%%% DNA: Construction method            
            %         assign class properties
        end
       
        function createID(obj)
            obj.ID = DataHash(obj.pulse);
        end
        
        function s = saveobj(obj)
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
