function fitness = calcFitness(CRC, MRC,theta,filter)
% ************************************************************************
% calcFitness:    Calculates the fitness value of the object            
%                 normalize MRC and CRC to increase effect of 
%                 side peaks
% ########################################################################
% INPUTS
% ########################################################################
% CRC:      matrix of the calculated rocking curve in the form 
%           [theta,time]
% MRC:      measured rocking curve in the form [theta,time]
% theta:    Vector containing the angles of the MRC
% filter:   Vector containing the angular area, that should be taken into 
%           acount to calculate the fintess 
%           [anglemin1,anglemax1,anglemin2,anglemax2]
%#########################################################################
% OUTPUT
% ########################################################################
% fitness:  positive number holding the deveation between CRC and MRC
% ########################################################################
% Philipp Krumey
% 15-03-2021, University of Duisburg-Essen
% *************************************************************************

    %calculate the indices, corresponding to anglemin1,anglemax1,anglemin2,
    %anglemax2
    [~, idx1] = min(abs(theta - filter(1)));
    [~, idx2] = min(abs(theta - filter(2)));
    [~, idx3] = min(abs(theta - filter(3)));
    [~, idx4] = min(abs(theta - filter(4))); 

    %filter MRC and save filtered data in filtMRC
    filtMRC = MRC([idx1:idx2 idx3:idx4],:);

    %filter CRC and save filtered data in filtCRC
    filtCRC = CRC([idx1:idx2 idx3:idx4],:);

    %recalculating deveation into fitness
    fitness = logNorm(filtMRC, filtCRC);
    
    if(fitness==Inf)
       error("Infinit fitness value");
    end
end

