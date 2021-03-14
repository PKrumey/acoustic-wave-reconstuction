function fitness = calcFitness(CRC, MRC,theta,filter)
    %%%% calcFitness: Calculates the fitness value of the object            
    %                 normalize MRC and CRC to increase effect of 
    %                 side peaks

    [~, idx1] = min(abs(theta - filter(1)));
    [~, idx2] = min(abs(theta - filter(2)));
    [~, idx3] = min(abs(theta - filter(3)));
    [~, idx4] = min(abs(theta - filter(4))); 

    filtMRC1 = MRC(idx1:idx2,:);
    filtMRC2 = MRC(idx3:idx4,:);
    for i = 1:size(MRC,2)
        filtMRC1(:,i) = filtMRC1(:,i);%smooth(filtMRC1(:,i));
        filtMRC2(:,i) = filtMRC2(:,i);%smooth(filtMRC2(:,i));
    end
    filtMRC = [filtMRC1; filtMRC2];

    filtCRC = CRC([idx1:idx2 idx3:idx4],:);

    %filtMRC0 = [MRC0(idx1:idx2); MRC0(idx3:idx4)];%[smooth(MRC0(idx1:idx2)); smooth(MRC0(idx3:idx4))];

    %recalculating deveation into fitness
    fitness = logNorm(filtMRC, filtCRC);
    
    if(fitness==Inf)
       error("Infinit fitness value");
    end
end

