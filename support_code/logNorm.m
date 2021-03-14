function fitness = logNormV2(MRC, CRC)
%logNorm: Normalization is done by taking the logarithm of both matrices
        normMRC = log(MRC+1e-4);
        normCRC = log(CRC+1e-4);
        [thetas,times] = size(CRC);
        % fitness is sum of deviation squares
        fitness = sum(sum((normMRC - normCRC).^2)/thetas)/times;
        
end

