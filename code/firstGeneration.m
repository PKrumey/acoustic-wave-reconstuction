function firstGen = firstGeneration(genSize,theta,time,CRC0)
%FIRSTGENERATION generate and save first generation
    firstGen = cell(genSize, 1);
    f0 = waitbar(0);
    for i = 1:genSize
        waitbar(i/genSize,f0,['Calculating Object Number ' num2str(i) ' of ' num2str(genSize)]);
        firstGen{i} = createDNAObject(theta,time,CRC0);  % creates DNA objects
    end
    close(f0);
end

function DNAObject = createDNAObject(theta,time,CRC0)
    %%%% firstGeneration: Generates a random pulse shape and
    %                     calculates its CRC
    DNAObject = DNA();
    pulse  = Thomsen(constants.ll)';
    x = -constants.ll:constants.ll;
    
    ks = 9;
    sinCoef = SinCoefs(pulse,ks,constants.ll);
    kc=3;
    cosCoef = CosCoefs(pulse,kc,constants.ll);
    pulse = fSeries(sinCoef,cosCoef,constants.ll);

    Area = trapz(x,abs(pulse));
    randA = (1.6*rand+0.2);
    pulse = randA*pulse/Area;
    sinCoef = randA*sinCoef/Area;
    cosCoef = randA*cosCoef/Area;
        
    DNAObject.pulse = pulse;
    DNAObject.sinCoef = sinCoef;
    DNAObject.cosCoef = cosCoef;
    CRC = norm2unp([CRC0 calcCRC(pulse,theta,time)]);
    DNAObject.CRC     = CRC(:,2:end);
    DNAObject.createID();
end

