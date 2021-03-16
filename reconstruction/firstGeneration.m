function firstGen = firstGeneration(genSize,theta,time,CRC0)
    % ************************************************************************
    % FIRSTGENERATION generate and save first generation
    % ########################################################################
    % INPUTS
    % ########################################################################
    % genSize:  number of elements per generation
    % theta:    Vector containing the angles for CRC
    % time:     Vector containing the time points for CRC
    % CRC0:     calculated unpumped rocking curve
    % ########################################################################
    % OUTPUT
    % ########################################################################
    % firstGen: cell-array conainting the first generation [genSize]   
    % ########################################################################
    % Philipp Krumey
    % 14-03-2021
    % University of Duisburg-Essen
    % ************************************************************************
    
    %calculate first generation and save into firstGen
    firstGen = cell(genSize, 1);
    f0 = waitbar(0);
    for i = 1:genSize
        waitbar(i/genSize,f0,['Calculating Object Number ' num2str(i) ' of ' num2str(genSize)]);
        firstGen{i} = createDNAObject(theta,time,CRC0);  % creates DNA objects
    end
    close(f0);
end

function DNAObject = createDNAObject(theta,time,CRC0)
    % ************************************************************************
    %%%% firstGeneration: Generates a random pulse shape and
    %                     calculates its CRC
    % ########################################################################
    % INPUTS
    % ########################################################################
    % theta:    Vector containing the angles for CRC
    % time:     Vector containing the time points for CRC
    % CRC0:     calculated unpumped rocking curve
    % ########################################################################
    % OUTPUT
    % ########################################################################
    % DNAObject: Object containing bipolar pulse data  
    % ************************************************************************    
    
    DNAObject = DNA();                  %inialize DNAObject
    pulse  = Thomsen(constants.ll)';    %generate random bipolar pulse
    x = -constants.ll:constants.ll;     
    
    %calculate Fourier coefficients and fourier expansion from pulse
    ks = 9;                             %number of sin coefficients
    sinCoef = SinCoefs(pulse,ks,constants.ll);
    kc=3;                               %number of cos coefficients
    cosCoef = CosCoefs(pulse,kc,constants.ll);
    pulse = fSeries(sinCoef,cosCoef,constants.ll);

    %calculate random strain level for pulse and fourier coefficients
    Area = trapz(x,abs(pulse));
    randA = (0.4*rand+0.8);              %random mean strain level between (0.8-1.2)/(2*constants.ll)
    pulse = randA*pulse/Area;
    sinCoef = randA*sinCoef/Area;
    cosCoef = randA*cosCoef/Area;
        
    %save bipulse data into DNA Object
    DNAObject.pulse = pulse;
    DNAObject.sinCoef = sinCoef;
    DNAObject.cosCoef = cosCoef;
    CRC = norm2unp([CRC0 calcCRC(pulse,theta,time)]);       %calculate normed rocking curve
    DNAObject.CRC     = CRC(:,2:end);
    DNAObject.createID();
end

