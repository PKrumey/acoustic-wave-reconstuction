classdef population < handle
    % ************************************************************************
    % POPULATION A class where the optimization takes place
    % ########################################################################
    % Functions
    % ########################################################################
    % population:       Constructor to initialize first generation
    % sortGeneration:   sort Generation currGen according to the fitness
    % calculate:        calculates new Generation out of best Object in 
    %                   currGen
    % createDNA:        create new DNAObject that holds new bipolar pulse
    %                   data
    % evaluate:         evaluates currGen according to the fitness
    % plotResult:       plots result of each iteration
    % saveBestDNA:      saves best result when the optimization finished
    % ########################################################################
    % Philipp Krumey
    % 14-03-2021
    % University of Duisburg-Essen
    % ************************************************************************
    
    properties
        genSize             % number of DNA objects created
        
        currGen             % cell holding the DNA objects of current generation
        avFit = []          % average of fitness values
        bestFit = []        % maximum fitness value
        
        finished = false    % boolean for finish check

        theta               % theta vector
        time                % time vector
        filter              % angular area in which rocking curves are compared
        MRC                 % measured rocking curves   
        MRC0                % unpumped measured rocking curve
        CRC0                % unpumped calculated rocking curve
        
        pulseOrigin         % synthetic pulse
    end
    
    methods
        function obj = population(genSize, filter, MRC, theta,time,varargin)
            % ************************************************************************
            %POPULATION Construct an instance of the class 'population'
            %   assign class properties and calculate first generation
            % ########################################################################
            % INPUTS
            % ########################################################################
            % genSize:          number of elements per generation
            % filter:           pixels of rocking curve which are compared
            % MRC:              Matrix of measured rocking curves in the form 
            %                   [angle,time]. The first column has to be an unpumped 
            %                   rocking curve
            % theta:            Vector containing the angles for MRC
            % time:             Vector containing the time points of MRC
            % varargin:         containing precalculated first generation
            % ########################################################################
            % OUTPUT
            % ########################################################################
            % obj               population object
            % ************************************************************************
            
            %Initialize and assign class variables
            obj.genSize = genSize;           
            obj.theta = theta;
            obj.time = time; 
            obj.filter = filter;
            obj.MRC = MRC(:, 2:end);
            obj.MRC0 = MRC(:,1);
            obj.CRC0 = calcCRC(zeros(constants.ll,1),theta,0);
            
            %Check for pre calculated first generation, otherwise calculate first 
            %generation
            if nargin > 1 && iscell(varargin{1})
                obj.currGen = varargin{1};
            else
                obj.currGen = firstGeneration(obj.genSize,obj.theta,obj.time,obj.CRC0);
            end
            
            %calculate fitness for first generation and sort it according to the 
            %fitness
            for i=1:length(obj.currGen)
                obj.currGen{i}.fitness = calcFitness(obj.currGen{i}.CRC,obj.MRC,obj.theta,obj.filter);
            end  
            obj.sortGeneration();
        end
        
        function sortGeneration(obj)
            % ************************************************************************
            %sortGeneration Sort DNAObjects in currGen according to their fitness
            %   assign class properties and calculate first generation
            % ************************************************************************
            
            %Define Vector that holds the fitness parameter for the whole generation
            fitness = cellfun(@(x) x.fitness,obj.currGen);       
            
            %Sort currGen according to the fitness
            [~,ix] = sort(fitness);
            currGenSort = obj.currGen(ix);
            obj.currGen = currGenSort;
            
            %Save best and average fitness
            obj.bestFit = [obj.bestFit obj.currGen{1}.fitness];
            obj.avFit = [obj.avFit sum(fitness)/obj.genSize];               
        end
        
        function calculate(obj)
            % ************************************************************************
            %calculate main optimization routine containing iteration loop 
            % ************************************************************************
            
            %Select Fourier Coefficients of the best DNAObject
            sinCoefN = length(obj.currGen{1}.sinCoef);
            cosCoefN = length(obj.currGen{1}.cosCoef);
            coefN = sinCoefN+cosCoefN;                      %Number of fourier coefficients
            generation = 0;                                 %initialize iteration variable
            
            %main loop
            while(true)
                %array containing random fourier coefficient number
                coefi = randperm(coefN-1)+1;
                %loop over all fourier coefficients
                for i=1:(coefN-1)
                    bestDNA = obj.currGen{1};               %select best DNAObject from last generation
                    nextGen = cell(obj.genSize,1);          %initialize cell array for next generation
                    nextGen{1} = bestDNA;                   %insert best DNAObject into the first space of the next generation
                    
                    %create arrays for all fourier coefficients of the next generation
                    sinCoef = ones(obj.genSize-1,sinCoefN).*bestDNA.sinCoef;
                    cosCoef = ones(obj.genSize-1,cosCoefN).*bestDNA.cosCoef;
                    
                    coefArea = 5e-3;                        %fourier coefficient area for the next generation
                    
                    %create new fourier coefficents in the area coefArea around the best fourier coefficient
                    if coefi(i)<=sinCoefN
                        coefArea = coefArea*(1-0.1*(coefi(i)-2)/(sinCoefN-2));
                        sinCoef(:,coefi(i)) = (bestDNA.sinCoef(coefi(i))-coefArea):2*coefArea/(obj.genSize-2):(bestDNA.sinCoef(coefi(i))+coefArea);
                    else
                        cosCoef(:,coefi(i)-sinCoefN) = (bestDNA.cosCoef(coefi(i)-sinCoefN)-coefArea):2*coefArea/(obj.genSize-2):(bestDNA.cosCoef(coefi(i)-sinCoefN)+coefArea);
                    end
                   
                    %create new generation out of calculated fourier coefficients
                    f0 = waitbar(0);
                    for j=1:(obj.genSize-1)
                        waitbar(j/obj.genSize,f0,['Calculating Object Number ' num2str(j) ' of ' num2str(obj.genSize)]);
                        nextGen{j+1} = obj.createDNA(sinCoef(j,:),cosCoef(j,:));
                    end
                    obj.currGen = nextGen;                  %save new generation
                    close(f0);
                    obj.evaluate(coefi(i),generation);      %evaluate new generation
                    
                    %check termination condition
                    if obj.finished==true
                       return
                    end
                end  
                generation = generation + 1;                %increase iteration variable
            end
        end
        
        function DNAObject = createDNA(obj,sinCoef,cosCoef)
            % ************************************************************************
            %DNAObject create new DNAObject from fourier coefficients 
            % ########################################################################
            % INPUTS
            % ########################################################################
            % sinCoef:          sinus coefficients for new bipolar pulse
            % cosCoef:          cosine coefficients for new bipolar pulse
            % ########################################################################
            % OUTPUT
            % ########################################################################
            % DNAObject:        new created DNAObject           
            % ************************************************************************
            DNAObject = DNA();                                                                %initialize DNAObject
            DNAObject.sinCoef = sinCoef;                                                      %save sine coefficients  
            DNAObject.cosCoef = cosCoef;                                                      %save coseine coefficients
            DNAObject.pulse = fSeries(DNAObject.sinCoef,DNAObject.cosCoef,constants.ll);      %create and save new bipolar pulse
            CRC = norm2unp([obj.CRC0 calcCRC(DNAObject.pulse,obj.theta,obj.time)]);           %calculate normed rocking curve for new bipolar pulse
            DNAObject.CRC     = CRC(:,2:end);                                                 %save normed calculated rocking curve
            DNAObject.fitness = calcFitness(DNAObject.CRC,obj.MRC,obj.theta,obj.filter);      %calculate fitness
            DNAObject.createID();                                                             %generate ID for new DNAObject
        end
              
        function evaluate(obj,coef,generation)
            % ************************************************************************
            %%%% evaluate: selects best bipolar pulse of the generation, prints out
            %           result and sets termination condition
            % ########################################################################
            % INPUTS
            % ########################################################################
            % coef:          %number which determines the optimized coefficient
            % generation:    %current generation number      
            % ************************************************************************
            
            %sort current generation
            obj.sortGeneration();
            
            %print coefficient and generation number to command window
            FormatSpec = 'Generation Number %i \n';
            fprintf(FormatSpec,generation);          
            
            FormatSpec = 'Coefficient Number %i \n';
            fprintf(FormatSpec,coef);
           
            % print best and average fitness to command window
            FormatSpec = 'Average Fitness = %0.5e \nBest Fitness = %0.5e \n\n';
            fprintf(FormatSpec, obj.avFit(end),obj.bestFit(end));

            % plot the average and best fitness and best pulse
            obj.plotResult;
        
            % check termination condition
            if length(obj.bestFit)>11
                if obj.bestFit(end-10)<=obj.bestFit(end)                   
                   obj.finished=true; 
                end                          
            end
        end
        
        function plotResult(obj)
                % ************************************************************************
                %%%%plotResult: plots the best average and worste pulse of the generation
                %           together with the best and average fitness and the fourier
                %           coefficients
                % ************************************************************************
                
                %initialize or select figure
                hfig1 = figure(1);
                set(hfig1,'Name','average, best and original pulse')
                
                %plot best and average fitness
                subplot(3,1,1);
                plot(obj.bestFit,'LineWidth',2,'DisplayName','Best Fitness');
                hold on;         
                plot(obj.avFit,'LineWidth',2,'DisplayName','Average Fitness');
                hold off
                set(gca,'fontsize',12)
                set(gca,'fontweight','bold')
                set(gca,'linewidth',2)
                xlabel('Generation')
                ylabel('Fitness')
                legend
                
                %plot best, average, worst bipolar pulse and if available the synthetic pulse 
                subplot(3,1,2);
                plot(-constants.ll:constants.ll,obj.currGen{1}.pulse,'LineWidth',2,'DisplayName','Best Pulse')
                hold on
                plot(-constants.ll:constants.ll,obj.currGen{floor(obj.genSize/2)}.pulse,'LineWidth',2,'DisplayName','Average Pulse');
                plot(-constants.ll:constants.ll,obj.currGen{end}.pulse,'LineWidth',2,'DisplayName','Worst Pulse');
                if ~isnan(obj.pulseOrigin)
                    plot(-constants.ll:constants.ll,obj.pulseOrigin,'LineWidth',2,'DisplayName','Original Pulse');
                end
                hold off
                set(gca,'fontsize',12)
                set(gca,'fontweight','bold')
                set(gca,'linewidth',2)
                xlabel('Length [nm]')
                ylabel('Strain')
                legend
                
                %plot best sine coefficients and if available synthetic sine coefficients
                subplot(3,1,3);
                plot(obj.currGen{1}.sinCoef,'.','MarkerSize',24,'DisplayName','Best Coefs')
                if ~isnan(obj.pulseOrigin)
                    hold on
                    plot(SinCoefs(obj.pulseOrigin,length(obj.currGen{1}.sinCoef),constants.ll),'.','MarkerSize',24,'DisplayName','Original Coefs');
                    hold off
                end
                set(gca,'fontsize',12)
                set(gca,'fontweight','bold')
                set(gca,'linewidth',2)
                xlabel('Sine Coefficient')
                ylabel('Strain')
                legend
        end
        
        function saveBestDNA(obj)
            % ************************************************************************
            %%%%saveBestDNA: save the best DNAObject of the reconstruction
            % ************************************************************************
        
            %select or make the directory where the results will be saved
            warning ('off','all');
            mkdir('Results');
            warning ('on','all');
            fileName = ['Results\' datestr(now, 'yyyy-mm-dd_HH-MM')];
            
            %write the optimization parameters into a struct
            inf.theta = obj.theta;
            inf.time = obj.time;
            inf.filter = obj.filter;
            inf.genSize =obj.genSize;
            
            %write the best result into a struct
            inf.bestDNA = obj.currGen{1}; 
            
            %save struct
            save(fileName,'inf');
            
            %plot the result for synthetic data
            if ~isnan(obj.pulseOrigin)
               plotGACRCResult(inf,obj.pulseOrigin); 
            end
        end
    end
end

