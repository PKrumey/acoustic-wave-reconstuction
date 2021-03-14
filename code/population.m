classdef population < handle
    %POPULATION A class that holds the DNA structure
    %   
    
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
        CRC0
        
        pulseOrigin         % Original Pulse
    end
    
    methods
        function obj = population(genSize, filter, MRC, theta,time,varargin)
            %POPULATION Construct an instance of the class 'population'
            %   assign class properties
            obj.genSize = genSize;
            
            obj.theta = theta;
            obj.time = time; 
            obj.filter = filter;
            obj.MRC = MRC(:, 2:end);
            obj.MRC0 = MRC(:,1);
            obj.CRC0 = calcCRC(zeros(constants.ll,1),theta,0);
            
            if nargin > 1 && iscell(varargin{1})
                obj.currGen = varargin{1};
            else
                obj.currGen = firstGeneration(obj.genSize,obj.theta,obj.time,obj.CRC0);
            end
            
            for i=1:length(obj.currGen)
                obj.currGen{i}.fitness = calcFitness(obj.currGen{i}.CRC,obj.MRC,obj.theta,obj.filter);
            end  
            obj.sortGeneration();
        end
        
        function sortGeneration(obj)
            fitness = cellfun(@(x) x.fitness,obj.currGen);
            
            figure(2);
            plot(fitness)            
            
            [~,ix] = sort(fitness);
            currGenSort = obj.currGen(ix);
            obj.currGen = currGenSort;
            
            obj.bestFit = [obj.bestFit obj.currGen{1}.fitness];
            obj.avFit = [obj.avFit sum(fitness)/obj.genSize];               
        end
        
        function calculate(obj)
            sinCoefN = length(obj.currGen{1}.sinCoef);
            cosCoefN = length(obj.currGen{1}.cosCoef);
            coefN = sinCoefN+cosCoefN;
            generation = 0;
            while(true)
                coefi = randperm(coefN-1)+1;
                for i=1:(coefN-1)
                    bestDNA = obj.currGen{1};
                    nextGen = cell(obj.genSize,1);
                    nextGen{1} = bestDNA;
                    
                    sinCoef = ones(obj.genSize-1,sinCoefN).*bestDNA.sinCoef;
                    cosCoef = ones(obj.genSize-1,cosCoefN).*bestDNA.cosCoef;
                    
                    coefArea = 5e-3;
                    if coefi(i)<=sinCoefN
                        coefArea = coefArea*(1-0.1*(coefi(i)-2)/(sinCoefN-2));
                        sinCoef(:,coefi(i)) = (bestDNA.sinCoef(coefi(i))-coefArea):2*coefArea/(obj.genSize-2):(bestDNA.sinCoef(coefi(i))+coefArea);
                    else
                        cosCoef(:,coefi(i)-sinCoefN) = (bestDNA.cosCoef(coefi(i)-sinCoefN)-coefArea):2*coefArea/(obj.genSize-2):(bestDNA.cosCoef(coefi(i)-sinCoefN)+coefArea);
                    end
                    
                    f0 = waitbar(0);
                    for j=1:(obj.genSize-1)
                        waitbar(j/obj.genSize,f0,['Calculating Object Number ' num2str(j) ' of ' num2str(obj.genSize)]);
                        nextGen{j+1} = obj.createDNA(sinCoef(j,:),cosCoef(j,:));
                    end
                    obj.currGen = nextGen;
                    close(f0);
                    obj.evaluate(coefi(i),generation);
                    if obj.finished==true
                       return
                    end
                end  
                generation = generation + 1;
            end
        end
        
        function DNAObject = createDNA(obj,sinCoef,cosCoef)
            DNAObject = DNA();
            DNAObject.sinCoef = sinCoef;
            DNAObject.cosCoef = cosCoef;
            DNAObject.pulse = fSeries(DNAObject.sinCoef,DNAObject.cosCoef,constants.ll);
            CRC = norm2unp([obj.CRC0 calcCRC(DNAObject.pulse,obj.theta,obj.time)]);
            DNAObject.CRC     = CRC(:,2:end);
            DNAObject.fitness = calcFitness(DNAObject.CRC,obj.MRC,obj.theta,obj.filter);
            DNAObject.createID();
        end
              
        function evaluate(obj,coef,generation)
            %%%% evaluate: Checks wether the objectiveScore is reached and
            %              stops the algorithm
            obj.sortGeneration();
            
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
                hfig1 = figure(1);
                set(hfig1,'Name','average, best and original pulse')
                
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
            warning ('off','all');
            mkdir('Results');
            warning ('on','all');
            fileName = ['Results\' datestr(now, 'yyyy-mm-dd_HH-MM')];
            
            inf.theta = obj.theta;
            inf.time = obj.time;
            inf.filter = obj.filter;
            inf.genSize =obj.genSize;
            
            inf.bestDNA = obj.currGen{1}; 
            
            save(fileName,'inf');
            
            if ~isnan(obj.pulseOrigin)
               plotGACRCResult(inf,obj.pulseOrigin); 
            end
        end
    end
end

