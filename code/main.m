function main(workspace,genSize,filter,firstGen)
    % ************************************************************************
    %MAIN routine to start optimization algorithm
    %   This is the main Function for running the optimization algorithm to retrieve the
    %   strain pulse out of the measured rocking curves for different time
    %   delays.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                               Input:
    % workspace:
    %    MRC:      Matrix of measured rocking curves in the form (angle)x(time). The
    %              first column has to be an unpumped rocking curve, which is used for
    %              the calculation of the perfect achievable score
    %    theta:    Vector containing the angles of the MRC
    %    time:     Vector containing the time points of the MRC. Note: t = 0 is
    %              not part of the vector time!
    % genSize:  number of elements per population
    % filter:   pixels of the rocking curve which are compared e.g. [-0.4 -0.1 0.1 0.4]
    % firstGen: (optional) gives first generation as an array
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Philipp Krumey
    % 14-03-2021 
    % University of Duisburg-Essen
    % *************************************************************************

    %% add workspace
    workspacecell = load(workspace);
    
    % print Simulation start
    disp(['Simulation started: ' datestr(now)]);
    tic;    
    
    %% create first generation
    pop = population(genSize, filter, workspacecell.MRC, workspacecell.theta, workspacecell.time, firstGen);
    
    %% check for bipolare pulse in workspace
    if isfield(workspacecell,'pulse')
       pop.pulseOrigin = workspacecell.pulse; 
    end
    
    %% start evolution
    pop.calculate;

    %save best result
    pop.saveBestDNA();
    
    % print Simulation end
    t=toc;
    disp(['Total computation time: ' datestr(datenum(0,0,0,0,0,t),'HH:MM:SS')]);
    disp(['Simulation finished: ' datestr(now)]);
    disp('DONE')   
end

