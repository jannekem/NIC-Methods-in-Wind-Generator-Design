function PSO_function(arrayTaskNumber, particleAmount, t_end, phi)
addpath('model');
%% INITIALIZTION OF SWARM VARIABLES
% initialize random generator
rng(sum((arrayTaskNumber+100)*clock));

% parameter limits, last value indicates if the value is considered integer
limits = [20 80 1; 
    35000 65000 0; 
    2000000 6000000 0;
    0.001 0.05 0;
    0.8 5 0; 
    0.6 0.95 0; 
    21000 48000 0;
    1.3 1.6 0;
    1 3 1;
    0.8 0.99 0;
    0.25 0.75 0;
    0.75 0.95 0;
    0.01 0.1 0;
    0.01 0.1 0];
parameterAmount = size(limits,1);       % amount of optimizable parameters
integerIndices = find(limits(:,3)==1);  % find the indices for the integer values
gamma1 = 2;
gamma2 = 2;


% generate the swarm and velocity matrices
swarm = zeros(particleAmount, parameterAmount);
velocity = zeros(particleAmount, parameterAmount);
bestpars = zeros(particleAmount, parameterAmount);  % best parameters for each particle 
bestresults = zeros(particleAmount,1);
bestidx = 0;        % holder for the global best index
bestval = 0;        % global best value
bestvals = zeros(t_end,1); % best value vector

% initialize swarm
for ind = 1:particleAmount 
    for par = 1:parameterAmount
       swarm(ind,par) = limits(par,1)+rand()*(limits(par,2)-limits(par,1));
    end
end


%% EVALUATE INITIAL VALUES
for ind = 1:particleAmount
    params = squeeze(swarm(ind,:));
    params(integerIndices) = round(params(integerIndices)); % round to the nearest integer
    [objects, constraints] = design_PMSM_generator(params);
    bestpars(ind,:) = swarm(ind,:);
    val = evalobjects(objects, constraints);
    bestresults(ind) = val;
end

[bestval, bestidx] = max(bestresults);  % evaluate the initial best results
bestvals(1) = bestval;
itercount = 1;
%% RUN THE PSO ALGORITHM
for t=2:t_end
    itercount = itercount + 1;
    
    
    % calculate velocity
    velocity = phi*velocity...
         + gamma1*rand(size(velocity)).*(bestpars-swarm)...
         + gamma2*rand(size(velocity))...
         .*(ones(particleAmount,1)*bestpars(bestidx,:)-swarm);
     
    % limit velocities
    velocityLimitFraction = 1;
    for ind = 1:particleAmount
        for par = 1:parameterAmount
            % limit velocity
            velocity(ind,par) = max(min(velocity(ind,par),velocityLimitFraction*(limits(par,2)-limits(par,1))),velocityLimitFraction*(limits(par,1)-limits(par,2)));
        end
    end

    % update swarm locations
    swarm = swarm + velocity;
    
    % check boundaries
    boundcheck = zeros(size(swarm));
    for par = 1:parameterAmount
        boundcheck(:,par) = swarm(:,par) - max(min(swarm(:,par),limits(par,2)),limits(par,1));
    end
     
    % simulate
    parfor ind = 1:particleAmount 
        if all(boundcheck(ind,:)==0) % check that none of the values are over limits
            params = swarm(ind,:); 
            % round the integer values
            params(integerIndices) = round(params(integerIndices));
            [objects, constraints] = design_PMSM_generator(params);
            val = evalobjects(objects, constraints);
        else
            val = -1000000000; % punish for exceeding the limits
        end
        if val > bestresults(ind)
            bestpars(ind,:) = swarm(ind,:);
            bestresults(ind) = val;
        end
    end
    
    % update global best value
    [bestvalcand, bestidxcand] = max(bestresults);
    if(bestvalcand > bestval)
       bestval = bestvalcand;
       bestidx = bestidxcand;
       
    end
    bestvals(t) = bestval;
    disp(['Iteration ' num2str(t)]);
    disp(['Best value: ' num2str(bestval)]);

    % Save results after every 1000 iterations
	if (itercount == 1000)
	    params = bestpars(bestidx,:);
	    params(integerIndices) = round(params(integerIndices));
	    filename = strcat('output-',int2str(arrayTaskNumber));
        save(filename, 'bestval','params','bestvals');
	    disp(sprintf('Saved values.'));
		itercount = 0;
	end
end


%% SAVE RESULTS TO DISK
params = bestpars(bestidx,:);
params(integerIndices) = round(params(integerIndices));
filename = strcat('output-',int2str(arrayTaskNumber));
save(filename, 'bestval','params','bestvals');
disp(sprintf('SUCCESS array task number %d',arrayTaskNumber));
exit
