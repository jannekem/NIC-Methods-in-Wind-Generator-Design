%% INITIALIZE THE SWARM AND VARIABLES
clear all
t_end = 40;
particleAmount = 40;
parameterAmount = 14;

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
integerIndices = find(limits(:,3)==1); % find the indices for the integer values
validity = zeros(particleAmount,1); % holder for the personal validities

% generate the swarm and velocity matrices
swarm = zeros(t_end, particleAmount, parameterAmount);
velocity = zeros(particleAmount, parameterAmount);
bestpars = zeros(particleAmount, parameterAmount); % best parameters for each particle 
bestresults = zeros(particleAmount,1);
bestidx = 0;        % holder for the global best index
bestval = 0;        % global best value

% initialize swarm
for ind = 1:particleAmount 
    for par = 1:parameterAmount
       swarm(1,ind,par) = limits(par,1)+rand()*(limits(par,2)-limits(par,1));
    end
end

%% EVALUATE THE INITIAL VALUES
for ind = 1:particleAmount
    params = squeeze(swarm(1,ind,:));
    params(integerIndices) = round(params(integerIndices)); % round to the nearest integer
    [objects, constraints] = design_PMSM_generator(params);
    bestpars(ind,:) = swarm(1,ind,:);
    val = evalobjects(objects, constraints);
    bestresults(ind) = val;
end

[bestval, bestidx] = max(bestresults);  % evaluate the initial best results

%% RUN THE PSO ALGORITHM
for t=2:t_end
    swarm_temp = squeeze(swarm(t-1,:,:)); % required for matrix subtraction  
    
    % calculate velocity
    velocity = velocity...
         + 2*rand(size(velocity)).*(bestpars-swarm_temp)...
         + 2*rand(size(velocity))...
         .*(ones(particleAmount,1)*bestpars(bestidx,:)-swarm_temp);
     
    % limit velocities
    velocityLimitFraction = 0.1;
    for ind = 1:particleAmount
        for par = 1:parameterAmount
            % limit velocity
            velocity(ind,par) = max(min(velocity(ind,par),velocityLimitFraction*(limits(par,2)-limits(par,1))),velocityLimitFraction*(limits(par,1)-limits(par,2)));
        end
    end

    % update swarm locations
    swarm(t,:,:) = swarm_temp + velocity;
    
    % check boundaries
    boundcheck = zeros(size(swarm_temp));
    for par = 1:parameterAmount
        boundcheck(:,par) = swarm(t,:,par) - max(min(swarm(t,:,par),limits(par,2)),limits(par,1));
    end
     
    % simulate
    parfor ind = 1:particleAmount 
        if size(boundcheck(ind)==0) % check that none of the values are over limits
            params = squeeze(swarm(t,ind,:)); 
            % round the integer values
            params(integerIndices) = round(params(integerIndices));
            [objects, constraints] = design_PMSM_generator(params);
            val = evalobjects(objects, constraints);
        else
            val = -1000000000; % punish for exceeding the limits
        end
        if val > bestresults(ind)
            bestpars(ind,:) = swarm(t,ind,:);
            bestresults(ind) = val;
        end
    end
    
    % update global best value
    [bestvalcand, bestidxcand] = max(bestresults);
    if(bestvalcand > bestval)
       bestval = bestvalcand;
       bestidx = bestidxcand;
    end
end

%% PLOT FIGURES

hf = figure('color','white');
hold on

xpar = 1;
ypar = 2;

x = swarm(1,:,xpar);
y = swarm(1,:,ypar);
ht = scatter(x,y);
set(ht,'XDataSource','x')
set(ht,'YDataSource','y')
drawnow

for t = 1:t_end
    pause(0.1)
    x = swarm(t,:,xpar);
    y = swarm(t,:,ypar);
    refreshdata(hf,'caller');
    axis([limits(xpar,1) limits(xpar,2) limits(ypar,1) limits(ypar,2)])
    drawnow
end
