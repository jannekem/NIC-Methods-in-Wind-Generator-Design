%% INITIALIZE THE SWARM AND VARIABLES
clear all
t_end = 100;
particleAmount = 50;
phi = 1; % inertia coefficient
gamma1 = 2;% acceleration constants
gamma2 = 2; 

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
validity = zeros(particleAmount,1);     % holder for the personal validities

% generate the swarm and velocity matrices
swarm = zeros(t_end, particleAmount, parameterAmount);
velocity = zeros(particleAmount, parameterAmount);
bestpars = zeros(particleAmount, parameterAmount);  % best parameters for each particle 
bestresults = zeros(particleAmount,1);
bestidx = 0;        % holder for the global best index
bestval = 0;        % global best value
bestvals = zeros(t_end,1); % best value vector

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
bestvals(1) = bestval;

%% RUN THE PSO ALGORITHM
for t=2:t_end
    swarm_temp = squeeze(swarm(t-1,:,:)); % required for matrix subtraction  
    
    % calculate velocity
    velocity = (0.5+rand()/2)*velocity...
         + gamma1*rand(size(velocity)).*(bestpars-swarm_temp)...
         + gamma2*rand(size(velocity))...
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
    bestvals(t) = bestval;
    disp(['Iteration ' num2str(t)]);
    disp(['Best value: ' num2str(bestval)]);
end

%% PLOT FIGURES

hf = figure('color','white');
hold on

xpar = 3;
ypar = 4;
zpar = 5;

x = swarm(1,:,xpar);
y = swarm(1,:,ypar);
z = swarm(1,:,zpar);
ht = scatter3(x,y,z);
set(ht,'XDataSource','x')
set(ht,'YDataSource','y')
set(ht,'ZDataSource','z')
grid
drawnow

for t = 1:t_end
    pause(0.05)
    grid
    x = swarm(t,:,xpar);
    y = swarm(t,:,ypar);
    z = swarm(t,:,zpar);
    refreshdata(hf,'caller');
    axis([limits(xpar,1) limits(xpar,2) limits(ypar,1) limits(ypar,2) limits(zpar,1) limits(zpar,2)])
    grid
    view(3)
    drawnow
end
% figure
% gplotmatrix(squeeze(swarm(t_end,:,:)),[],swarm(t_end,:,3))

%% PRINT RESULTS TO CONSOLE

params = bestpars(bestidx,:);
params(integerIndices)=round(params(integerIndices));
[objects, constraints] = design_PMSM_generator(params);
disp(['Best parameters: ', num2str(params(1)) ', ' num2str(params(2)) ', ' num2str(params(3)) ', ' num2str(params(4)) ...
    ', ' num2str(params(5)) ', ' num2str(params(6)) ...
    ', ' num2str(params(7)) ', ' num2str(params(8)) ...
    ', ' num2str(params(9)) ', ' num2str(params(10)) ...
    ', ' num2str(params(11)) ', ' num2str(params(12)) ...
    ', ' num2str(params(13)) ', ' num2str(params(14))]);
disp(['Constraints: ', num2str(constraints(1)), ' ' num2str(constraints(2)), ' ' num2str(constraints(3))]);





