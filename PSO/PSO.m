%% INITIALIZE THE SWARM AND VARIABLES
clear all
t_end = 100;
particleAmount = 50;
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
    val = -swarm(1,ind,1)^2-swarm(1,ind,2)^2;
    bestpars(ind,:) = swarm(1,ind,:);
    bestresults(ind) = val;
end

[bestval, bestidx] = max(bestresults);

%% 