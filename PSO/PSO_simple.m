%initialize parameters
clear all
t_end = 100;
particleAmount = 50;
parameterAmount = 2;
limits = [-1 1; -1 1]; % parameter limits


% generate the swarm and velocity matrices
swarm = zeros(t_end, particleAmount, parameterAmount);
swarmd = zeros(t_end, particleAmount, 2);
velocity = zeros(particleAmount, parameterAmount);
velocityd = zeros(particleAmount, 2);
bestpars = zeros(particleAmount, parameterAmount); % best parameters for each particle 
bestparsd = zeros(particleAmount, 2); % best discrete parameters ""
bestresults = zeros(particleAmount,1);
bestidx = 0;        % holder for the global best index
bestval = 0;

% initialize swarm
for ind = 1:particleAmount 
    for par = 1:parameterAmount
       swarm(1,ind,par) = limits(par,1)+rand()*(limits(par,2)-limits(par,1));
    end
end

for ind = 1:particleAmount
    params = swarm(1,ind,:);
    val = -swarm(1,ind,1)^2-swarm(1,ind,2)^2;
    bestpars(ind,:) = swarm(1,ind,:);
    bestparsd(ind,:) = swarmd(1,ind,:);
    bestresults(ind) = val;
end

[bestval, bestidx] = max(bestresults);




paratemp = zeros(particleAmount, 14); % DEBUG VARIABLE FOR FINDING TIME INTENSIVE PARAMETER COMBINATIONS

    meas_t = zeros(t_end,1);
for t=2:t_end
    
    swarm_temp = squeeze(swarm(t-1,:,:)); % required for matrix subtraction
    % create parameter validity matrix to ignore bad personal bests
    validity = bestresults;
    validity(validity(:,1)==-100) = 0;
    validity(validity(:,1)~=0) = 1;
    sparvalidity = sparse(diag(validity));
        
    % calculate velocity
    
%     velocity = velocity...
%         + 2*rand(size(velocity)).*(bestpars-swarm_temp)...
%         + 2*rand(size(velocity))...
%         .*(ones(particleAmount,1)*bestpars(bestidx,:)-swarm_temp)
    velocityLimitFraction = 0.01;
    for ind = 1:particleAmount
        for par = 1:parameterAmount
            velocity(ind,par) = velocity(ind,par)+2*rand()*(bestpars(ind,par)-swarm_temp(ind,par))+2*rand()*(bestpars(bestidx,par)-swarm_temp(ind,par));
            % limit velocity
            velocity(ind,par) = max(min(velocity(ind,par),velocityLimitFraction*(limits(par,2)-limits(par,1))),velocityLimitFraction*(limits(par,1)-limits(par,2)));
        end
    end

    % update swarm locations
    
    swarm(t,:,:) = swarm_temp + velocity;
    
    % check boundaries
    boundcheck = zeros(size(swarm_temp));
    boundcheck(:,par) = swarm(t,:,par) - max(min(swarm(t,:,par),limits(par,2)),limits(par,1));
 
    %boundcheck(boundcheck~=0) = 1;
 
    % simulate
    for ind = 1:particleAmount 
        if size(boundcheck(ind)==0) % check that none of the values is over limits
            params = swarm(t,ind,:); % TODO: discrete value handling
            result = -swarm(t,ind,1)^2-swarm(t,ind,2)^2;
        else
            result = -200;
        end
        if result > bestresults(ind)
            bestpars(ind,:) = swarm(t,ind,:);
            bestparsd(ind,:) = swarmd(1,ind,:);
            bestresults(ind) = result;
        end
    end
    
    % update global best value
    [bestvalcand, bestidxcand] = max(bestresults);
    if(bestvalcand > bestval)
       bestval = bestvalcand;
       bestid = bestidxcand;
    end
    
   
    
    t
end

% % display figure
hf = figure('color','white');
hold on
x = swarm(1,:,1);
y = swarm(1,:,2);
ht = scatter(x,y);
set(ht,'XDataSource','x')
set(ht,'YDataSource','y')
drawnow


for t = 1:t_end
    pause(0.01)
    x = swarm(t,:,1);
    y = swarm(t,:,2);
    refreshdata(hf,'caller');
    axis([limits(1,1) limits(1,2) limits(2,1) limits(2,2)])
    drawnow
end


%params = combineparams(bestpars(bestidx,:)', bestparsd(bestidx,:));
%[objects, constraints] = design_PMSM_generator(params)

