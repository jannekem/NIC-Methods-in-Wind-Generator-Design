%initialize parameters
clear all
t_end = 50;
particleAmount = 100;
limits = [35000 65000; 2000000 6000000;
    0.001 0.05; 0.8 5; 0.6 0.95; 21000 48000;
    1.3 1.6; 0.8 0.99; 0.25 0.75;
    0.75 0.95; 0.01 0.1; 0.01 0.1]; % parameter limits
discrete_limits = [20 80; 1 3]; % indices for integer parameters

% generate the swarm and velocity matrices
swarm = zeros(t_end, particleAmount, 12);
swarmd = zeros(t_end, particleAmount, 2);
velocity = zeros(particleAmount, 12);
velocityd = zeros(particleAmount, 2);
bestpars = zeros(particleAmount, 12); % best parameters for each particle 
bestparsd = zeros(particleAmount, 2); % best discrete parameters ""
bestresults = zeros(particleAmount,1);
bestidx = 0;        % holder for the global best index
bestval = 0;

% initialize swarm
for ind = 1:particleAmount 
    for par = 1:12
       swarm(1,ind,par) = limits(par,1)+rand()*(limits(par,2)-limits(par,1));
    end
    swarmd(1,ind,1) = 50;% randi([20 80],1);
    swarmd(1,ind,2) = 2; %randi([1 3],1);
end

parfor ind = 1:particleAmount
    params = combineparams(swarm(1,ind,:), swarmd(1,ind,:));
    [objects, constraints] = design_PMSM_generator(params);
    bestpars(ind,:) = swarm(1,ind,:);
    bestparsd(ind,:) = swarmd(1,ind,:);
    bestresults(ind) = evalobjects(objects, constraints);
end

[bestval, bestidx] = max(bestresults);

hf = figure('color','white');
hold on
x = swarm(1,:,1);
y = swarm(1,:,2);
ht = scatter(x,y);
set(ht,'XDataSource','x')
set(ht,'YDataSource','y')
drawnow

% fig = figure('Visible','off');
% x = swarm(1,:,1);
% y = swarm(1,:,2);
% scatter(x,y);

% % save video frame
% aviobj = avifile('out.avi','FPS',1);
% img = hardcopy(fig, '-dzbuffer', '-r0'); 
% aviobj = addframe(aviobj, im2frame(img));

paratemp = zeros(particleAmount, 14); % DEBUG VARIABLE FOR FINDING TIME INTENSIVE PARAMETER COMBINATIONS

    meas_t = zeros(t_end,1);
for t=2:t_end
    swarm_temp = squeeze(swarm(t-1,:,:)); % required for matrix subtraction
    % create parameter validity matrix to ignore bad personal bests
    validity = bestresults;
    validity(validity==-100) = 0;
    validity(validity~=0) = 1;
    
        
    % calculate velocity
    
    velocity = velocity...
        + 2*rand(size(velocity)).*(bestpars-swarm_temp)...
        + 2*rand(size(velocity))...
        .*(ones(particleAmount,1)*bestpars(bestidx,:)-swarm_temp);
    
    % update swarm locations
    swarm(t,:,:) = swarm_temp + velocity;
    
    % check boundaries
    boundcheck = zeros(size(swarm_temp));
    for par = 1:12
        boundcheck(:,par) = swarm(t,:,par) - max(min(swarm(t,:,par),limits(par,2)),limits(par,1));
    end
    boundcheck(boundcheck~=0) = 1;
    
    % simulate
    parfor ind = 1:particleAmount 
        if  sum(boundcheck(ind,:))==0 % check that none of the values is over limits
            params = combineparams(swarm(t,ind,:), swarmd(1,ind,:)); % TODO: discrete value handling
            [objects, constraints] = design_PMSM_generator(params);
            result = evalobjects(objects, constraints);
        else % if values were over limits, don't simulate
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
    x = swarm(t,:,1);
    y = swarm(t,:,2);
    %scatter(x,y);
    
%     % save to video
%     img = hardcopy(fig, '-dzbuffer', '-r0'); 
%     aviobj = addframe(aviobj, im2frame(img));
    t
    refreshdata(hf,'caller');
    axis([limits(1,1) limits(1,2) limits(2,1) limits(2,2)])
    drawnow
end

params = combineparams(bestpars(bestidx,:)', bestparsd(bestidx,:));
[objects, constraints] = design_PMSM_generator(params)

% aviobj = close(aviobj);