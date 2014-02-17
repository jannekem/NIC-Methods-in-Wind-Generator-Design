%initialize parameters
t_end = 10;
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
    swarmd(1,ind,1) = randi([20 80],1);
    swarmd(1,ind,2) = randi([1 3],1);
end
tic
parfor ind = 1:particleAmount
    params = combineparams(swarm(1,ind,:), swarmd(1,ind,:));
    [objects, constraints] = design_PMSM_generator(params);
    bestpars(ind,:) = swarm(1,ind,:);
    bestparsd(ind,:) = swarmd(1,ind,:);
    bestresults(ind) = evalobjects(objects, constraints);
end
toc
[bestval, bestidx] = max(bestresults);

hf = figure('color','white');
hold on
x = swarm(1,:,1);
y = swarm(1,:,2);
ht = scatter(x,y);
set(ht,'XDataSource','x')
set(ht,'YDataSource','y')
drawnow

    % slice limits to reduce parfor-function overhead
    limits_slice1 = limits(:,1);
    limits_slice2 = limits(:,2);

    meas_t = zeros(t_end,1);
for t=2:t_end
    swarm_temp = squeeze(swarm(t-1,:,:)); % required for matrix subtraction
    % create parameter validity matrix to ignore bad personal bests
    validity = bestresults;
    validity(validity(:,1)==-100) = 0;
    validity(validity(:,1)~=0) = 1;
        
    % calculate velocity
    
    velocity = velocity...
        + diag(validity)*2*rand(size(velocity)).*(bestpars-swarm_temp)...
        + 2*rand(size(velocity))...
        .*(ones(particleAmount,1)*bestpars(bestidx,:)-swarm_temp);
   
    % update swarm locations
    swarm(t,:,:) = swarm_temp + velocity;
    
    % enforce boundaries
    
    for par = 1:12
        swarm(t,:,par) = max(min(swarm(t,:,par),limits_slice2(par)),limits_slice1(par));
    end
    tic
    % simulate
    for ind = 1:particleAmount 
        params = combineparams(swarm(t,ind,:), swarmd(1,ind,:)); % TODO: discrete value handling
        [objects, constraints] = design_PMSM_generator(params);
        result = evalobjects(objects, constraints);
        if result > bestresults(ind)
            bestpars(ind,:) = swarm(t,ind,:);
            bestparsd(ind,:) = swarmd(1,ind,:);
            bestresults(ind) = result;
        end
    end
    toc
    % update global best value
    [bestvalcand, bestidxcand] = max(bestresults);
    if(bestvalcand > bestval)
       bestval = bestvalcand;
       bestid = bestidxcand;
    end
    bestval
    t
    x = swarm(t,:,1);
    y = swarm(t,:,2);
    refreshdata(hf,'caller');
    drawnow
end

params = combineparams(bestpars(bestidx,:)', bestparsd(bestidx,:));
[objects, constraints] = design_PMSM_generator(params)


