%% 

timerDelay = 18000;   % set the run time for the script
count = 0;
stat=true;
ttimer = timer('TimerFcn','stat=false; disp(''STOP!'')','StartDelay',timerDelay);
start(ttimer);

parameters = zeros(1,14);
bestValues = zeros(1,1);
while stat==true
    count = count + 1;
    [params, bestvals] = PSO_function(150,40,0.7,2,2);
    parameters(count, :)=params;
    bestValues(count) = bestvals(end);
end

delete(ttimer)
clear ttimer

save('savedparameters.mat');