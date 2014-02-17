% this function is used to evaluate the parameters
% the higher the return value, the better system
function val = evalobjects(objects, constraints)
Pout_target = 3000000;

% check constraints
if(constraints(1)<=0 || constraints(2)<1 || constraints(3)>=100)
    val = -100;
    return;
end

val = abs(Pout_target-objects(1))/100000 + objects(2) ...
    - objects(3)/100000 + 10*objects(4) + 10*objects(5) ...
    - objects(6)/1000000;