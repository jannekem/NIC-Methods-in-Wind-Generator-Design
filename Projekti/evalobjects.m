% this function is used to evaluate the parameters
% the higher the return value, the better system
function val = evalobjects(objects, constraints)
Pout_target = 3000000;

% evaluate fitness function
val = - abs(Pout_target-objects(1))+ objects(2)*100 - objects(4)/1000+ objects(5)*1000;

% punish for exceeding the constraints
if(constraints(1)<=0 || constraints(2)<1 || constraints(3)>=100)
    val = val - 10000;
end


%val = 1000/abs(Pout_target-objects(1)) + objects(2) ...
%     - objects(3)/100000 + 10*objects(4) + 10*objects(5) ...
%     - objects(6)/100000;