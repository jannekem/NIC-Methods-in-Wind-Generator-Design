function genetic_func(arrayTaskNumber)
addpath('model');
%{
Step 0 - Initializing some constants.
%}

% initialize random generator
rng(sum((arrayTaskNumber+100)*clock));

display = true;
population_size = 30;
param_size = 14;
survivals_size = round(population_size/3);
mutation_chance = 0.05;
crossover_chance = 0.7;
iter_count = 1;
iter_limit = 10000;
iter_converged = 20000;
iter_converged_count = 0;
condition_fitness = inf;

fitness = zeros(1,population_size);
population = zeros(population_size, param_size);
survivals = zeros(survivals_size, param_size);
elite = zeros(1,param_size);
elite_fitness = zeros(1,iter_limit);

% State space limitations for inputs:
% NOTE: Some parameters must be integers and thus will be treated different
% from the other parameters. Integer parameters are indicated by
% 'integer_indexes'

max_limits = [80,65000, 6000000, 0.05, 5, 0.95, 48000, 1.6, 3, 0.99, 0.75, 0.95, 0.1, 0.1];
min_limits = [20, 35000, 2000000, 0.001, 0.8, 0.6, 21000, 1.3, 1, 0.8, 0.25, 0.75, 0.01, 0.01];
integer_indexes = [1 9];

%%
%[objects, constraints] = design_PMSM_generator(params);

%{
Step 1 - Start: Generate random population.
%}

for i = 1:population_size
    for a = 1:param_size
       % Check wheter the parameter is an integer. If it is, round the
       % value to the nearest integer.
       if(sum(integer_indexes == a)~= 0)
          population(i,a) = round(rand*(max_limits(a) - min_limits(a)) + min_limits(a));  
       else
          population(i,a) = rand*(max_limits(a) - min_limits(a)) + min_limits(a);  
       end
    end
end

savecounter = 1;

%%
while(iter_count <= iter_limit)
    if(display)
        disp(['iteration: ' num2str(iter_count)])
    end
    %{
    Step 2 - Fitness: Evaluating the fitness of wind turbine design based on
    the given simulation model.
    %}
    
    % For finding elite for steps 3 and 4
    
    elite_idx = 0;
    elite_fitness(iter_count) = -inf;
    
    % Worst fitness is needed to handle negative fitness values later on.
    
    worst_fitness = inf;
    
    % Calculating fitness:
    
    for i = 1:population_size
        [objectives, constraints] = design_PMSM_generator(population(i,:)');
        fitness(i) = evalobjects(objectives,constraints);     % Fitness function
    end
    
    for i = 1:population_size
       if(fitness(i) > elite_fitness(iter_count))
            elite_idx = i;
            elite_fitness(iter_count) = fitness(i);
        end
        
        if(fitness(i) < worst_fitness)
            worst_fitness = fitness(i);
        end 
    end
    
    elite = population(elite_idx,:);
    avr_fitness(iter_count) = sum(fitness)/population_size;   % For performance evaluation
    fitness(elite_idx) = [];
    
    % Scaling negative fitness values to positive.
    
    fitness = fitness + abs(worst_fitness);
    
    %{
    Step 3 - Test: Test whether satisfactory solution is found
    %}
    
    if(elite_fitness > condition_fitness)
       if(display)
            disp(['Candidate with sufficent fitness of ' ...
             num2str(elite_fitness) ' found. Candidate index: ' num2str(elite_idx)])
       end
       break; 
    end    
    
    %{
    Step 4 - New population: Generate a new population using selection,
    crossover and mutation. Selection is based on "roulette wheel
    selection", where 1/3 of total population survives. Additionally Elitism is used, that is, the best candidate
    is passed on to the next generation as is.
    %}
    
    population(elite_idx,:) = [];
    fitness_sum = sum(fitness);
    
    % Roulette wheel selection:
    
    for i = 1:survivals_size
        winner = rand*fitness_sum;
        total = 0;
        for a = 1:population_size-1        
            total = total+fitness(a);
            if(total >= winner)
                survivals(i,:) = population(a,:);
                fitness_sum = fitness_sum - fitness(a);
                fitness(a) = [];
                population(a,:) = [];
                break;
            end
            
        end
    end
    
    population = zeros(population_size,param_size);
    population(1,:) = elite;    % Elitism
    
    % Create new population
    
    for i = 2:population_size
        
       % Select parents from the survivals
       
       parent_idx = randperm(survivals_size, 2);
       
       % Crossover
       if(rand < crossover_chance)      % Crossover happens
           cut_idx = ceil(rand*(param_size-1));
           
           population(i,1:cut_idx) = survivals(parent_idx(1),1:cut_idx);
           population(i,(cut_idx+1):end) = survivals(parent_idx(2),(cut_idx+1):end);   
       else                             % No crossover.
           population(i,:) = survivals(parent_idx(1),:);
           %{
           % Alternative solution (take average):
           for ind = 1:param_size
               population(i,ind) = (survivals(parent_idx(1),ind)+ ...
                   survivals(parent_idx(2),ind))/2;
               if(sum(integer_indexes == ind) ~= 0)
                  population(i,ind) = round(population(i,ind)); 
               end
           end
           %}
       end
       
       % Mutations
       mutations = rand(1,param_size);
       mutation_idx = [];
       for ind = 1:param_size
           if(mutations(ind) < mutation_chance)
               mutation_idx(end+1) = ind;
           end
       end
       
       num_mutations = length(mutation_idx);
       if(num_mutations >= 1)       % Mutation happens
           
           for ind = 1:num_mutations   
               param_idx = mutation_idx(ind);
               if(rand >= 0.5) % Positive amount
                   mutation = (max_limits(param_idx)-population(i,param_idx))*rand;
               else            % Negative amount
                   mutation = (min_limits(param_idx)-population(i,param_idx))*rand;
               end
               
               % If an integer parameter is mutated...
               if(sum(integer_indexes == param_idx) ~= 0)
                   mutation = round(mutation);
               end
               
               population(i,param_idx) = population(i,param_idx) + mutation;
           end         
       end
    end
    
    % Break iteration if elite fitness hasen't improved on 200
    % consecutive iterations.
    if(iter_count >= 2)
       if(elite_fitness(iter_count) == elite_fitness(iter_count-1))
           iter_converged_count = iter_converged_count+1; 
       else
           iter_converged_count = 0;
       end
    end
    if(iter_converged_count >= iter_converged)
        if(display)
            disp(['Elite fitness has not improved on ' num2str(iter_converged) ...
                'consecutive iterations. breaking at iteration ' num2str(iter_count)])
        end
        break;
    end
    
    iter_count = iter_count + 1;
	
	savecounter = savecounter + 1;
	if (savecounter == 500)
		filename = strcat('genetic-output-',int2str(arrayTaskNumber));
		save(filename, 'elite_fitness','avr_fitness','elite','iter_count');
		savecounter = 0;
	end
end

if(display)
    figure(1);
    plot(elite_fitness);
    xlabel('Iter');
    ylabel('Fitness');
    title('Elite fitness development');

    figure(2);
    plot(avr_fitness);
    xlabel('Iter');
    ylabel('Fitness');
    title('Average fitness development');
end

%% SAVE RESULTS TO DISK
filename = strcat('genetic-output-',int2str(arrayTaskNumber));
save(filename, 'elite_fitness','avr_fitness','elite','iter_count');
disp(sprintf('SUCCESS array task number %d',arrayTaskNumber));
exit
