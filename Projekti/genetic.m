%{
Step 0 - Initializing some constants.
%}
population_size = 20;
param_size = 14;
survivals_size = round(population_size/3);
mutation_chance = 0.05;
crossover_chance = 0.7;
iter_count = 1;
iter_limit = 15;
condition_fitness = inf;
global_best = -inf; % For performance evaluation

fitness = zeros(1,population_size);
population = zeros(population_size, param_size);
survivals = zeros(survivals_size, param_size);
elite = zeros(1,param_size);

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

%%
while(iter_count < iter_limit)
    disp(['iteration: ' num2str(iter_count)])
    %{
    Step 2 - Fitness: Evaluating the fitness of wind turbine design based on
    the given simulation model.
    %}
    
    % For finding elite for steps 3 and 4
    
    elite_idx = 0;
    elite_fitness = -inf;
    
    % Worst fitness is needed to handle negative fitness values later on.
    
    worst_fitness = 0;
    
    % Calculating fitness:
    
    for i = 1:population_size
        [objectives, constraints] = design_PMSM_generator(population(i,:)');
        fitness(i) = evalobjects(objectives,constraints);     % Fitness function
        if(fitness(i) > elite_fitness)
            elite_idx = i;
            elite_fitness = fitness(i);
        end
        
        if(fitness(i) < worst_fitness)
            worst_fitness = fitness(i);
        end
    end
    
    elite = population(elite_idx,:);
    
    % Scaling negative fitness values to positive.
    
    fitness = fitness + abs(worst_fitness);
    
    %{
    Step 3 - Test: Test whether satisfactory solution is found
    %}
    
    if(elite_fitness > global_best)
        global_best = elite_fitness;
    end
    
    if(elite_fitness > condition_fitness)
       disp(['Candidate with sufficent fitness of ' ...
           num2str(elite_fitness) ' found. Candidate index: ' num2str(elite_idx)])
       break; 
    end    
    
    %{
    Step 4 - New population: Generate a new population using selection,
    crossover and mutation. Selection is based on "roulette wheel
    selection", where 1/3 of total population survives. Additionally Elitism is used, that is, the best candidate
    is passed on to the next generation as is.
    %}
    
    % Roulette wheel selection:
    
    fitness_sum = sum(fitness);
    avr_fitness(iter_count) = fitness_sum/population_size;   % For performance evaluation
    
    for i = 1:survivals_size
        winner = rand*fitness_sum;
        total = 0;
        len = size(population,1);
        for a = 1:len
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
    
    % Elitism:
    
    population(1,:) = elite;
        
    % Crossover and mutation
    
    for i = 2:population_size
       if(rand < crossover_chance)      % Crossover happens
           cut_idx = ceil(rand*(param_size-1));
           survival_idx1 = ceil(rand*survivals_size);
           survival_idx2 = ceil(rand*survivals_size);
           
           population(i,1:cut_idx) = survivals(survival_idx1,1:cut_idx);
           population(i,(cut_idx+1):end) = survivals(survival_idx2,(cut_idx+1):end);   
       else                             % No crossover
           survival_idx1 = ceil(rand*survivals_size);
           population(i,:) = survivals(survival_idx1,:);
       end
       
       if(rand < mutation_chance)       % Mutation happens
           mutation_idx = ceil(rand*param_size);
           range = max_limits(mutation_idx) - min_limits(mutation_idx);
           mutation = (rand*2 -1)*range*0.1;
           
           % If an integer parameter is mutated...
           if(sum(integer_indexes == mutation_idx) ~= 0)
              if(mutation < 0)
                  mutation = floor(mutation);
              else
                  mutation = ceil(mutation);
              end
           end
           
           if(population(i,mutation_idx)+mutation < min_limits(mutation_idx))
               population(i,mutation_idx) = population(i,mutation_idx) - mutation;
               continue;
           end
           
           if(population(i,mutation_idx)+mutation > max_limits(mutation_idx))
               population(i,mutation_idx) = population(i,mutation_idx) - mutation;
               continue;
           end
           
           population(i,mutation_idx) = population(i,mutation_idx) + mutation;
           
       end
    end
    
    iter_count = iter_count + 1;
end





