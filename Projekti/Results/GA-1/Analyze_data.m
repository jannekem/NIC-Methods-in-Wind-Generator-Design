clear all
elite_fitnesses=[];
iter_counts=[];
elites=[];
avr_fitnesses=[];
iterations = 200;

for index = 1:iterations
    % read output from the jobs
    filename = strcat( 'genetic-output-', int2str( index ) );
    try
        load ( filename );
        elite_fitnesses(end+1) = elite_fitness;
        elites(end+1) = elite;
        avr_fitnesses(end+1) = avr_fitness;
        elite_fitnesses(end+1) = elite_fitness;
    catch
        disp ( sprintf ( 'FAILURE no file %s', filename ) );
    end
    
   
end
mean_value = mean ( bestvals )
variance = var ( bestvals )
hist ( bestvals )

