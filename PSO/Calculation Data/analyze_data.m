clear all
bestvals=[];
parameters=[];
iterations = 200;

for index = 1:iterations
    % read output from the jobs
    filename = strcat( 'output-', int2str( index ) );
    try
        load ( filename );
        parameters = [ parameters; params ];
        bestvals = [ bestvals, bestval ];
    catch
        disp ( sprintf ( 'FAILURE no file %s', filename ) );
    end
    
   
end
mean_value = mean ( bestvals )
variance = var ( bestvals )
hist ( bestvals )

