clear all
bestvals=[];
parameters=[];
iterations = 100;

for index = 1:iterations
    % read output from the jobs
    filename = strcat( 'output-', int2str( index ) );
    try
        load ( filename );
    catch
        disp ( sprintf ( 'FAILURE no file %s', filename ) );
    end
    parameters = [ parameters; params ];
    bestvals = [ bestvals, bestval ];
end