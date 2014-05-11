clear all
PSObestvalues=[];  % best value of each run
PSOparameters=[];  % best parameters of each run
PSOvalues=[];    % the development of the values

GAbestvalues=[];
GAparameters=[];
GAvalues=[];

iterations = 1000;
failures = 0;

for index = 1:iterations
    % read output from the jobs
    filename = strcat( 'genetic-output-', int2str( index ) );
    try
        load ( filename );
        GAparameters = [ GAparameters; elite ];
        GAbestvalues = [ GAbestvalues, max(elite_fitness) ];
        GAvalues = [GAvalues, elite_fitness'];
    catch
        disp ( sprintf ( 'FAILURE no file %s', filename ) );
        failures = failures + 1;
    end
    
   
end

for index = 1:iterations
    % read output from the jobs
    filename = strcat( 'output-', int2str( index ) );
    try
        load ( filename );
        PSOparameters = [ PSOparameters; params ];
        PSObestvalues = [ PSObestvalues, bestval ];
        PSOvalues = [PSOvalues, bestvals];
    catch
        disp ( sprintf ( 'FAILURE no file %s', filename ) );
        failures = failures + 1;
    end
    
   
end

GAmean_value = mean ( GAbestvalues );
GAvariance = var ( GAbestvalues );

PSOmean_value = mean ( PSObestvalues );
PSOvariance = var ( PSObestvalues );
%hist ( bestvalues )

%% fill zeros
A = GAvalues;
l = A == 0;
ii = bsxfun(@plus,size(A,1)*(0:size(A,2)-1),cumsum(~l));
out = A;
out(l) = A(ii(l));
% remove all negative efficiencies
%values = max(0,values);
out = max(0,out);
GAvalues = out;

%% (ii) Compute mean consumption
M = length(GAbestvalues);
N = length(elite_fitness);
stdv=zeros(1,N);
mv = zeros(1,N);
for indn = 1:N
    GAmv(indn)=sum(GAvalues(indn,:))/M;
    GAstdv(indn)=std(GAvalues(indn,:));
end
%%
M = length(PSObestvalues);
N = length(bestvals);
PSOstdv=zeros(1,N);
PSOmv = zeros(1,N);
for indn = 1:N
    PSOmv(indn)=sum(PSOvalues(indn,:))/M;
    PSOstdv(indn)=std(PSOvalues(indn,:));
end

%%
figure(2);
hold on;
plot(GAmv,'r');
plot(GAmv+GAstdv,'r--');
plot(PSOmv,'b');
plot(PSOmv+PSOstdv,'b--');
legend('GA','GA std', 'PSO','PSO std');

plot(GAmv-GAstdv,'r--');
plot(PSOmv-PSOstdv,'b--');

hold off;
title('Algorithm comparison');
xlabel('iteration');ylabel('efficiency');
axis([0 10000 0.98 0.99])