%% 
clear all
bestvalues=[];  % best value of each run
parameters=[];  % best parameters of each run

values=[];      % the development of the values

iterations = 1000;
failures = 0;

for index = 1:iterations
    % read output from the jobs
    filename = strcat( 'genetic-output-', int2str( index ) );
    try
        load ( filename );
        parameters = [ parameters; elite ];
        bestvalues = [ bestvalues, max(elite_fitness) ];
        values = [values, elite_fitness'];
    catch
        disp ( sprintf ( 'FAILURE no file %s', filename ) );
        failures = failures + 1;
    end
    
   
end
mean_value = mean ( bestvalues )
variance = var ( bestvalues )
%hist ( bestvalues )

%% fill zeros
A = values;
l = A == 0;
ii = bsxfun(@plus,size(A,1)*(0:size(A,2)-1),cumsum(~l));
out = A;
out(l) = A(ii(l));
% remove all negative efficiencies
%values = max(0,values);
out = max(0,out);
values = out;

%% (ii) Compute mean consumption
M = length(bestvalues);
N = length(elite_fitness);
stdv=zeros(1,N);
mv = zeros(1,N);
for indn = 1:N
    mv(indn)=sum(values(indn,:))/M;
    stdv(indn)=std(values(indn,:));
end


figure(2);
% create custom grayscale colormap
customColormap=[];
for ind = 1:70
    val = 250 - ind*1;
    customColormap = [customColormap; val,val,val];
end
customColormap = customColormap./255;


set(gca,'NextPlot','replacechildren','ColorOrder',customColormap)
plot(values);
hold on;
plot(mv,'r','Linewidth',4,'Color', 'b');
plot(mv+stdv,'k','Linewidth',2.5);
plot(mv-stdv,'k','Linewidth',2.5);
hold off;
xlabel('iteration');ylabel('efficiency');
%axis([0 10 0 30])

%% PLOT best values

customColormap=[];
for ind = 1:50
    val = 250 - ind*5;
    customColormap = [customColormap; val,val,val];
end
customColormap = customColormap./255;

figure(3)
colormap(customColormap)
scatter(parameters(:,1),parameters(:,2),[],values(end,:))
%xlabel('Relative slot opening')
%ylabel('Relative slot width')


figure(4)
colormap(customColormap)
scatter(parameters(:,3),parameters(:,4),[],values(end,:))

figure(5)
colormap(customColormap)
scatter(parameters(:,5),parameters(:,6),[],values(end,:))

figure(6)
colormap(customColormap)
scatter(parameters(:,7),parameters(:,8),[],values(end,:))

figure(7)
colormap(customColormap)
scatter(parameters(:,9),parameters(:,10),[],values(end,:))

figure(8)
colormap(customColormap)
scatter(parameters(:,11),parameters(:,12),[],values(end,:))

figure(9)
colormap(customColormap)
scatter(parameters(:,13),parameters(:,14),[],values(end,:))


%%
figure(10)
hist(values(end,:))
