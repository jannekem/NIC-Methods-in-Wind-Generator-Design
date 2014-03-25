threshold = 1700;
limits = [20 80 1; 
    35000 65000 0; 
    2000000 6000000 0;
    0.001 0.05 0;
    0.8 5 0; 
    0.6 0.95 0; 
    21000 48000 0;
    1.3 1.6 0;
    1 3 1;
    0.8 0.99 0;
    0.25 0.75 0;
    0.75 0.95 0;
    0.01 0.1 0;
    0.01 0.1 0];

hf = figure('color','white');
hold on

xpar = 3;
ypar = 4;
zpar = 5;

thresholdIndices = find(bestValues > threshold);

x = parameters(thresholdIndices,xpar);
y = parameters(thresholdIndices,ypar);
z = parameters(thresholdIndices,zpar);
ht = scatter3(x,y,z);
set(ht,'XDataSource','x')
set(ht,'YDataSource','y')
set(ht,'ZDataSource','z')
axis([limits(xpar,1) limits(xpar,2) limits(ypar,1) limits(ypar,2) limits(zpar,1) limits(zpar,2)])
grid
view(3)
drawnow
