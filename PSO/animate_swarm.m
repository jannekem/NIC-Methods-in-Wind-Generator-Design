function val = animate_swarm(swarm, limits)
hf = figure('color','white');
hold on

xpar = 3;
ypar = 4;
zpar = 5;

x = swarm(1,:,xpar);
y = swarm(1,:,ypar);
z = swarm(1,:,zpar);
ht = scatter3(x,y,z);
set(ht,'XDataSource','x')
set(ht,'YDataSource','y')
set(ht,'ZDataSource','z')
grid
drawnow

for t = 1:t_end
    pause(0.05)
    grid
    x = swarm(t,:,xpar);
    y = swarm(t,:,ypar);
    z = swarm(t,:,zpar);
    refreshdata(hf,'caller');
    axis([limits(xpar,1) limits(xpar,2) limits(ypar,1) limits(ypar,2) limits(zpar,1) limits(zpar,2)])
    grid
    view(3)
    drawnow
end
