function [] = spring()
clear
clc
close all;

%variables related to the springs
spring.k = 1; %k constant
spring.x0 = 1; %equilibrium distance of each spring

%variables related to the simulation
sim.dt = 1; %time step
sim.delay = 0.001; %time in seconds between executing steps
sim.steps = 10000; %number of iterations to perform

%variables related to the particles
particles.N = 10; %Number of particles
particles.m = 1; %mass of each particle
particles.fixed = boolean(zeros(1,particles.N)); %logical array of which particles are fixed
particles.fixed(1) = true;
particles.fixed(particles.N) = true;
%particles.fixed(5) = true;
particles.indexToExamine = 3;%index of particle to exam in detail (1 - N)

%initialize other variables and create figure for plotting
figure('units','normalized','outerposition',[0 0 1 1]) 
[particles] = generateInitialPositions(particles);
std = zeros(1,sim.steps);
examine.x = zeros(1,sim.steps);
examine.v = zeros(1,sim.steps);
examine.a = zeros(1,sim.steps);

%loop through for all avaiable steps
for n = 1:10000
    [particles] = getA(particles,spring); %calculate the acceleration
    std(n) = sqrt(sum((diff(particles.x) - spring.x0).^2)./(particles.N-1)); %determine the standard deviation of the particle locations
    examine = examineSingleParticle(particles,examine,n); %add the current values to the array storage for the particle being examined
    particles = getNewPositions(particles,sim); %update the particle positions
    plotPositions(particles,examine,std,n); %plot results
    pause(sim.delay); %simulation pause
end
end

function [examine] = examineSingleParticle(particles,examine,currStep)
examine.x(currStep) = particles.x(particles.indexToExamine);
examine.v(currStep) = particles.v(particles.indexToExamine);
examine.a(currStep) = particles.a(particles.indexToExamine);
end

function [particles] = generateInitialPositions(particles)
particles.x = 1:particles.N;
particles.v = zeros(1,particles.N);
particles.a = zeros(1,particles.N);
particles.v(3) = 0.5;
particles.v(particles.fixed) = 0;
end

function [] = plotPositions(particles,examine,std,currStep)
axisVec = [min(particles.x)-1,max(particles.x)+1,-0.5,0.5];
y = zeros(1,particles.N);

h = subplot(4,2,1);
cla(h)
plot(particles.x(particles.fixed),y(particles.fixed),'sb');
hold on
plot(particles.x(~particles.fixed),y(~particles.fixed),'*b');
quiver([particles.x(particles.indexToExamine),particles.x(particles.indexToExamine)],[0.5,-0.5],[0,0],[-0.4,0.4],0,'b');
grid on
axis(axisVec)
title('Position of Particles');
xlabel('x');

h = subplot(4,2,3);
cla(h)
quiver(particles.x,y,particles.v,y,0.25);
hold on
quiver([particles.x(particles.indexToExamine),particles.x(particles.indexToExamine)],[0.5,-0.5],[0,0],[-0.4,0.4],0,'b');
grid on
axis(axisVec)
title('Velocity of Each Particle');
xlabel('x');

h = subplot(4,2,5);
cla(h)
quiver(particles.x,y,particles.a,y,0.25);
hold on
quiver([particles.x(particles.indexToExamine),particles.x(particles.indexToExamine)],[0.5,-0.5],[0,0],[-0.4,0.4],0,'b');
grid on
axis(axisVec)
title('Acceleration of Each Particle');
xlabel('x');

subplot(4,2,2)
plot(1:currStep,examine.x(1:currStep));
grid on
title(['Position of Particle ',num2str(particles.indexToExamine)]);
xlabel('time step');
ylabel('Position');

subplot(4,2,4)
plot(1:currStep,examine.v(1:currStep));
grid on
title(['Velocity of Particle ',num2str(particles.indexToExamine)]);
xlabel('time step');
ylabel('Velocity');

subplot(4,2,6)
plot(1:currStep,examine.a(1:currStep));
grid on
title(['Acceleration of Particle ',num2str(particles.indexToExamine)]);
xlabel('time step');
ylabel('Force');

subplot(4,2,8)
plot(1:currStep,std(1:currStep));
grid on
title('Standard Deviation of Particle Positions From Stabilization');
xlabel('time step');
ylabel('Standard Deviation');
end

function [particles] = getA(particles,spring)
particles.a = zeros(1,particles.N);
for n = 1:particles.N
    if n == 1
        particles.a(n) = calcForce(particles.x(n),particles.x(n+1),spring);
    elseif n == particles.N
        particles.a(n) = calcForce(particles.x(n),particles.x(n-1),spring);
    else
        particles.a(n) = calcForce(particles.x(n),particles.x(n-1),spring) + calcForce(particles.x(n),particles.x(n+1),spring);
    end
end
particles.a(particles.fixed) = 0;
particles.a = particles.a / particles.m;
end


function out = calcForce(x,x2,spring)
if(x < x2)
    tmp = x2 - spring.x0;
else
    tmp = x2 + spring.x0;
end
out = spring.k*(tmp - x)*abs(tmp - x);
end

function [particles] = getNewPositions(particles,sim)
dv = particles.a * sim.dt;
particles.v = particles.v + dv;
dx = particles.v * sim.dt + particles.a * sim.dt^2 / 2;
particles.x = particles.x + dx;
end