clc
clear all 
close all

% include folder containing the TensegrityStucture class:
addpath('../tensegrityObjects')

% Structure dimensions
barLength = 1;            %[m]
barSpacing = barLength/2;   %[m]
barRadius = 0.010;          %[m]
stringRadius = 0.005;       %[m]
CoMheight = 0.5;            %[m]

nodes = [-barSpacing     barLength*0.5  CoMheight;
         -barSpacing    -barLength*0.5  CoMheight;
          barSpacing     barLength*0.5  CoMheight;
          barSpacing    -barLength*0.5  CoMheight;
          0             -barSpacing     barLength*0.5+CoMheight;
          0             -barSpacing    -barLength*0.5+CoMheight;
          0              barSpacing     barLength*0.5+CoMheight;
          0              barSpacing    -barLength*0.5+CoMheight;        
          barLength*0.5  0             -barSpacing+CoMheight;
         -barLength*0.5  0             -barSpacing+CoMheight;
          barLength*0.5  0              barSpacing+CoMheight;
         -barLength*0.5  0              barSpacing+CoMheight];
     
bars = [1:2:11; 
        2:2:12];
    
strings = [1  1   1  1  2  2  2  2  3  3  3  3  4  4  4  4  5  5  6  6  7  7  8  8;
           7  8  10 12  5  6 10 12  7  8  9 11  5  6  9 11 11 12  9 10 11 12  9 10];

% Pretension
pretension = 15; %[%]
stringRestLength = ((100-pretension)/100)*ones(24,1)*norm(nodes(1,:)-nodes(7,:));
stringRestLengthRelease = stringRestLength;
% more pretension on specific strings (Greg actuation)
overtension = 55; % 55% is the maximum
%ActuatedStrings = [21 22];
ActuatedStrings = [3 7 5 6 19 20 9 10 21 22 12 16];
stringRestLength(ActuatedStrings) = ((100-overtension)/100)*norm(nodes(1,:)-nodes(7,:));

% Physical properties

K = 10e3;              % String stiffness (N/m)
stringStiffness = K*ones(24,1);
c = 80;             % damping constant
stringDamping = c*ones(24,1);  %string damping vector
barStiffness = 100000*ones(6,1);
nodalMass = 0.42*ones(12,1);    % for a 5kg robot    
delT = 0.001; % Set the physics simulation timestep (s).

% Creation of the object 
superBall = TensegrityStructure(nodes, strings, bars, zeros(12,3), stringStiffness,...
    barStiffness, stringDamping, nodalMass, delT, delT, stringRestLength);

% Plot the object
superBallDynamicsPlot = TensegrityPlot(nodes, strings, bars, barRadius, stringRadius);
f = figure('units','normalized','outerposition',[0 0 1 1]);
generatePlot(superBallDynamicsPlot, gca);
updatePlot(superBallDynamicsPlot);
%settings to make it pretty
axis equal
view(90, 0); % X-Z view
view(3)
grid on
light('Position',[0 0 10],'Style','local')
lighting flat
colormap([0.8 0.8 1; 0 1 1])
% plot the ground
hold on
[x, y] = meshgrid(-3*barLength:0.1:3*barLength); % Generate x and y data
z = -barRadius*ones(size(x, 1)); % Generate z data
C = 2*x.*y;
ground = surf(x, y, z); % Plot the surface
ground.EdgeColor = 'none';

drawnow; % Draw and hold initial conditions
pause(0);

% Update the dynamics
displayTimespan = 1/20; % 20 fps, increase if systems can't keep up
myDynamicsUpdate(superBall, superBallDynamicsPlot, displayTimespan);
for i = 1:200 
    
    if (i==50) % release 
        % update the structure with normal pretension on every strings
        fprintf('release\n');
        minForceDensity = 1; %?
        staticTensions = getStaticTensions(superBall,minForceDensity);
        lengths = getLengths(superBall,superBall.nodePoints);
        q = staticTensions./lengths(1:length(staticTensions));
        setStringRestLengths(superBall,q,stringRestLengthRelease);
    end
    myDynamicsUpdate();
end

