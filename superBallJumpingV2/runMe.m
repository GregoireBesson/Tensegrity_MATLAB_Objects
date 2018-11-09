%% Simple dynamics example with 6-bar tensegrity.
clc
clear all;
close all;

addpath('../tensegrityObjects')

%global restLengthDataStore;
%global membersLengthDataStore;

%% Define tensegrity structure

% Physical parameters
barLength = 1;                  % SUPERball length, (m)
barSpacing = barLength/2;       % space between bars, usually l/2 (m)
bar_radius = 0.010;             % (m)
string_radius = 0.005;          % (m) minimum 5mm
nodalMass = 0.42*ones(12,1);    % target: a 5kg robot
pretension = 15;                % tension on strings at rest, (%)
maxTension = 50;                % max tension on actuated strings, (%)
K = 1000;                       % String stiffness, (N/m)
c = 80;                         % viscous friction coef, (Ns/m)
stringStiffness = K*ones(24,1); % String stiffness (N/m)
barStiffness = 100000*ones(6,1);% Bar stiffness (N/m)
stringDamping = c*ones(24,1);   % string damping vector
F = zeros(12, 3);               % n by 3 matrix nodal forces

% nodes location
nodes = [-barSpacing*0.5     barLength*0.5  0;
         -barSpacing*0.5    -barLength*0.5  0;
          barSpacing*0.5     barLength*0.5  0;
          barSpacing*0.5    -barLength*0.5  0;
          0             -barSpacing*0.5     barLength*0.5;
          0             -barSpacing*0.5    -barLength*0.5;
          0              barSpacing*0.5     barLength*0.5;
          0              barSpacing*0.5    -barLength*0.5;        
          barLength*0.5  0             -barSpacing*0.5;
         -barLength*0.5  0             -barSpacing*0.5;
          barLength*0.5  0              barSpacing*0.5;
         -barLength*0.5  0              barSpacing*0.5];
     
% bar connectivity
bars = [1:2:11; 
        2:2:12];
    
% string connectivity
strings = [1  1   1  1  2  2  2  2  3  3  3  3  4  4  4  4  5  5  6  6  7  7  8  8;
           7  8  10 12  5  6 10 12  7  8  9 11  5  6  9 11 11 12  9 10 11 12  9 10];
       
% vector containing which strings are going to be pulled
%actuatedStrings = [3 7 5 6 19 20 9 10 21 22 12 16]; %to show an upright mvt
 actuatedStrings = [9 10 11 15 23 24 17 18 4 8 5 6]; %to start in good pos
%actuatedStrings = [9 10];                  %1 lower actuation
%actuatedStrings = [9 10 11 15];            %2 lower actuations
%actuatedStrings = [9 10 11 15 23 24];      %3 lower actuations

% Compute rest lengths from a certain pretension
l0 = norm(nodes(1,:)-nodes(7,:));       % initial string length
stringRestLength = ((100-pretension)/100)*ones(24,1)*l0;
     
% rotate the structure so that it land on a triangle base  
Mz  = makehgtform('zrotate', -45*pi/180);   %rot Matrix 45deg around z
My  = makehgtform('yrotate', 55*pi/180);    %rot Matrix 55deg around y
nodes = (Mz(1:3,1:3)*nodes')';
nodes = (My(1:3,1:3)*nodes')';

% set the droping height
CoMz = 0.5;                                 % (m)
nodes(:,3) = nodes(:,3) + CoMz;             % shift all the nodes in z

%% Creation of the structure

% Simulation and plot timesteps
delT = 0.001;                   % timestep for dynamic sim in seconds
delTUKF  = 0.001;               % timestep for the U Kalman Filter

% creation of the object superBall
superBall = TensegrityStructure(nodes, strings, bars, F, stringStiffness,...
    barStiffness, stringDamping, nodalMass, delT, delTUKF, stringRestLength);

%% Create dynamics display

superBallDynamicsPlot = TensegrityPlot(nodes, strings, bars, ...
    bar_radius, string_radius);
f = figure('units','normalized','outerposition',[0 0 1 1]);
% use a method within TensegrityPlot class to generate a plot of the
% structure
generatePlot(superBallDynamicsPlot,gca);
updatePlot(superBallDynamicsPlot);

%settings to make it pretty
axis equal
view(90, 0); % X-Z view
view(3)
grid on
light('Position',[0 0 10],'Style','local')
lighting flat
colormap([0.8 0.8 1; 0 1 1]);
lims = 1.2*barLength;
%xlim([-1.2*lims 1.2*lims])
%ylim([-1.2*lims 1.2*lims])
zlim(1.6*[-0.01 lims])
% plot the ground
hold on
[x, y] = meshgrid(-3*barLength:0.1:3*barLength); % Generate x and y data
z = -bar_radius*ones(size(x, 1)); % Generate z data
C = 2*x.*y;
ground = surf(x, y, z); % Plot the surface
ground.EdgeColor = 'none';

drawnow; % Draw and hold initial conditions
pause(0);

%% Run dynamics

displayTimespan = 1/20;     % 20fps
% set the dynamics parameters
myDynamicsUpdate(superBall, superBallDynamicsPlot, displayTimespan, ...
    actuatedStrings, pretension, maxTension, l0);

nbLoop = 300;

% Simulation loop
for i = 1:nbLoop
    myDynamicsUpdate();
    
    %real time plot of data
    plotData(superBall,i,nbLoop,'PostSim');
end

% figure(2)
% plot(stringTensionsDataStore)
% figure(3)
% plot(membersLengthDataStore)
% figure(4)
% plot(stringRestLengthDataStore)