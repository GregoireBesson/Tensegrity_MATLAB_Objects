%This just cleans up any extraneous timers that the function didn't clean
%up
out = timerfind;
delete(out);

%General matlab clearing
clc; close all; clear all;
addpath('../tensegrityObjects')
r = 0.10;             % Radius of top tetrahedron ring in meters
h = 0.15;             % Height of tetrahedrn in meters
rad = 0.0075;         % Radius of bars in plot in meters
N = 3;                % Number of tetrahedrons

tspan =0.05;          % time between plot updates in seconds
delT = 0.0005;         % timestep for dynamic sim in seconds
K = 2000;             %outer rim string stiffness in Newtons/meter
c = 1000;             % damping constant, too lazy to figure out units.
% set by hand waving until model looks reasonable
lims = 2*( 4/30*3.5*N/5);  % Axes limits for plotting some arbitrary
% scaling to fit axes as you increase the
% number of tetrahedrons
minQ = 100*N;          %  minimum force density in N/meter
% (amount of force per unit length in a member)

%This is used for setting the stiffness of cables in the model, here I set
%saddle cables to be 50 times stiffer than the 6 outer rim cables
stringMultiples = repmat([1 1 1 1 1 1 25 25 25 25 25 25],1,2*(N-1))';

stringStiffness = K*stringMultiples; %String stiffness vector in Newtons/meter
%needs to be ss by 1 where ss is the
%number of strings


stringDamping = c*ones(2*(N-1)*12,1);  %string damping vector
%needs to be ss by 1 where ss is the
%number of strings
nodalMass = 0.1*ones(2*N*7 - 6,1); % point mass is placed at each node with this magnitude in kg
%nodalMass([1:7 ((1:7)+N*7)]) = 0;  % set fixed nodes mass to zero

g = 9.81; %m/s^2 acelration due to gravity

%Newtons Loading of the tetrahedrons used in inverse kinematics,
%currently configured just to hold up their own mass
%Should be 3 by n where n is the number of nodes
F = [zeros(2*N*7-6,1), zeros(2*N*7-6,1), [zeros(7,1); -nodalMass(8:end)*9.81]];

%Generate reaction forces for vertical columns
F(1:7,3) = sum(F(:,3))*ones(7,1)/7;
F = zeros(size(F));
barStiffness = 100000*ones(2*N*18-12,1); %bar stiffness vector
%needs to be bb by 1 where bb is the
%number of bars

%%%  Create some global variables to share with the refreshPlot function
angle = 0;   %angle of bending about axis in x-y plane in radians
axisRot = 0; %angle that the bending axis makes with the x-axis in radians
twist = 0;   %Amount of torsional twist in radians
NR = 0.8;    %Nested Ratio for how stacked each tetrahedron is inside the next in the range 0 to 1

options = optimoptions('quadprog','Algorithm',  'interior-point-convex','Display','off');

%%%%%%%%%%%%%%%%% Create quaternions to implement angles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Note that matlab uses quaternions of the form [cos(a/2), sin(a/2)*[ux uy uz]]
% where a is the rotation and ux uy uz is a unit vector axis about which to
% rotate
quatTwist= [cos(twist) 0 0 sin(twist)];
quatBend = [cos(angle), sin(angle)*[sin(axisRot) cos(axisRot) 0]];
% Multiply quaternions is like applying successive rotations
quatt = quatmultiply(quatBend,quatTwist);

%create an N by 4 matrix of quaternons to pass to getSpineNodes
quats = [0.71 0 0.71 0;
    repmat(quatt,N-1,1)];

%%%%%%%%% Get all nodal positions for spine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spineNodes = getHexSpineNodes(r,h,NR,quats);

%%%%%% order of whch nodes make up which string and bar members %%%%%%%%%%%

% We need to create two matrices which are 2 by ss and 2 by bb where ss is
% the number of strings and bb is the number of bars. each column of this
% matrix corresponds to a string or bar and the top and bottom entries are
% the node numbers that the string or bar spans i.e. string one below goes
% between node #'s 2 and 5. By convention the first row should be the lower
% node # but the TensegrityStructure class will fix this if you mess up.

%String #'s:  1 2 3 4 5 6 7 8 9
stringList = [3 3  5  5  7  7 9 10 11 12 13 14;
              9 11 11 13 13 9 1  1  1  1  1  1];
stringNodes = zeros(2,(N-1)*9);
for i = 0:(N-2)
    %Here we are replicating the pattern of string connectivity for an N
    %number spine
    stringNodes(:,(1:12)+12*(i)) = stringList + 7*i;
end
stringNodes = [stringNodes (stringNodes+N*7)];
disp(stringNodes)
%Bar #'s:  1 2 3 4 5 6
barList = [1 1 1 1 1 1 2 3 4 5 6 7 2 4 6 3 5 7;
           2 3 4 5 6 7 3 4 5 6 7 2 4 6 2 5 7 3];
barNodes = zeros(2,(N)*18);
plotBarNodes = zeros(2,(N)*12);
for i = 0:(N-1)
    %Here we are replicating the pattern of bar connectivity for an N
    %number spine
    barNodes(:,(1:18)+18*(i)) = barList + 7*i;
    plotBarNodes(:,(1:12)+12*(i)) = barList(:,1:end-6) + 7*i;
end
barNodes = [barNodes (barNodes+N*7)];
plotBarNodes = [plotBarNodes (plotBarNodes+N*7)];
for i = 2:7
barNodes(barNodes == i+N*7) = i;
plotBarNodes(plotBarNodes == i+N*7) = i;
stringNodes(stringNodes == i+N*7) = i;
end
stringNodes(stringNodes>N*7+1) = stringNodes(stringNodes>N*7+1)-6;
barNodes(:,(7:18)+N*18) = [];
barNodes(barNodes>N*7+1) = barNodes(barNodes>N*7+1)-6;
plotBarNodes(:,(7:12)+N*12) = [];
plotBarNodes(plotBarNodes>N*7+1) = plotBarNodes(plotBarNodes>N*7+1)-6;
disp(size(barNodes))

%%%%%%%%%%%%%% Create objects for spine IK/dynamics, and plotting objects %%%%%%
% Pass all variables to the TensegrityStructure constructor to create the
% object which we call spine. This object contains methods for inverse
% kinematics as well as dynamics
spine = TensegrityStructure(spineNodes, stringNodes, barNodes, F,stringStiffness,...
    barStiffness,stringDamping,nodalMass,delT,[],[]);

%Pass variables to the constructor to create two different TensegrityPlot
%objects for plotting the IK command as well as the dynamic response of the
%structure to that command
spineCommandPlot = TensegrityPlot(spineNodes, stringNodes, plotBarNodes  , 0.005,0.001);
spineDynamicsPlot = TensegrityPlot(spineNodes, stringNodes, plotBarNodes , 0.005,0.001);

%Make sure to call this function once before initializing the dynamics
%because it will update the rest lengths
tensions = getStaticTensions(spine, minQ);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Create Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = figure('units','normalized','outerposition',[0 0 1 1]);

%%%%%%%% IK Subplot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax = subplot(1,2,1,'Parent',f,'units','normalized','outerposition',...
    [0.01 0.1 0.48 0.9]);

% use a method within TensegrityPlot class to generate a plot of the
% structure
generatePlot(spineCommandPlot,ax)
%settings to make it pretty
axis equal
view(3)
grid on
light('Position',[0 0 1],'Style','local')
%lighting flat
lighting gouraud
colormap cool% winter
xlim([-lims lims])
ylim([-lims lims])
zlim([-0.01 1.6*lims])
title('Static IK Command');


%%%%%% Dynamics Subplot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax = subplot(1,2,2,'Parent',f,'units','normalized','outerposition',...
    [0.51 0.1 0.48 0.9]);

% use a method within TensegrityPlot class to generate a plot of the
% structure
generatePlot(spineDynamicsPlot,ax)

%settings to make it pretty
axis equal
view(3)
grid on
light('Position',[0 0 1],'Style','local')
lighting flat
xlim([-lims lims])
ylim([-lims lims])
zlim([-0.01 1.6*lims])
title('Dynamics Simulation');
set(gca,'xtick',[-100:0.2:100])
set(gca,'ytick',[-100:0.2:100])
% Create an object calls from the class TensegrityCallbackFunctions which is
% just a function wrapper to hold a vector of persistent values, and some update functions
% It also holds a function for creating sliders
calls = TensegrityCallbackFunctions([angle,axisRot,twist,NR]);

%Use the function to create some slider below and attach some callback
%functions to update angle, bed axis, twist, and NR
calls.makeSlider(f,0,0,@(es,ed) updateVal(calls,es.Value,1),[-pi/12 pi/12],'angle')
calls.makeSlider(f,475,0,@(es,ed) updateVal(calls,es.Value,2),[-pi pi],'bend axis')
calls.makeSlider(f,2*475,0,@(es,ed) updateVal(calls,es.Value,3),[-pi/12 pi/12],'axial twist')
calls.makeSlider(f,3*475,0,@(es,ed) updateVal(calls,es.Value,4),[0.1 0.98],'Nested Ratio')

%A custom function needed for each tensegrity structure to update the nodes
%however you see fit, essentially the first time I call this function I set
%some persistent objects/structures (all the items after the ... on the second
%line) this speeds things up since less memory is passed to the function
% each call see the actual function for more details
pStruct = struct('minQ',minQ,'tspan',tspan,'r',r,'h',h,'N',N);
vec1 = [angle, axisRot, twist, NR];
hexSpineUpdate(vec1,...
    spine, spineCommandPlot,spineDynamicsPlot,pStruct);

%Create a function handle which only passes the vector of values we will be
%updating
spineUpdates = @(vec) hexSpineUpdate(vec);
%Create a timer to update everything, 20 fps should
%look smooth prob best not to go below this
% for i = 1:100
%    timerUpdate(calls,spineUpdates)
% end
t = timer;
t.TimerFcn = @(myTimerObj, thisEvent) timerUpdate(calls,spineUpdates);
t.Period = tspan;
t.ExecutionMode = 'fixedRate';
start(t);

