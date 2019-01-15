%% Simple dynamics example with 6-bar tensegrity.
clc
clear all;
close all;

addpath('../tensegrityObjects')

%% Define tensegrity structure

% Physical parameters
barLength = 0.5;                % SUPERball length, (m)
barSpacing = barLength/2;       % space between bars, usually l/2 (m)
bar_radius = 0.01;              % (m)
string_radius = 0.005;          % (m) minimum 5mm
nodalMass = 0.2*ones(12,1);     % target: a 5kg robot                      % Mockup is 120g (10g per node) but equations cannot be solved, minimum mass is 3x larger so I divided g by 3
pretension = 15;                % tension on strings at rest, (%)
maxTension = 60;                % max tension on actuated strings, (%)
K = 500;                        % String stiffness, (N/m)
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

% each column contain the pair of colum numbers in strings to actuate
actuators = [1  5  9 13 17 19 21 23 11  3 12  4;
             2  6 10 14 18 20 22 24 15  7 16  8];

% Compute rest lengths from a certain pretension
l0 = norm(nodes(1,:)-nodes(7,:));       % initial string length
stringRestLength = ((100-pretension)/100)*ones(24,1)*l0;
     
% rotate the structure so that it land on a triangle base (nodes 3,8 and 9)
Mz  = makehgtform('zrotate', -45*pi/180);   %rot Matrix 45deg around z
My  = makehgtform('yrotate', 55*pi/180);    %rot Matrix 55deg around y
nodes = (Mz(1:3,1:3)*nodes')';
nodes = (My(1:3,1:3)*nodes')';

% set the droping height
CoMz = barLength/2;                         % (m)
nodes(:,3) = nodes(:,3) + CoMz;             % shift all the nodes in z

%% Evolution parameters

displaySimulation = 0;          % boolean to display every simulation
selectionMode = 'random';       % selectionMode can be 'manual' or 'random'
nbActuators = 6;                % should be in [1;12]
delayAct = 0;                   % in ms
nbIndividuals = 5;              % size of population
nbGeneration = 5;               % number of generation
nbActuations = 2;               % size of actuation sequence
k = 2;                          % selection parameter (remove k worst ind.)
fitness = 'distMax';               % perf to be evaluated (jump, dist or both)              

traveledDist = zeros(nbIndividuals,1);
distMax = zeros(nbIndividuals,1);
zmax = zeros(nbIndividuals,1);
actuatedStrings = zeros(nbIndividuals, nbActuations, 2*nbActuators);
performance = zeros(nbIndividuals, nbGeneration);

%Create the random number stream for reproducibility:
rngParameters = RandStream('mlfg6331_64','Seed','Shuffle');

for g = 1:nbGeneration
    for i = 1:nbIndividuals
        actuationCounter = 0;
        
        % randomm actuators for evolution
        if ( g == 1 || strcmpi(selectionMode,'random'))
            [actuatedStrings,actuationCounter] = randomStrings(actuators,nbActuators,rngParameters,actuatedStrings,i,actuationCounter);
        end
        
        % manual actuators for testing
        if (strcmpi(selectionMode,'manual'))
            %actuatedStrings(i,1,:) = [13 14 9 10 3 7 19 20 23 24 1 2];
            %actuatedStrings(i,2,:) = [3 7 5 6 19 20 9 10 21 22 12 16];
            
            %actuatedStrings = [3 7 5 6 19 20 9 10 21 22 12 16]; %to show an upright mvt
            %actuatedStrings = [9 10 11 15 23 24 17 18 4 8 5 6]; %to start in good pos
            %actuatedStrings = [9 10];                  %1 lower actuation
            %actuatedStrings = [9 10 11 15];            %2 lower actuations
            %actuatedStrings = [9 10 11 15 23 24];      %3 lower actuations
            %actuatedStrings = [9 10 19 20];            %1 lower and 1 upper actuation
            %actuatedStrings = [17 18];                 %1 actuation to test stiffness
            %actuatedStrings = [ 9 10;
            %                   23 24;
            %                    3  7];                 % serial actuation (row by row)
        end
        
        %% Creation of the structure
        
        % Simulation and plot timesteps
        delT = 0.001;                   % timestep for dynamic sim in seconds
        delTUKF  = 0.001;               % timestep for the U Kalman Filter
        
        % creation of the object superBall
        superBall = TensegrityStructure(nodes, strings, bars, F, stringStiffness,...
            barStiffness, stringDamping, nodalMass, delT, delTUKF, ...
            stringRestLength, selectionMode);
        
        %% Create dynamics display
        
        superBallDynamicsPlot = TensegrityPlot(nodes, strings, bars, ...
            bar_radius, string_radius);
        
        if (displaySimulation)
        
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
            lims = 2*barLength;
            xlim([-1.2*lims 1.2*lims])
            ylim([-1.2*lims 1.2*lims])
            zlim(1*[-0.01 lims])
            % plot the ground
            hold on
            [x, y] = meshgrid(-3*barLength:0.1:3*barLength); % Generate x and y data
            z = -bar_radius*ones(size(x, 1)); % Generate z data
            C = 2*x.*y;
            ground = surf(x, y, z); % Plot the surface
            ground.EdgeColor = 'none';
            
            drawnow; % Draw and hold initial conditions
            %pause(0);
        
        end
        
        %% Run dynamics
        
        displayTimespan = 1/20;     % 20fps
        % set the dynamics parameters
        [~,~] = myDynamicsUpdate(superBall, superBallDynamicsPlot,...
            displayTimespan, actuatedStrings, pretension, maxTension, l0,...
            actuators, nbActuators,rngParameters,i,actuationCounter,...
            nbActuations, displaySimulation);
        
        nbLoop = round((100/100)*180*nbActuations); %2000 -> 100sec de simulation
        
        % Simulation loop
        for l = 1:nbLoop
            [actuatedStrings,actuationCounter] = myDynamicsUpdate();
            
            %display Data 'NoPlot', 'PostSim' or 'RealTime' (make simulation much slower!)
            plotData(superBall,superBallDynamicsPlot,displayTimespan,...
                l,nbLoop,'NoPlot');
        end
        
        %% Simulation results
        
        traveledDist(i) = superBall.TraveledDist;
        distMax(i) = superBall.DistMax;
        zmax(i) = superBall.Zmax;
        
    end
    
    %% Select which performance to be evaluated
    
    if (strcmpi(fitness,'jump'))
        performance(:,g) = zmax;
    elseif (strcmpi(fitness,'dist'))
        performance(:,g) = traveledDist;
    elseif (strcmpi(fitness,'distMax'))
        performance(:,g) = distMax;
    else
        performance(:,g) = traveledDist.*zmax;
    end
    
    % sort the perfomances
    [sortedPerformance, Indexes] = sort(performance(:,g),'ascend');
    
    %% Mutation
    % send the previous sequences to the next
    pastActuatedStrings = actuatedStrings(:,1:nbActuations,:);
    nextActuatedStrings = actuatedStrings(:,1:nbActuations,:);
    % replace the k worst indiv with mutated version of the k best
    % OR with random strings
    for i = 1:k
        for a = 0:(nbActuations-1)
            [nextActuatedStrings,~] = randomStrings(actuators,nbActuators,...
                rngParameters,nextActuatedStrings,Indexes(i),a);
        end
    end
    
    %update the strings to actuate for next generation
    actuatedStrings = nextActuatedStrings;
    
    % force the next selection mode to manual cause we have now defined
    % sequence
    selectionMode = 'manual';

end

plotEvolution(performance,fitness,nbActuators,nbGeneration,nbIndividuals,nbActuations,k);