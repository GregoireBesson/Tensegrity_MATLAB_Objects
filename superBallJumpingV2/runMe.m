% Author:     ???
% ModifiedBy: Gregoire Besson, EPFL student, gregoire.besson@alumni.epfl.ch
% Date:       Winter 2018-2019
%
% Define the tensegrity structure in a first place, set up evolution
% parameters, define and run simulation for each individual

%% Simple dynamics example with 6-bar tensegrity.
clc
clear all;
close all;

addpath('../tensegrityObjects')

%% Define tensegrity structure

% Physical parameters
barLength = 0.85;               % SUPERball length, (m)
barSpacing = barLength/2;       % space between bars, usually l/2 (m)
bar_radius = 0.01;              % (m)
string_radius = 0.005;          % (m) minimum 5mm
weightEmptyEndCap = 0.100;       % (kg)
weightActEndCap =0.200;          % (kg)
pretension = 10;                % tension on strings at rest, (%)
maxTension = 60;                % max tension on actuated strings, (%)
K = 400;                        % String stiffness, (N/m)
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
% the column number correspond to the motor number in the NetworkMap
actuators = [1  5  9 13 17 19 21 23 11  3 12  4;
             2  6 10 14 18 20 22 24 15  7 16  8];

% one motor on each end cap so 12 possible motors on a 6-bars structure       
[~,nbMotorsMax] = size(actuators);

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

%plot sim Data 'NoPlot', 'PostSim' or 'RealTime' (make sim much slower!)
displayData = 'NoPlot'; 
displaySimulation = false;      % boolean to display every simulation
saveResults = true;             % boolean to save results in a mat file
initMode = 'random';            % actuators selection, 'manual' or 'random'
filename = 'output/distToGoal4over6MotI15G30c3t2.mat';
% number of end cap actuated on the robot 
nbAvlbleActuators = 6; %1+round(rand*(nbMotorsMax-1));  % should be in [1;12] 
% number of end cap actuated simultaneously should be in [1:nbAvlbleActuators]
nbUsedActuators = 4; %1+round(rand*(nbAvailableActuators-1)); 
% update the weight of the robot considering the number of actuators
nodalMass = (nbAvlbleActuators*weightActEndCap +  (12-nbAvlbleActuators)*weightEmptyEndCap)/12 *ones(12,1);
nbIndividuals = 15;              % size of population
nbGeneration = 30;               % number of generation
nbActuationCycle = 3;           % size of actuation sequence
delayAct = 0;                   % in ms

% Fitness function
fitness = 'DistToGoal';         % Jump, Dist, Jump*Dist or DistToGoal
goal = [0 -2];                  % X Y coordinates of the wanted goal

% Selection parameters
selectionMode = 'tournament';   % ranking or tournament
k = 5;                          % ranking number (ranking selection)
t = 2;                          % tournament size (tournament selection)
elitism = true;                 % copy e elites without mutation if true
e = 1;                          % number of elites to be copied 


% Mutation parameters
if (nbUsedActuators == nbAvlbleActuators)
    mutation = 'soft';
else
    mutation = 'soft';    % hard(all genes mutated), soft(1/cycle), 
end                             % ConstNbAct (change 1act with 1available
p = 0.2;                        % probability of mutation

%Create the random number stream for reproducibility:
rngParameters = RandStream('mlfg6331_64','Seed','Shuffle');
%Draw a random set of location for the actuators, same for all the indiv
AvlblActuators = datasample(rngParameters,1:12,nbAvlbleActuators,'Replace',false);
% OR user choose the end cap that are actuated (comment if you want rndm)
%AvlblActuators = [1 4 5 8 10 11]; % 1 motor on each rod

% initialization of the actuator genome, 1=actuated, 0=not
genes = zeros(nbGeneration, nbIndividuals, nbActuationCycle, 12);
traveledDist = zeros(nbIndividuals,1);
distMax = zeros(nbIndividuals,1);
zmax = zeros(nbIndividuals,1);
dist2goal = zeros(nbIndividuals,1);
performance = zeros(nbIndividuals, nbGeneration);
nbAvgActuatorsBestIndiv = zeros(1,nbGeneration);

for g = 1:nbGeneration
    
     % Status Flag
     fprintf('------------------- Generation    %d/%d -------------------\n',g,nbGeneration);
    
    for i = 1:nbIndividuals
        
        % Status Flag
        fprintf('Individual     %d/%d \n',i,nbIndividuals);
        
        % reset the actuation counter for each individual
        actuationCycleCounter = 1;
        
        % Init actuators at the first generation
        if ( g == 1 )
            % random initialization
            if (strcmpi(initMode,'random'))
                %[actuatedStrings, genes, actuationCycleCounter] = randomStrings(actuators,nbActuators,rngParameters,actuatedStrings,genes,i,g,actuationCycleCounter);
                genes(1,i,:,:) = initRandomGenome(nbActuationCycle,...
                             AvlblActuators,nbUsedActuators,rngParameters);
            % manual initialization for testing
            elseif (strcmpi(initMode,'manual'))
                % load the best individual and the genome from previous run
                load('output/jump6over12MotI100G1c1T2.mat','BestIndividual','savedGenes','fitness','oldPerf','oldPerfIndexes');
                % to just show the BestIndividual from a specific run
                %genes(1,1,:,:) = BestIndividual;
                % to continue an evolution from the last run
                %genes(1,:,:,:) = savedGenes(end,:,:,:);
                % to start evolution from the best individuals from last
                bestIndexes = oldPerfIndexes(end-nbIndividuals+1:end);
                genes(1,:,:,:) = savedGenes(end,bestIndexes,:,:);
                
                %     motors:     1  2  3  4  5  6  7  8  9 10 11 12
%                 genes(1,1,1,:) = [0  0  1  0  0  0  0  0  0  0  0  0];
%                 genes(1,1,2,:) = [0  0  0  0  0  1  0  0  0  0  0  0];
%                 genes(1,1,3,:) = [0  1  0  0  0  0  0  0  0  0  0  0];
                
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
                
            else
                error('Initialization mode unknown, enter random or manual');
            end
        end
        
        % compute the first set of strings of the Cycle of actuation
        firstStringsToActuate = genes2strings(genes,g,i,1,actuators);
        
        %% Creation of the structure
        
        % recreate structure and run simulation only on new individuals
        if ( g==1 || ...
             g>1 && (strcmpi(selectionMode,'ranking')) && ismember(i,perfIndexes(1:end-k)) || ...
             g>1 && (strcmpi(selectionMode,'tournament')) && elitism==false || ...
             g>1 && (strcmpi(selectionMode,'tournament')) && elitism && i>e ) 
        
            % Simulation and plot timesteps
            delT = 0.001;                   % timestep for dynamic sim in seconds
            delTUKF  = 0.001;               % timestep for the U Kalman Filter
            
            % creation of the object superBall
            superBall = TensegrityStructure(nodes, strings, bars, F, stringStiffness,...
                barStiffness, stringDamping, nodalMass, delT, delTUKF, ...
                stringRestLength, initMode);
            
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
                lims = barLength;
                xlim([-3*lims 3*lims])
                ylim([-3*lims 3*lims])
                zlim(1*[-0.01 lims])
                % plot the ground
                hold on
                [x, y] = meshgrid(-4*barLength:0.1:4*barLength); % Generate x and y data
                z = -bar_radius*ones(size(x, 1)); % Generate z data
                C = 2*x.*y;
                ground = surf(x, y, z); % Plot the surface
                ground.EdgeColor = 'none';
                % draw the goal as a small sphere
                if (strcmpi(fitness,'DistToGoal'))
                    r = 0.05; 
                    [X,Y,Z]= sphere;
                    surf(X*r+goal(1), Y*r+goal(2), Z*r+barLength*0.5);
                end
                %axis equal;
                
                drawnow; % Draw and hold initial conditions
                %pause(0);
                
            end
            
            %% Run dynamics
            
            displayTimespan = 1/20;     % 20fps
            % set the dynamics parameters
            myDynamicsUpdate(superBall, superBallDynamicsPlot,...
                displayTimespan, pretension, maxTension, l0,...
                actuators,i, nbActuationCycle, displaySimulation,genes,g, firstStringsToActuate);
            
            nbLoop = round(280*nbActuationCycle); %2000 -> 100sec de simulation
            
            % Simulation loop
            for l = 1:nbLoop
                myDynamicsUpdate();
                
                plotData(superBall,superBallDynamicsPlot,displayTimespan,...
                    l,nbLoop,displayData,goal);
            end
            
            %% Simulation results
            
            traveledDist(i) = superBall.TraveledDist;
            distMax(i) = superBall.DistMax;
            zmax(i) = superBall.Zmax;
            dist2goal(i) = superBall.Dist2goal;
       
        end
    end
    
    %% Evaluation
    
    % sorting mode by default
    sortingMode = 'ascend';
    
    %Select which performance to be evaluated
    if (strcmpi(fitness,'Jump'))
        performance(:,g) = zmax;
    elseif (strcmpi(fitness,'Dist'))
        performance(:,g) = traveledDist;
    elseif (strcmpi(fitness,'DistMax'))
        performance(:,g) = distMax;
    elseif (strcmpi(fitness,'Dist*Jump'))
        performance(:,g) = traveledDist.*zmax;
    elseif (strcmpi(fitness,'DistToGoal'))
        performance(:,g) = dist2goal;
        % change the sorting mode, cause we want to minimize this fitness
        sortingMode = 'descend';
    else 
        error('Unknown fitness performance');
    end
    
    % copy the perf of the unchanged indiv from last generation
    if (g > 1)
        if (strcmpi(selectionMode,'ranking'))
            performance(perfIndexes(k+1:end),g) = performance(perfIndexes(k+1:end),g-1);
        elseif ((strcmpi(selectionMode,'tournament')) && elitism)           % TODO correct BUG if e>2
            performance(1:e,g) = sortedPerformance(end);
        end
    end
    
    %% Selection
      
    % sort the perfomances depending if we want to maximize or minimize Fit
    [sortedPerformance, perfIndexes] = sort(performance(:,g), sortingMode);
    
    % compute the average number of actuator per cycle used by BestInd
    if strcmpi(initMode,'random')
        nbAvgActuatorsBestIndiv(g) = length( find(genes(g,perfIndexes(end),:,:)==1) )/nbActuationCycle;
    end
    
    if ( (g < nbGeneration) && (strcmpi(selectionMode,'tournament')) )
        % fill the next generation with tournament method
        genes = tournamentSelection(genes,nbIndividuals,performance,g,t,elitism,e,sortingMode);
    end

    %% Mutation
    
    % mutate the individuals that need to be mutated depending on the
    % selection method
    if (g < nbGeneration)
        genes = genesMutation(mutation,genes,g,nbIndividuals,perfIndexes,AvlblActuators,nbActuationCycle,k,p,selectionMode,elitism,e);
    end
    
    
    % save partial results after each generation, overwrite 
    if (saveResults)
        if (strcmpi(sortingMode,'ascend'))
            [~,BestIndividualIndex] = max(performance(:,g));
        else
            [~,BestIndividualIndex] = min(performance(:,g));
        end
        BestIndividual = genes(nbGeneration,BestIndividualIndex,:,:);
        savedGenes = genes;
        oldPerf = performance;
        oldPerfIndexes = perfIndexes;
        save(filename,'BestIndividual','savedGenes','fitness','oldPerf','oldPerfIndexes')
        fprintf('Results saved \n');
    end
end

%% Save and plot 

fprintf('Evolution complete !\n');

if (saveResults) 
    plotEvolution(performance,sortedPerformance,fitness,nbAvlbleActuators,nbUsedActuators,nbGeneration,nbIndividuals,nbActuationCycle,t,p,nbAvgActuatorsBestIndiv);
end