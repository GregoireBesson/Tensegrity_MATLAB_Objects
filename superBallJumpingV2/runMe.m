% Author:     Jeff Friesen
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
barLength = 0.50;               % SUPERball length, (m)
passiveCableRestLength = 0.31;  % (m)
activeCableRestLength = 0.34;   % (m)
barSpacing = barLength/2;       % space between bars, usually l/2 (m)
bar_radius = 0.020;             % (m)
string_radius = 0.005;          % (m) minimum 5mm
weightEmptyEndCap = 0.25;      % (kg)
weightActEndCap = 0.25;        % (kg)
pretension = 5;                 % tension on strings at rest, (%)
maxTension = 55;                % max tension on actuated strings, (%)
%K = 700;                       % String stiffness, (N/m)
Kactive = 1000;  %800 for loc    % Stiffness of active cables, (N/m)
Kpassive = 1000; %800 for loc    % Stiffness of passive cables, (N/m)
c = 80;                         % viscous friction coef, (Ns/m)
stringStiffness = Kpassive*ones(24,1); % String stiffness (N/m)
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
strings = [1 1  1  1 2 2  2  2 3 3 3  3 4 4  4  4  5  5  6  6  7  7  8  8;
           7 8 10 12 5 6 10 12 7 8 9 11 5 6  9 11 11 12  9 10 11 12  9 10];

% each column contain the pair of colum numbers in strings to actuate
% the column number correspond to the motor number in the NetworkMap
actuators = [1  5  9 13 17 19 21 23 11  3 12  4;
             2  6 10 14 18 20 22 24 15  7 16  8];

% new actuators pull on 1 string, table contain the string         
newActuators = [6 17 4 24 15 9];

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
displayData = 'PostSim'; 
displaySimulation = true;      % boolean to display every simulation
saveResults = false;             % boolean to save results in a mat file
initMode = 'manual';            % actuators selection, 'manual' or 'random'
oldFilename = 'output/SMALLERdist9over12MotI25G25k5.mat';
newFilename = 'output/SMALLERdist9over12MotI25G25k5.mat';
% draw rndm available actuators for every indiv if true
changeAvlActuators = false;
% number of end cap actuated on the robot 
nbAvlbleActuators = 6; %1+round(rand*(nbMotorsMax-1));%should be in [1;12] 
%number of endcap actuated simultaneously should be in[1:nbAvlbleActuators]
nbUsedActuators = 6; %1+round(rand*(nbAvailableActuators-1)); 
% update the weight of the robot considering the number of actuators
nodalMass = (nbAvlbleActuators*weightActEndCap +  ...
    (12-nbAvlbleActuators)*weightEmptyEndCap)/12 *ones(12,1);
nbIndividuals = 1;           % size of population
nbGeneration = 1;               % number of generation
nbActuationCycle = 1;           % size of actuation sequence
delayAct = 0;                   % in ms

% Fitness function
fitness = 'jump';               % Jump, Dist, Jump*Dist or DistToGoal
goal = [0 -2];                  % X Y coordinates of the wanted goal

% Selection parameters
selectionMode = 'tournament';   % ranking or tournament
k = 5;                          % ranking number (ranking selection)
t = 2;                          % tournament size (tournament selection)
elitism = true;                 % copy e elites without mutation if true
e = 1;                          % number of elites to be copied 

% Mutation parameters
mutation = 'soft';              % hard(all genes mutated), soft(1/cycle), 
                                % ConstNbAct (change 1act with 1available
p = 0.2;                        % probability of mutation

%Create the random number stream for reproducibility:
rngParameters = RandStream('mlfg6331_64','Seed','Shuffle');

%Draw a random set of location for the actuators, same for all the indiv
AvlblActuators = datasample(rngParameters,1:12,nbAvlbleActuators,...
    'Replace',false);
% OR user choose the end cap that are actuated (comment if you want rndm)
%AvlblActuators = [1 4 5 8 9 12]; % 1 motor on each rod

% initialization of the actuator genome, 1=actuated, 0=not
genes = zeros(nbGeneration, nbIndividuals, nbActuationCycle, 6);
traveledDist = zeros(nbIndividuals,1);
distMax = zeros(nbIndividuals,1);
zmax = zeros(nbIndividuals,1);
dist2goal = zeros(nbIndividuals,1);
performance = zeros(nbIndividuals, nbGeneration);
nbAvgActuatorsBestIndiv = zeros(1,nbGeneration);
AvlActuatorsPerInd = zeros(nbGeneration,nbIndividuals,nbAvlbleActuators);

for g = 1:nbGeneration
    
     % Status Flag
     fprintf('------------------ Generation %d/%d ------------------\n',...
         g,nbGeneration);
    
    for i = 1:nbIndividuals
        
        % Status Flag
        fprintf('Individual     %d/%d \n',i,nbIndividuals);
        
        % reset the actuation counter for each individual
        actuationCycleCounter = 1;
        
        % Init actuators at the first generation
        if ( g == 1 )
            if (changeAvlActuators)
                %Draw a random set of location for the actuators, 
                %same for all the indiv
                AvlblActuators = datasample(rngParameters,1:12,...
                                 nbAvlbleActuators,'Replace',false);
                AvlActuatorsPerInd(g,i,:) = AvlblActuators;
            end
            % random initialization
            if (strcmpi(initMode,'random'))
                genes(1,i,:,:) = initRandomGenome(nbActuationCycle,...
                             AvlblActuators,nbUsedActuators,rngParameters);
                AvlActuatorsPerInd(g,i,:) = AvlblActuators;
            % manual initialization for testing
            elseif (strcmpi(initMode,'manual'))
                if (i==1)
                    % load the best individual and the genome from previous run
                    %load(oldFilename,'BestIndividual','savedGenes',...
                    %    'fitness','oldPerf','oldPerfIndexes','oldAvlbAct');
                    % to continue an evolution from the last run
                    %genes(1,:,:,:) = savedGenes(end,:,:,:);
                    % to start evolution from the best individuals from last
                    %bestIndexes = oldPerfIndexes(end-nbIndividuals+1:end);
                    %genes(1,:,:,:) = savedGenes(end,bestIndexes,:,:);
                    %AvlActuatorsPerInd(1,:,:) = oldAvlbAct(end,bestIndexes,:);
                    %     motors:     1  2  3  4  5  6  
                    genes(1,1,1,:) = [1  1  1  1  1  1];% compressed vertexes
                   
                    
                    %genes(1,1,1,:) = [0  0  0  0  0  0  0  0  0  0  0  0];
                    %genes(1,1,1,:) = [0  0  1  0  0  0  0  1  1  0  0  0];      % intuitive directional jump
                    %genes(1,1,2,:) = [0  1  0  0  1  0  0  0  0  0  1  0];
                    %genes(1,1,2,:) = [1  0  0  1  1  0  0  1  1  0  0  1];
                    %genes(1,1,3,:) = [1  0  0  1  1  0  0  1  1  0  0  1];
                end
                
                %     motors:     1  2  3  4  5  6  7  8  9 10 11 12
%                 genes(1,1,1,:) = [0  0  1  0  0  0  0  0  0  0  0  0];
%                 genes(1,1,2,:) = [0  0  0  0  0  1  0  0  0  0  0  0];
%                 genes(1,1,3,:) = [0  1  0  0  0  0  0  0  0  0  0  0];
                
                %actuatedStrings(i,1,:) = [13 14 9 10 3 7 19 20 23 24 1 2];
                %actuatedStrings(i,2,:) = [3 7 5 6 19 20 9 10 21 22 12 16];
                %to show an upright mvt:
                %actuatedStrings = [3 7 5 6 19 20 9 10 21 22 12 16]; 
                %to start in good pos:
                %actuatedStrings = [9 10 11 15 23 24 17 18 4 8 5 6]; 
                %actuatedStrings = [9 10];              %1 lower actuation
                %actuatedStrings = [9 10 11 15];        %2 lower actuations
                %actuatedStrings = [9 10 11 15 23 24];  %3 lower actuations
                %actuatedStrings = [9 10 19 20];  %1 lower and 1 upper act
                %actuatedStrings = [17 18];       %1 actuation to test stif
                %actuatedStrings = [ 9 10;
                %                   23 24;
                %                    3  7]; % serial actuation (row by row)
                
            else
                error('Initialization mode unknown, random or manual');
            end
        end
        
        %activeStrings = reshape(newActuators(:,genes(g,i,1,:)==1),...
        %                [1,length(newActuators(:,genes(g,i,1,:)==1))]);
        %TEST TO ACTUATE 1 string per actuator (vertical string)
        %activeStrings = [6 17 4 15 9 24];
        stringStiffness(newActuators) = Kactive;
       
        
        % compute the first set of strings of the Cycle of actuation
        firstStringsToActuate = genes2strings(genes,g,i,1,newActuators);
        %firstStringsToActuate = [6 17 4 15 9 24];   
        % Needed if active cables have differente length than passive ones:
        stringRestLength(newActuators) = activeCableRestLength;
        
        %% Creation of the structure
        
        % recreate structure and run simulation only on new individuals
        if ( g==1 || ...
             g>1 && (strcmpi(selectionMode,'ranking')) && ...
             ismember(i,perfIndexes(1:end-k)) || ...
             g>1 && (strcmpi(selectionMode,'tournament')) && ...
             elitism==false || ...
             g>1 && (strcmpi(selectionMode,'tournament')) && ...
             elitism && i>e ) 
        
            % Simulation and plot timesteps
            delT = 0.001;           % timestep for dynamic sim in seconds
            delTUKF  = 0.001;       % timestep for the U Kalman Filter
            
            % creation of the object superBall
            superBall = TensegrityStructure(nodes, strings, bars, F, ...
                stringStiffness,barStiffness, stringDamping, nodalMass, ...
                delT, delTUKF, stringRestLength, initMode);
            
            %% Create dynamics display
            
            superBallDynamicsPlot = TensegrityPlot(nodes, strings, bars, ...
                bar_radius, string_radius);
            
            if (displaySimulation)
                
                f = figure('units','normalized','outerposition',[0 0 1 1]);
                % use a method within TensegrityPlot class to generate 
                % a plot of the structure
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
                lims = barLength*2;
                xlim([-2*lims 2*lims])
                ylim([-2*lims 2*lims])
                zlim(1*[-0.01 1*lims])
                % plot the ground
                hold on
                [x, y] = meshgrid(-4*barLength:0.1:4*barLength); 
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
                pause(0);
                
            end
            
            %% Run dynamics
            
            displayTimespan = 1/20;     % 20fps
            % set the dynamics parameters
            myDynamicsUpdate(superBall, superBallDynamicsPlot,...
                displayTimespan, pretension, maxTension, activeCableRestLength,...
                newActuators,i, nbActuationCycle, displaySimulation,genes,...
                g,firstStringsToActuate);
            
            % will define the simulation duration
            nbLoop = round(280*nbActuationCycle);
        
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
        jump = zmax - barLength/2.5 %(zmax of CoM - initial CoM z)
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
            performance(perfIndexes(k+1:end),g) = ...
                performance(perfIndexes(k+1:end),g-1);
        elseif ((strcmpi(selectionMode,'tournament')) && elitism)           % TODO correct BUG if e>2
            performance(1:e,g) = sortedPerformance(end);
        end
    end
    
    %% Selection
      
    % sort the perfomances depending if we want to maximize or minimize Fit
    [sortedPerformance, perfIndexes] = sort(performance(:,g), sortingMode);
    
    % compute the average number of actuator per cycle used by BestInd
    nbAvgActuatorsBestIndiv(g) = ...
          length(find(genes(g,perfIndexes(end),:,:)==1) )/nbActuationCycle;
    
    if ( (g < nbGeneration) && (strcmpi(selectionMode,'tournament')) )
        % fill the next generation with tournament method
        [genes,AvlActuatorsPerInd] = tournamentSelection(genes,...
            nbIndividuals,performance,g,t,elitism,e,sortingMode,...
            AvlActuatorsPerInd);
    end

    %% Mutation
    
    % mutate the individuals that need to be mutated depending on the
    % selection method
    if (g < nbGeneration)
        genes = genesMutation(mutation,genes,g,nbIndividuals,...
            perfIndexes,AvlActuatorsPerInd,nbActuationCycle,k,p,...
            selectionMode,elitism,e);
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
        oldAvlbAct = AvlActuatorsPerInd;
        save(newFilename,'BestIndividual','savedGenes','fitness',...
            'oldPerf','oldPerfIndexes','oldAvlbAct')
        fprintf('                                     Results saved \n\n');
    end
end

%% Save and plot 

fprintf('Evolution complete !\n');

if (saveResults) 
    plotEvolution(performance,sortedPerformance,fitness,...
        nbAvlbleActuators,nbUsedActuators,nbGeneration,nbIndividuals,...
        nbActuationCycle,k,p,nbAvgActuatorsBestIndiv);
end