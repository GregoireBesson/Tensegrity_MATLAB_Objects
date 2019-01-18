% Author:   Gregoire Besson, EPFL student, gregoire.besson@alumni.epfl.ch
% Date:     Winter 2018-2019
%
% Subplot structure data such as String tensions, forces,
% lengths and rest lengths in real time
%
% Inputs: - obj (TensegrityStructure)
%         - plotObj (TensegrityPlot)
%         - tspan (display timespan to derive the time)
%         - i   (silulation loop index)
%         - nbLoop (total nb of simulation loops)
%         - plotCmd (Plot command: 'RealTime' (super slow) or 'PostSim')
%
% Output: - suplot of the data of interest
%         - TraveledDist (Euclidean Dist of CoM from start to end
%         - Zmax (maximum jump height achieved)
%         - DistMax (maximum distance from start reached along the sim)

function plotData(obj,plotObj,tspan,i,nbLoop,plotCmd)

persistent X Y Z f2 timeVector membersLengthDataStore Zmax DistMax... 
    stringRestLengthDataStore stringTensionsDataStore stringToPlot CoM goal

% initialisation
if (i==1)
    % constants definition
    X = 1;
    Y = 2;
    Z = 3;
    
    if (strcmpi(plotCmd,'RealTime'))
        f2 = figure('Name','Data','NumberTitle','off');
    end
    timeVector = (tspan:tspan:nbLoop*tspan);
    stringRestLengthDataStore = zeros(24,nbLoop);
    membersLengthDataStore = zeros(30,nbLoop);
    stringTensionsDataStore = zeros(24,nbLoop);
    CoM = zeros(nbLoop,3);
    stringToPlot = [9 10 11 15 23 24];
    Zmax = 0;
    DistMax = 0;
end

% store Data
stringRestLengthDataStore(:,i) = obj.stringRestLength;
membersLengthDataStore(:,i) = obj.memberLength;
stringTensionsDataStore(:,i) = obj.stringTensions;
CoM(i,:) = mean(plotObj.nodePoints);
% check max height
if (CoM(i,Z) > Zmax)
   Zmax = CoM(i,Z); 
end
% no sqrt to save computation time
currentDist = (CoM(i,X)-CoM(1,X))^2 + (CoM(i,Y)-CoM(1,Y))^2;
if ( currentDist > DistMax )
   DistMax =  currentDist;
end

% Real Time Plot (time consuming)
if (strcmpi(plotCmd,'RealTime'))
    figure(f2);
    subplot(3,1,1);
    plot(timeVector(1:i),stringRestLengthDataStore(stringToPlot,1:i))
    title('Strings Rest Lengths')
    xlabel('Time [s]')
    ylabel('[m]')
    grid on
    subplot(3,1,2);
    plot(timeVector(1:i),membersLengthDataStore(stringToPlot,1:i))
    title('Members Lengths')
    xlabel('Time [s]')
    ylabel('[m]')
    grid on
    subplot(3,1,3);
    plot(timeVector(1:i),stringTensionsDataStore(stringToPlot,1:i))
    title('String Tensions')
    xlabel('Time [s]')
    ylabel('[Newtons]')
    grid on
    hold on
    
% Postsimulation Plot 
elseif strcmpi(plotCmd,'PostSim') && (i==nbLoop)
    obj.TraveledDist = sqrt((CoM(end,X)-CoM(1,X))^2 + (CoM(end,Y)-CoM(1,Y))^2);
    obj.Zmax = Zmax;
    obj.DistMax = DistMax;
    
    f2 = figure('Name','Data','NumberTitle','off');
    % Set Size and position
    set(f2, 'Position', [20 20 1000 600]);
    figure(f2);
    subplot(3,2,1);
    %plot(timeVector,stringRestLengthDataStore(stringToPlot,:))
    plot(timeVector,stringRestLengthDataStore)
    title('Strings Rest Lengths')
    xlabel('Time [s]')
    ylabel('[m]')
    grid on
    subplot(3,2,3);
    %plot(timeVector,membersLengthDataStore(stringToPlot,:))
    plot(timeVector,membersLengthDataStore)
    title('Members Lengths')
    xlabel('Time [s]')
    ylabel('[m]')
    grid on
    subplot(3,2,5);
    %plot(timeVector,stringTensionsDataStore(stringToPlot,:))
    plot(timeVector,stringTensionsDataStore)
    title('String Tensions')
    xlabel('Time [s]')
    ylabel('[Newtons]')
    grid on
    hold on
    % Trajectories
    subplot(3,2,2);
    plot3(CoM(:,X),CoM(:,Y),CoM(:,Z))
    title('3D CoM Trajectory')
    xlabel('X [m]')
    ylabel('Y [m]')
    zlabel('Z [m]')
    grid on
    hold on
    plot3(CoM(1,X),CoM(1,Y),CoM(1,Z),'go','DisplayName','Start')
    plot3(CoM(end,X),CoM(end,Y),CoM(end,Z),'ro','DisplayName','End')
    lims = 1;
    xlim([-1.2*lims 1.2*lims])
    ylim([-1.2*lims 1.2*lims])
    zlim(1*[-0.01 lims])
    subplot(3,2,4);
    plot(CoM(:,X),CoM(:,Y),'DisplayName','CoM Trace')
    title('XY CoM Trajectory')
    xlabel('X [m]')
    ylabel('Y [m]')
    legend show
    grid on
    hold on
    plot(CoM(1,X),CoM(1,Y),'go','DisplayName','Start')
    plot(CoM(end,X),CoM(end,Y),'ro','DisplayName','End')
    xlim([-1.2*lims 1.2*lims])
    ylim([-1.2*lims 1.2*lims])
    subplot(3,2,6);
    plot(timeVector,CoM(:,Z))
    title('CoM height')
    xlabel('Time [s]')
    zlabel('Z [m]')
    grid on
    
    %No plot, just store Dist and Zmax
    elseif strcmpi(plotCmd,'NoPlot') && (i==nbLoop)
    obj.TraveledDist = sqrt((CoM(end,X)-CoM(1,X))^2 + (CoM(end,Y)-CoM(1,Y))^2);
    obj.Zmax = Zmax;
    obj.DistMax = DistMax;
end