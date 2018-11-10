% Function that subplot structure data such as String tensions, forces,
% lengths and rest lengths in real time
%
% Inputs: - obj (TensegrityStructure)
%         - i   (silulation loop index)
% Output: - suplot of the data of interest

function plotData(obj,i,nbLoop,plotCmd)

persistent f2 timeVector membersLengthDataStore ... 
    stringRestLengthDataStore stringTensionsDataStore

% initialisation
if (i==1)
    if (strcmpi(plotCmd,'RealTime'))
        f2 = figure('Name','Data','NumberTitle','off');
    end
    timeVector = (obj.delT:obj.delT:nbLoop*obj.delT);
    stringRestLengthDataStore = zeros(24,nbLoop);
    membersLengthDataStore = zeros(30,nbLoop);
    stringTensionsDataStore = zeros(24,nbLoop);
end

% store Data
stringRestLengthDataStore(:,i) = obj.stringRestLength;
membersLengthDataStore(:,i) = obj.memberLength;
stringTensionsDataStore(:,i) = obj.stringTensions;

% Real Time Plot (time consuming)
if (strcmpi(plotCmd,'RealTime'))
    figure(f2);
    subplot(3,1,1);
    plot(timeVector(1:i),stringRestLengthDataStore(:,1:i))
    title('Strings Rest Lengths')
    xlabel('Time [s]')
    ylabel('[m]')
    grid on
    subplot(3,1,2);
    plot(timeVector(1:i),membersLengthDataStore(:,1:i))
    title('Members Lengths')
    xlabel('Time [s]')
    ylabel('[m]')
    grid on
    subplot(3,1,3);
    plot(timeVector(1:i),stringTensionsDataStore(:,1:i))
    title('String Tension')
    xlabel('Time [s]')
    ylabel('[Newtons]')
    grid on
    hold on
% Postsimulation Plot 
elseif (i==nbLoop)
    f2 = figure('Name','Data','NumberTitle','off');
    figure(f2);
    subplot(3,1,1);
    plot(timeVector,stringRestLengthDataStore)
    title('Strings Rest Lengths')
    xlabel('Time [s]')
    ylabel('[m]')
    grid on
    subplot(3,1,2);
    plot(timeVector,membersLengthDataStore)
    title('Members Lengths')
    xlabel('Time [s]')
    ylabel('[m]')
    grid on
    subplot(3,1,3);
    plot(timeVector,stringTensionsDataStore)
    title('String Tension')
    xlabel('Time [s]')
    ylabel('[Newtons]')
    grid on
    hold on
end