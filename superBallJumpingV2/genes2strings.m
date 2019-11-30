% Author:   Gregoire Besson, EPFL student, gregoire.besson@alumni.epfl.ch
% Date:     Winter 2018-2019
%
% Convert the genome into command to actuate corresponding strings

function stringsToActuate = genes2strings(genes,genNumber,individual,actCycleCounter,actuators)
    
    % find wich actuators has a positive gene
    motorIndexes =  genes(genNumber,individual,actCycleCounter,:) == 1 ;
    % deduce the corresponding strings to actuate
    stringIndexes = actuators(:,motorIndexes);
    % reshape to the output format
    stringsToActuate = reshape(stringIndexes,1,length(find(motorIndexes==1)));

end