% Author:   Gregoire Besson, EPFL student, gregoire.besson@alumni.epfl.ch
% Date:     Winter 2018-2019
%
% Copy the genome of the current genome to the next one and run a mutation 
% process on the k worst indivudual mutation here consists only on toggling
% a motor between ON and OFF
%
% Inputs:   - mutation method: hard or soft
%           - genesIn(nbGen,nbIndiv,nbActCycle,nbMotors), current genome
%           - g, generation counter
%           - nbIndiv, population size
%           - indexes, ascending sorted individuales regarding their perf
%           - AvlblActuators(g,nbInd,nbAvlAct), list of available actuators
%           - nbActuationCycle
%           - k, selection parameter: worst k indiv mutate
%           - p, probability of mutation 0 < p < 1
%           - selectionMode, ranking or tournament
%           - elitism, bool
%           - e, number of elites to be copied
%
% Output:   - genesOut(nbGen,nbIndiv,nbActCycle,nbMotors), updated genome
%             containing all previous generations and the new one

function genesOut = genesMutation(mutation,genesIn,g,nbIndiv,indexes,AvlblActuators,nbActuationCycle,k,p,selectionMode,elitism,e)

if (strcmpi(selectionMode,'ranking'))

    % copy current generation to next
    genesIn(g+1,:,:,:) = genesIn(g,:,:,:);
    
    % for each k bad individuals
    for i = 1:k
        
        % extract his genes (for all the actuation cycles)
        genesToMutate = genesIn(g,indexes(i),:,:);
        
        if strcmpi(mutation,'hard')
            % generate a random number value for each gene
            x = rand(nbActuationCycle,12);
            % all value under proba threshold fixed by user is mutated
            mutateThis = x <= p;
            genesToMutate(mutateThis) = ~genesToMutate(mutateThis);
            % put back the new genes in next generation
            genesIn(g+1,indexes(i),:,:) = genesToMutate;
        elseif strcmpi(mutation,'soft')
            for c = 1:nbActuationCycle
                % generate a random number value for each cycle
                x = rand();
                % probability of mutation of this cycle 
                if (x <= p)
                    % toggle a random gene in this cycle (between 1 to 12)
                    randomGene = 1 + round(rand()*11);
                    genesIn(g+1,i,c,randomGene) = ~genesIn(g+1,i,c,randomGene);
                end
            end       
        else 
            error('Unknown mutation method ! \n');
        end
    end
    
elseif (strcmpi(selectionMode,'tournament'))
    
    start = 1;
    [~,~,nbAvlAct] = size(AvlblActuators);
    
    % don't mutate the elites if elitism enabled
    if (elitism)
        start = e+1;
    end
    
    for i = start:nbIndiv
        % extract genes of the indiv (for all the actuation cycles)
        genesToMutate = genesIn(g+1,i,:,:);
        
        if strcmpi(mutation,'hard')
            % generate a random number value for each gene
            x = rand(nbActuationCycle,12);
            % all value under proba threshold fixed by user is mutated
            mutateThis = x <= p;
            genesToMutate(mutateThis) = ~genesToMutate(mutateThis);
            % put back the new genes in the genome
            genesIn(g+1,i,:,:) = genesToMutate;
            
        elseif strcmpi(mutation,'soft')
            for c = 1:nbActuationCycle
                % generate a random number for each cycle
                x = rand();
                % probability of mutation of this cycle 
                if (x <= p)
                    % toggle a random gene (btwn 1 to nbAvAct)
                    randomGeneIndex = 1 + round(rand()*(nbAvlAct-1));
                    randomGene = AvlblActuators(g+1,i,randomGeneIndex);
                    genesIn(g+1,i,c,randomGene) = ~genesIn(g+1,i,c,randomGene);
                end
            end
            
        elseif strcmpi(mutation,'ConstNbAct')
            for c = 1:nbActuationCycle
                % generate a random number for each cycle
                x = rand();
                % probability of mutation of this cycle 
                if (x <= p)
                    
                    %find wich actuators were actuates at that cycle
                    prevAct = find(genesIn(g+1,i,c,:)==1);
                    nbAct = length(prevAct);
                    
                    % pick one randomly to turn off
                    indexToTurnOff = 1+round(rand*(nbAct-1));
                    genesIn(g+1,i,c,prevAct(indexToTurnOff)) = 0;
                    
                    % turn on an available motor
                    indexesNotUsedAct = ~ismember(AvlblActuators(g+1,i,:),prevAct);
                    notUsedAct = AvlblActuators(g+1,i,indexesNotUsedAct);
                    nbNotUsedAct = length(notUsedAct);
                    indexToTurnOn = 1+round(rand*(nbNotUsedAct-1));
                    genesIn(g+1,i,c,notUsedAct(indexToTurnOn)) = 1;                   
                end
            end  
            
        else 
            error('Unknown mutation method ! \n');
        end
    end
else 
    error('Unknown selection method');
end

% output the mutated genome
genesOut = genesIn;
   
end