function genesOut = tournamentSelection(genesIn,nbIndiv,performance,g,t,elitism,e,sortingMode)
    
%Create the random number stream for reproducibility:
rngParameters = RandStream('mlfg6331_64','Seed','Shuffle');

start = 1;
    
% copy the best individual if elitism is enable
if (elitism)
    % find the best indiv from current generation
    if (strcmpi(sortingMode,'ascend'))                                     % TODO: what happens if e>1 ????
        [~,BestIndiv] = max(performance(:,g));
    else 
        [~,BestIndiv] = min(performance(:,g));
    end
            
    % copy it as first indivi of next generation
    genesIn(g+1,1,:,:) = genesIn(g,BestIndiv,:,:);
    
    start = e+1;
end

% need to perform as much tournament as individuals to fill the next gen
for i=start:nbIndiv

    % Draw t different random individuals
    contestantsIndexes = datasample(rngParameters,1:nbIndiv,t,'Replace',false);
    
    % select champion between this t
    if (strcmpi(sortingMode,'ascend'))
        [~,championIndex] = max(performance(contestantsIndexes,g));
    else
        [~,championIndex] = min(performance(contestantsIndexes,g));
    end
    
    % copy it in next generation
    genesIn(g+1,i,:,:) = genesIn(g,contestantsIndexes(championIndex),:,:);
   
end

% output the mutated genome
genesOut = genesIn;

end
