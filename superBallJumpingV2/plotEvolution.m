% Function that plot the evolution of the performance with respect to the
% generation number
%
% Inputs: - performance (nbIndiv x nbGeneration matric containing perf value)
%         - Fitness (Performance to be evaluated: jump, dist or jumpxdist)
%         - nbActuators (at the same time)
%         - nbGeneration 
%         - nbIndividuals (population size)
%         - nbActuations (in the sequence)
%         - k (slection parameter, nb of bad individuals replaced each gen)
%
% Output: - plot of avg perf over the generations
%         - box plot of perf over the generations

function plotEvolution(performance,Fitness,nbActuators,nbGeneration,nbIndividuals,nbActuations,k)

    generationsVector = 1:nbGeneration;
    avgPerformance = mean(performance);
    
    figure();
    plot(generationsVector, avgPerformance)
    hold on
    boxplot(performance);
    title(['Evolution parameters: ', num2str(nbActuators),' actuators, ', num2str(nbActuations), ' actuations, ', num2str(nbIndividuals), ' individuals, Selection parameter = ',num2str(k)]);
    xlabel('Generation');
    ylabel(['Performance: ', Fitness, ' (m)']);
    grid on
    
end

% tip: ['Temperature is ',num2str(c),' C']