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

function plotEvolution(performance,Fitness,nbActuators,nbGeneration,nbIndividuals,nbActuations,k,p)

    generationsVector = 1:nbGeneration;
    avgPerformance = mean(performance);
    
    figure();
    plot(generationsVector, avgPerformance)
    hold on
    boxplot(performance);
    title([num2str(nbActuations), ' actuation cycles , ', num2str(nbIndividuals), ' individuals, k = ',num2str(k),', p = ',num2str(p),', ' num2str(nbActuators),' init actuators']);
    xlabel('Generation');
    ylabel(['Performance: ', Fitness, ' (m)']);
    grid on
    
end

% tip: ['Temperature is ',num2str(c),' C']