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

function plotEvolution(performance,Fitness,nbActuators,nbGeneration,nbIndividuals,nbActuations,t,p,nbAvgActuatorsBestIndiv)

    generationsVector = 1:nbGeneration;
    avgPerformance = mean(performance);
    if strcmpi(Fitness,'DistToGoal')
        bestPerformance = min(performance);
    else    
        bestPerformance = max(performance);
    end
    
    figure();
    yyaxis left
    plot(generationsVector, bestPerformance,'--','LineWidth',1.5)
    hold on
    plot(generationsVector, avgPerformance,'r','LineWidth',1.5)
    boxplot(performance);
    title([num2str(nbIndividuals), ' Indiv., ' ,num2str(nbActuators), ' actuators, ' ,num2str(nbActuations), ' actuation cycle(s), ', 'Tournmt size = ',num2str(t),', p_{mut} = ',num2str(p)]);
    xlabel('Generation');
    ylabel(['Performance: ', Fitness, ' (m)']);
    yyaxis right
    plot(generationsVector,nbAvgActuatorsBestIndiv,'--','LineWidth',1)
    ylabel('Average number of motors per cycle (Best Indiv.)')
    grid on
    legend('Best performance','Avg performance','Nmb of motors')
    set(gca,'fontsize', 14);
    
end