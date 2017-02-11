function [Position, NoofINs, INFR, Stats] = IntroNoteFRPreMotorWindowAnalysis(Neural_INR, PreMotorLatencies, PlotOption)

for Latency = 1:length(PreMotorLatencies),
    
    PreMotorLatency = PreMotorLatencies(Latency);
    INIndex = 0;

    for i = 1:length(Neural_INR.NoofINs),
        if (Neural_INR.NoofINs(i) > 0)
            BoutSpikeTimes = Neural_INR.BoutDetails(i).SpikeTimes;
            INs = Neural_INR.INs{i};
            for j = 1:length(INs),
                INIndex = INIndex + 1;

                Position{Latency}(INIndex,1) = j - length(INs) - 1;
                Position{Latency}(INIndex,2) = j;
                NoofINs{Latency}(INIndex,1) = length(INs);

                INOnset = Neural_INR.BoutDetails(i).onsets(INs(j)) - PreMotorLatency;
                INOffset = Neural_INR.BoutDetails(i).offsets(INs(j)) - PreMotorLatency;

                INSpikeIndices = find((BoutSpikeTimes >= INOnset) & (BoutSpikeTimes < INOffset));
                INFR{Latency}(INIndex, 1) = length(INSpikeIndices)/(INOffset - INOnset);
            end
        end
    end

    for i = 1:size(Neural_INR.WithinBoutNoofINs, 1),
        if (Neural_INR.WithinBoutNoofINs(i,1) > 0)
            BoutSpikeTimes = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).SpikeTimes;
            INs = Neural_INR.WithinBoutINs{i};
            for j = 1:length(INs),
                INIndex = INIndex + 1;

                Position{Latency}(INIndex,1) = j - length(INs) - 1;
                Position{Latency}(INIndex,2) = j;
                NoofINs{Latency}(INIndex,1) = length(INs);

                INOnset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(INs(j)) - PreMotorLatency;
                INOffset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).offsets(INs(j)) - PreMotorLatency;
                
                INSpikeIndices = find((BoutSpikeTimes >= INOnset) & (BoutSpikeTimes < INOffset));
                INFR{Latency}(INIndex, 1) = length(INSpikeIndices)/(INOffset - INOnset);
            end
        end
    end
end

for i = 1:length(PreMotorLatencies),
    if (strfind(PlotOption, 'on'))
        subplot(2, round(length(PreMotorLatencies)/2), i);
        hold on;
        for j = min(Position{i}(:,1)):1:max(Position{i}(:,1)),
            errorbar(j, mean(INFR{i}(find(Position{i}(:,1) == j))), std(INFR{i}(find(Position{i}(:,1) == j))), 'ko');
        end
    end
    [r, p] = corrcoef(Position{i}(:,1), INFR{i});
    Stats(i, 1) = r(1,2);
    Stats(i, 2) = p(1,2);
    if (p(1,2) < 0.05)
        Stats(i,3) = 1;
    else
        Stats(i,3) = 0;
    end
    if (strfind(PlotOption, 'on'))
        if (p(1,2) < 0.05)
            plot(unique(Position{i}(:,1)), polyval(polyfit(Position{i}(:,1), INFR{i}, 1), unique(Position{i}(:,1))), 'r');
        else
            plot(unique(Position{i}(:,1)), polyval(polyfit(Position{i}(:,1), INFR{i}, 1), unique(Position{i}(:,1))), 'k');
        end
    end
    [Rsq, F] = CalculateGoodnessofLinearFit(polyfit(Position{i}(:,1), INFR{i}, 1), Position{i}(:,1), INFR{i});
    disp(['PreMotorLatency is ', num2str(PreMotorLatencies(i)), ' and Rsq is ', num2str(Rsq), ' and corr. coeff ^2 is ', num2str(r(1,2).^2), ' and the p-value is ', num2str(p(1,2)), ' for INs']);
    if (strfind(PlotOption, 'on'))
        axis tight;
        temp = axis;
        axis([(min(Position{i}(:,1)) - 0.5) -0.75 0 temp(4)*1.1]);
        set(gca, 'Box', 'off');
    end
end
disp('Finished Analysis');