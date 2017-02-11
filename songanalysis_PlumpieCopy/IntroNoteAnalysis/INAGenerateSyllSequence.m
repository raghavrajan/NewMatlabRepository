function [Labels] = INAGenerateSyllSequence(TransProb, Syllables, NoofTrials)

for i = 1:NoofTrials,
    SimSyll = 1;
    SimSyllNo = 1;
    while (SimSyll < (length(Syllables)))
        PossibleStates = find(TransProb(SimSyll,:));
        if (length(PossibleStates) == 1)
            SimSyll = PossibleStates;
            Labels{i}(SimSyllNo) = Syllables(SimSyll);
            SimSyllNo = SimSyllNo + 1;
        else
            Probs = [0 cumsum(TransProb(SimSyll, PossibleStates))];
            RandNo = rand(1);
            for j = 2:length(Probs),
                if ((RandNo >= Probs(j-1)) && (RandNo < Probs(j)))
                    SimSyll = PossibleStates(j-1);
                    Labels{i}(SimSyllNo) = Syllables(SimSyll);
                    SimSyllNo = SimSyllNo + 1;
                    break;
                end
            end
        end
    end
end
