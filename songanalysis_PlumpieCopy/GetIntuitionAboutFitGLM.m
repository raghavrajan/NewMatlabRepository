function [] = GetIntuitionAboutFitGLM()

% A script to understand different relations between variables and how that
% comes out in a GLM model

% Ok, let me first define the variables. I'm going to mimic the conditions
% I have for Harini's data analysis, where I analyze different song bout
% parameters as the response variable and two predictor variables, one of
% which is the distance from the female, a continuous predictor variable
% but recorded at discrete distances. The other predictor variable is a
% categorical variable for the context of the song - directed or undirected
% First I will try response variables that are discrete integers

NumTrials = 50;
Distances = [0 20 60 110 165 200]; % representing distance from femae in cm
Context = [1 2]; % representing DIR and UNDIR
INNumRange = 10;
Colours = 'rb';

figure;
% Response variables is # of INs - so an integer >= 1
% Predictor variable 1 - distance from female using 6 distances

% Case 1: Song does not change in # of INs for either context and for any
% of the distances. #INs is normally distributed with mean of 5 and std of
% 1.5 for all distances and contexts
INNum_Mean = 5;
INNum_STD = 1.5;

for j = Context(:)',
    Means_BootCIs{j} = zeros(length(Distances), 2);
end
Response_PredictorVar = [];
for i = 1:length(Distances),

    for j = Context(:)',
        ResponseVar = ceil(normrnd(INNum_Mean, INNum_STD, NumTrials, 1));
        ResponseVar(find(ResponseVar < 1)) = 1;
        Response_PredictorVar = [Response_PredictorVar; [ResponseVar(:) ones(NumTrials,1)*i ones(NumTrials,1)*j]];
        Means(j,i) = mean(ResponseVar);
        Means_BootCIs{j}(i,:) = bootci(10000, @mean, ResponseVar);
    end
end
subplot(1,3,1);
hold on;
for i = 1:length(Context),
    plot(Distances, Means(i,:), [Colours(i), 's-'], 'LineWidth', 1);
    plot(repmat(Distances(:), 1, 2)', Means_BootCIs{i}', Colours(i))
end
xlabel('Distance from female (cm)');
ylabel('# of INs');
axis auto;

subplot(1,3,2);
hold on;
for i = 1:length(Distances),
    Indices = find(Response_PredictorVar(:,2) == i);
    DistanceMeans(i) = mean(Response_PredictorVar(Indices,1));
    DistanceCIs(i,:) = bootci(10000, @mean, Response_PredictorVar(Indices,1));
end
plot(Distances, DistanceMeans, 'ko');
plot(repmat(Distances(:), 1, 2)', DistanceCIs', 'k');
xlabel('Distance from female (cm)');
ylabel('# of INs');
title('Case 1');
axis auto;

subplot(1,3,3);
hold on;
for i = 1:length(Context),
    Indices = find(Response_PredictorVar(:,3) == i);
    ContextMeans(i) = mean(Response_PredictorVar(Indices,1));
    ContextCIs(i,:) = bootci(10000, @mean, Response_PredictorVar(Indices,1));
end
plot(Context, ContextMeans, 'ko');
plot(repmat(Context(:), 1, 2)', ContextCIs', 'k');
xlabel('Context');
ylabel('# of INs');
axis auto;

disp('Case 1: # of INs same with distance and context');
mdl = fitglm([Distances(Response_PredictorVar(:,2))' Response_PredictorVar(:,3)], Response_PredictorVar(:,1), 'interactions', 'CategoricalVars', [2], 'VarNames', {'Distance', 'Context', '# of INs'}, 'Distribution', 'poisson')


% Case 2: Song changes in # of INs only for context, but does not change
% with distance. With context, the mean is 1.5 times higher.

clear Means_BootCIs Response_PredictorVar Means
INNum_Mean = 5;
INNum_STD = 1.5;

for j = Context(:)',
    Means_BootCIs{j} = zeros(length(Distances), 2);
end
Response_PredictorVar = [];
for i = 1:length(Distances),
    for j = Context(:)',
        if (j == 1)
            ResponseVar = ceil(normrnd(INNum_Mean*1.5, INNum_STD, NumTrials, 1));
        else
            ResponseVar = ceil(normrnd(INNum_Mean, INNum_STD, NumTrials, 1));
        end
        ResponseVar(find(ResponseVar < 1)) = 1;
        Response_PredictorVar = [Response_PredictorVar; [ResponseVar(:) ones(NumTrials,1)*i ones(NumTrials,1)*j]];
        Means(j,i) = mean(ResponseVar);
        Means_BootCIs{j}(i,:) = bootci(10000, @mean, ResponseVar);
    end
end
figure;
subplot(1,3,1);
hold on;
for i = 1:length(Context),
    plot(Distances, Means(i,:), [Colours(i), 's-'], 'LineWidth', 1);
    plot(repmat(Distances(:), 1, 2)', Means_BootCIs{i}', Colours(i))
end
xlabel('Distance from female (cm)');
ylabel('# of INs');
axis auto;

subplot(1,3,2);
hold on;
for i = 1:length(Distances),
    Indices = find(Response_PredictorVar(:,2) == i);
    DistanceMeans(i) = mean(Response_PredictorVar(Indices,1));
    DistanceCIs(i,:) = bootci(10000, @mean, Response_PredictorVar(Indices,1));
end
plot(Distances, DistanceMeans, 'ko');
plot(repmat(Distances(:), 1, 2)', DistanceCIs', 'k');
xlabel('Distance from female (cm)');
ylabel('# of INs');
title('Case 2');
axis auto;

subplot(1,3,3);
hold on;
for i = 1:length(Context),
    Indices = find(Response_PredictorVar(:,3) == i);
    ContextMeans(i) = mean(Response_PredictorVar(Indices,1));
    ContextCIs(i,:) = bootci(10000, @mean, Response_PredictorVar(Indices,1));
end
plot(Context, ContextMeans, 'ko');
plot(repmat(Context(:), 1, 2)', ContextCIs', 'k');
xlabel('Context');
ylabel('# of INs');
axis auto;

disp('Case 2: # of INs same with distance but mean is 1.5 times greater in one context');
mdl = fitglm([Distances(Response_PredictorVar(:,2))' Response_PredictorVar(:,3)], Response_PredictorVar(:,1), 'interactions', 'CategoricalVars', [2], 'VarNames', {'Distance', 'Context', '# of INs'}, 'Distribution', 'poisson')

% Case 3: Song changes linearly with distance but is same across contexts

clear Means_BootCIs Response_PredictorVar Means
INNum_Mean = 8;
INNum_STD = 1.5;

for j = Context(:)',
    Means_BootCIs{j} = zeros(length(Distances), 2);
end

Response_PredictorVar = [];
for i = 1:length(Distances),
    for j = Context(:)',
        ResponseVar = ceil(normrnd(INNum_Mean - i*0.5, INNum_STD, NumTrials, 1));
        ResponseVar(find(ResponseVar < 1)) = 1;
        Response_PredictorVar = [Response_PredictorVar; [ResponseVar(:) ones(NumTrials,1)*i ones(NumTrials,1)*j]];
        Means(j,i) = mean(ResponseVar);
        Means_BootCIs{j}(i,:) = bootci(10000, @mean, ResponseVar);
    end
end
figure;
subplot(1,3,1);
hold on;
for i = 1:length(Context),
    plot(Distances, Means(i,:), [Colours(i), 's-'], 'LineWidth', 1);
    plot(repmat(Distances(:), 1, 2)', Means_BootCIs{i}', Colours(i))
end
xlabel('Distance from female (cm)');
ylabel('# of INs');
axis auto;

subplot(1,3,2);
hold on;
for i = 1:length(Distances),
    Indices = find(Response_PredictorVar(:,2) == i);
    DistanceMeans(i) = mean(Response_PredictorVar(Indices,1));
    DistanceCIs(i,:) = bootci(10000, @mean, Response_PredictorVar(Indices,1));
end
plot(Distances, DistanceMeans, 'ko');
plot(repmat(Distances(:), 1, 2)', DistanceCIs', 'k');
xlabel('Distance from female (cm)');
ylabel('# of INs');
title('Case 3');
axis auto;

subplot(1,3,3);
hold on;
for i = 1:length(Context),
    Indices = find(Response_PredictorVar(:,3) == i);
    ContextMeans(i) = mean(Response_PredictorVar(Indices,1));
    ContextCIs(i,:) = bootci(10000, @mean, Response_PredictorVar(Indices,1));
end
plot(Context, ContextMeans, 'ko');
plot(repmat(Context(:), 1, 2)', ContextCIs', 'k');
xlabel('Context');
ylabel('# of INs');
axis auto;

NumINs = Response_PredictorVar(:,1);
Distance = Distances(Response_PredictorVar(:,2));
Context_Values = Response_PredictorVar(:,3);
disp('Case 3: # of INs same with context but decreases by 0.5 for each successive distance');
mdl = fitglm([Distances(Response_PredictorVar(:,2))' Response_PredictorVar(:,3)], Response_PredictorVar(:,1), 'interactions', 'CategoricalVars', [2], 'VarNames', {'Distance', 'Context_Values', 'NumINs'}, 'Distribution', 'poisson')

% Case 4: Song changes linearly with distance only in one context - however
% means are same at the first distance across context

clear Means_BootCIs Response_PredictorVar Means
INNum_Mean = 8;
INNum_STD = 1.5;

for j = Context(:)',
    Means_BootCIs{j} = zeros(length(Distances), 2);
end

Response_PredictorVar = [];
for i = 1:length(Distances),
    for j = Context(:)',
        if (j == 1)
            ResponseVar = ceil(normrnd(INNum_Mean - i*0.5, INNum_STD, NumTrials, 1));
        else
            ResponseVar = ceil(normrnd(INNum_Mean, INNum_STD, NumTrials, 1));
        end
        ResponseVar(find(ResponseVar < 1)) = 1;
        Response_PredictorVar = [Response_PredictorVar; [ResponseVar(:) ones(NumTrials,1)*i ones(NumTrials,1)*j]];
        Means(j,i) = mean(ResponseVar);
        Means_BootCIs{j}(i,:) = bootci(10000, @mean, ResponseVar);
    end
end
figure;
subplot(1,3,1);
hold on;
for i = 1:length(Context),
    plot(Distances, Means(i,:), [Colours(i), 's-'], 'LineWidth', 1);
    plot(repmat(Distances(:), 1, 2)', Means_BootCIs{i}', Colours(i))
end
xlabel('Distance from female (cm)');
ylabel('# of INs');
axis auto;

subplot(1,3,2);
hold on;
for i = 1:length(Distances),
    Indices = find(Response_PredictorVar(:,2) == i);
    DistanceMeans(i) = mean(Response_PredictorVar(Indices,1));
    DistanceCIs(i,:) = bootci(10000, @mean, Response_PredictorVar(Indices,1));
end
plot(Distances, DistanceMeans, 'ko');
plot(repmat(Distances(:), 1, 2)', DistanceCIs', 'k');
xlabel('Distance from female (cm)');
ylabel('# of INs');
title('Case 4');
axis auto;

subplot(1,3,3);
hold on;
for i = 1:length(Context),
    Indices = find(Response_PredictorVar(:,3) == i);
    ContextMeans(i) = mean(Response_PredictorVar(Indices,1));
    ContextCIs(i,:) = bootci(10000, @mean, Response_PredictorVar(Indices,1));
end
plot(Context, ContextMeans, 'ko');
plot(repmat(Context(:), 1, 2)', ContextCIs', 'k');
xlabel('Context');
ylabel('# of INs');
axis auto;

NumINs = Response_PredictorVar(:,1);
Distance = Distances(Response_PredictorVar(:,2));
Context_Values = Response_PredictorVar(:,3);
disp('Case 4: # of INs decreases by 0.5 for each successive distance in context 1, means same at distance 1 across context');
mdl = fitglm([Distances(Response_PredictorVar(:,2))' Response_PredictorVar(:,3)], Response_PredictorVar(:,1), 'interactions', 'CategoricalVars', [2], 'VarNames', {'Distance', 'Context_Values', 'NumINs'}, 'Distribution', 'poisson')

% Case 5: Song changes linearly with distance and also there is a
% difference in means across context at the shortest distance

clear Means_BootCIs Response_PredictorVar Means
INNum_Mean = 8;
INNum_STD = 1.5;

for j = Context(:)',
    Means_BootCIs{j} = zeros(length(Distances), 2);
end

Response_PredictorVar = [];
for i = 1:length(Distances),
    for j = Context(:)',
        if (j == 1)
            ResponseVar = ceil(normrnd(INNum_Mean - i*0.5, INNum_STD, NumTrials, 1));
        else
            ResponseVar = ceil(normrnd(INNum_Mean-3, INNum_STD, NumTrials, 1));
        end
        ResponseVar(find(ResponseVar < 1)) = 1;
        Response_PredictorVar = [Response_PredictorVar; [ResponseVar(:) ones(NumTrials,1)*i ones(NumTrials,1)*j]];
        Means(j,i) = mean(ResponseVar);
        Means_BootCIs{j}(i,:) = bootci(10000, @mean, ResponseVar);
    end
end
figure;
subplot(1,3,1);
hold on;
for i = 1:length(Context),
    plot(Distances, Means(i,:), [Colours(i), 's-'], 'LineWidth', 1);
    plot(repmat(Distances(:), 1, 2)', Means_BootCIs{i}', Colours(i))
end
xlabel('Distance from female (cm)');
ylabel('# of INs');
axis auto;

subplot(1,3,2);
hold on;
for i = 1:length(Distances),
    Indices = find(Response_PredictorVar(:,2) == i);
    DistanceMeans(i) = mean(Response_PredictorVar(Indices,1));
    DistanceCIs(i,:) = bootci(10000, @mean, Response_PredictorVar(Indices,1));
end
plot(Distances, DistanceMeans, 'ko');
plot(repmat(Distances(:), 1, 2)', DistanceCIs', 'k');
xlabel('Distance from female (cm)');
ylabel('# of INs');
title('Case 5');
axis auto;

subplot(1,3,3);
hold on;
for i = 1:length(Context),
    Indices = find(Response_PredictorVar(:,3) == i);
    ContextMeans(i) = mean(Response_PredictorVar(Indices,1));
    ContextCIs(i,:) = bootci(10000, @mean, Response_PredictorVar(Indices,1));
end
plot(Context, ContextMeans, 'ko');
plot(repmat(Context(:), 1, 2)', ContextCIs', 'k');
xlabel('Context');
ylabel('# of INs');
axis auto;

NumINs = Response_PredictorVar(:,1);
Distance = Distances(Response_PredictorVar(:,2));
Context_Values = Response_PredictorVar(:,3);
disp('Case 5: # of INs decreases by 0.5 for each successive distance in context 1, and means are different clat distance 1 across context');
mdl = fitglm([Distances(Response_PredictorVar(:,2))' Response_PredictorVar(:,3)], Response_PredictorVar(:,1), 'interactions', 'CategoricalVars', [2], 'VarNames', {'Distance', 'Context_Values', 'NumINs'}, 'Distribution', 'poisson')

% Case 6: Song is different between the contexts and in both contexts, the
% song changes with distances to the same extent

clear Means_BootCIs Response_PredictorVar Means
INNum_Mean = 8;
INNum_STD = 1.5;

for j = Context(:)',
    Means_BootCIs{j} = zeros(length(Distances), 2);
end

Response_PredictorVar = [];
for i = 1:length(Distances),
    for j = Context(:)',
        if (j == 1)
            ResponseVar = ceil(normrnd(INNum_Mean - i*0.5, INNum_STD, NumTrials, 1));
        else
            ResponseVar = ceil(normrnd(INNum_Mean/1.5 - i*0.5, INNum_STD, NumTrials, 1));
        end
        ResponseVar(find(ResponseVar < 1)) = 1;
        Response_PredictorVar = [Response_PredictorVar; [ResponseVar(:) ones(NumTrials,1)*i ones(NumTrials,1)*j]];
        Means(j,i) = mean(ResponseVar);
        Means_BootCIs{j}(i,:) = bootci(10000, @mean, ResponseVar);
    end
end
figure;
subplot(1,3,1);
hold on;
for i = 1:length(Context),
    plot(Distances, Means(i,:), [Colours(i), 's-'], 'LineWidth', 1);
    plot(repmat(Distances(:), 1, 2)', Means_BootCIs{i}', Colours(i))
end
xlabel('Distance from female (cm)');
ylabel('# of INs');
axis auto;

subplot(1,3,2);
hold on;
for i = 1:length(Distances),
    Indices = find(Response_PredictorVar(:,2) == i);
    DistanceMeans(i) = mean(Response_PredictorVar(Indices,1));
    DistanceCIs(i,:) = bootci(10000, @mean, Response_PredictorVar(Indices,1));
end
plot(Distances, DistanceMeans, 'ko');
plot(repmat(Distances(:), 1, 2)', DistanceCIs', 'k');
xlabel('Distance from female (cm)');
ylabel('# of INs');
title('Case 6');
axis auto;

subplot(1,3,3);
hold on;
for i = 1:length(Context),
    Indices = find(Response_PredictorVar(:,3) == i);
    ContextMeans(i) = mean(Response_PredictorVar(Indices,1));
    ContextCIs(i,:) = bootci(10000, @mean, Response_PredictorVar(Indices,1));
end
plot(Context, ContextMeans, 'ko');
plot(repmat(Context(:), 1, 2)', ContextCIs', 'k');
xlabel('Context');
ylabel('# of INs');
axis auto;

NumINs = Response_PredictorVar(:,1);
Distance = Distances(Response_PredictorVar(:,2));
Context_Values = Response_PredictorVar(:,3);
disp('Case 5: # of INs decreases by 0.5 for each successive distance in context 1, and means are different clat distance 1 across context');
mdl = fitglm([Distances(Response_PredictorVar(:,2))' Response_PredictorVar(:,3)], Response_PredictorVar(:,1), 'interactions', 'CategoricalVars', [2], 'VarNames', {'Distance', 'Context_Values', 'NumINs'}, 'Distribution', 'poisson')

% Case 7: Song is different between the contexts and in both contexts, it
% decreases with distance, but the rate of decrease is different for
% different contexts

clear Means_BootCIs Response_PredictorVar Means
INNum_Mean = 12;
INNum_STD = 1.5;

for j = Context(:)',
    Means_BootCIs{j} = zeros(length(Distances), 2);
end

Response_PredictorVar = [];
for i = 1:length(Distances),
    for j = Context(:)',
        if (j == 1)
            ResponseVar = ceil(normrnd(INNum_Mean - i*1.5, INNum_STD, NumTrials, 1));
        else
            ResponseVar = ceil(normrnd(INNum_Mean/1.33 - i*0.5, INNum_STD, NumTrials, 1));
        end
        ResponseVar(find(ResponseVar < 1)) = 1;
        Response_PredictorVar = [Response_PredictorVar; [ResponseVar(:) ones(NumTrials,1)*i ones(NumTrials,1)*j]];
        Means(j,i) = mean(ResponseVar);
        Means_BootCIs{j}(i,:) = bootci(10000, @mean, ResponseVar);
    end
end
figure;
subplot(1,3,1);
hold on;
for i = 1:length(Context),
    plot(Distances, Means(i,:), [Colours(i), 's-'], 'LineWidth', 1);
    plot(repmat(Distances(:), 1, 2)', Means_BootCIs{i}', Colours(i))
end
xlabel('Distance from female (cm)');
ylabel('# of INs');
axis auto;

subplot(1,3,2);
hold on;
for i = 1:length(Distances),
    Indices = find(Response_PredictorVar(:,2) == i);
    DistanceMeans(i) = mean(Response_PredictorVar(Indices,1));
    DistanceCIs(i,:) = bootci(10000, @mean, Response_PredictorVar(Indices,1));
end
plot(Distances, DistanceMeans, 'ko');
plot(repmat(Distances(:), 1, 2)', DistanceCIs', 'k');
xlabel('Distance from female (cm)');
ylabel('# of INs');
title('Case 7');
axis auto;

subplot(1,3,3);
hold on;
for i = 1:length(Context),
    Indices = find(Response_PredictorVar(:,3) == i);
    ContextMeans(i) = mean(Response_PredictorVar(Indices,1));
    ContextCIs(i,:) = bootci(10000, @mean, Response_PredictorVar(Indices,1));
end
plot(Context, ContextMeans, 'ko');
plot(repmat(Context(:), 1, 2)', ContextCIs', 'k');
xlabel('Context');
ylabel('# of INs');
axis auto;

NumINs = Response_PredictorVar(:,1);
Distance = Distances(Response_PredictorVar(:,2));
Context_Values = Response_PredictorVar(:,3);
disp('Case 7: # of INs decreases by different rates for each successive distance in context 1, and means are different clat distance 1 across context');
mdl = fitglm([Distances(Response_PredictorVar(:,2))' Response_PredictorVar(:,3)], Response_PredictorVar(:,1), 'interactions', 'CategoricalVars', [2], 'VarNames', {'Distance', 'Context_Values', 'NumINs'}, 'Distribution', 'poisson')

% Case 7: Song is different between the contexts and in both contexts, it
% decreases with distance, but the rate of decrease is different for
% different contexts and we don't have values for all distances for both
% contexts. For Context 1, we have only the first 3 distances and for
% context 2, we have the last 4 distances - so basically one point of
% overlap. I will try to get the means to be equal at the point of overlap

clear Means_BootCIs Response_PredictorVar Means
INNum_Mean = 12;
INNum_STD = 1.5;

for j = Context(:)',
    Means_BootCIs{j} = zeros(length(Distances), 2);
end

Response_PredictorVar = [];
for i = 1:length(Distances),
    for j = Context(:)',
        if (j == 1)
            if (i < 4)
                ResponseVar = ceil(normrnd(INNum_Mean - i*1.5, INNum_STD, NumTrials, 1));
            else
                ResponseVar = [];
            end
        else
            if (i < 4)
                ResponseVar = [];
            else
                ResponseVar = ceil(normrnd(INNum_Mean/1.33 - i*0.5, INNum_STD, NumTrials, 1));
            end
        end
        ResponseVar(find(ResponseVar < 1)) = 1;
        if (~isempty(ResponseVar))
            Response_PredictorVar = [Response_PredictorVar; [ResponseVar(:) ones(NumTrials,1)*i ones(NumTrials,1)*j]];
            Means(j,i) = mean(ResponseVar);
            Means_BootCIs{j}(i,:) = bootci(10000, @mean, ResponseVar);
        else
            Means(j,i) = NaN;
            Means_BootCIs{j}(i,:) = ones(1,2)*NaN;
        end
    end
end
figure;
subplot(1,3,1);
hold on;
for i = 1:length(Context),
    plot(Distances, Means(i,:), [Colours(i), 's-'], 'LineWidth', 1);
    plot(repmat(Distances(:), 1, 2)', Means_BootCIs{i}', Colours(i))
end
xlabel('Distance from female (cm)');
ylabel('# of INs');
axis auto;

subplot(1,3,2);
hold on;
for i = 1:length(Distances),
    Indices = find(Response_PredictorVar(:,2) == i);
    DistanceMeans(i) = mean(Response_PredictorVar(Indices,1));
    DistanceCIs(i,:) = bootci(10000, @mean, Response_PredictorVar(Indices,1));
end
plot(Distances, DistanceMeans, 'ko');
plot(repmat(Distances(:), 1, 2)', DistanceCIs', 'k');
xlabel('Distance from female (cm)');
ylabel('# of INs');
title('Case 8');
axis auto;

subplot(1,3,3);
hold on;
for i = 1:length(Context),
    Indices = find(Response_PredictorVar(:,3) == i);
    ContextMeans(i) = mean(Response_PredictorVar(Indices,1));
    ContextCIs(i,:) = bootci(10000, @mean, Response_PredictorVar(Indices,1));
end
plot(Context, ContextMeans, 'ko');
plot(repmat(Context(:), 1, 2)', ContextCIs', 'k');
xlabel('Context');
ylabel('# of INs');
axis auto;

NumINs = Response_PredictorVar(:,1);
Distance = Distances(Response_PredictorVar(:,2));
Context_Values = Response_PredictorVar(:,3);
disp('Case 8: # of INs decreases by different rates for each successive distance in context 1, and means are different clat distance 1 across context');
mdl = fitglm([Distances(Response_PredictorVar(:,2))' Response_PredictorVar(:,3)], Response_PredictorVar(:,1), 'interactions', 'CategoricalVars', [2], 'VarNames', {'Distance', 'Context_Values', 'NumINs'}, 'Distribution', 'poisson')

% Case 7: Song is different between the contexts and in both contexts, it
% decreases with distance, but the rate of decrease is different for
% different contexts and we don't have values for all distances for both
% contexts. For Context 1, we have only the first 3 distances and for
% context 2, we have the last 4 distances - so basically one point of
% overlap. I will try to get the means to be equal at the point of overlap

clear Means_BootCIs Response_PredictorVar Means
INNum_Mean = 12;
INNum_STD = 1.5;

for j = Context(:)',
    Means_BootCIs{j} = zeros(length(Distances), 2);
end

Response_PredictorVar = [];
for i = 1:length(Distances),
    for j = Context(:)',
        if (j == 1)
            if (i < 4)
                ResponseVar = ceil(normrnd(1.25*INNum_Mean - i*1.5, INNum_STD, NumTrials, 1));
            else
                ResponseVar = [];
            end
        else
            if (i < 4)
                ResponseVar = [];
            else
                ResponseVar = ceil(normrnd(INNum_Mean/1.33 - i*0.5, INNum_STD, NumTrials, 1));
            end
        end
        ResponseVar(find(ResponseVar < 1)) = 1;
        if (~isempty(ResponseVar))
            Response_PredictorVar = [Response_PredictorVar; [ResponseVar(:) ones(NumTrials,1)*i ones(NumTrials,1)*j]];
            Means(j,i) = mean(ResponseVar);
            Means_BootCIs{j}(i,:) = bootci(10000, @mean, ResponseVar);
        else
            Means(j,i) = NaN;
            Means_BootCIs{j}(i,:) = ones(1,2)*NaN;
        end
    end
end
figure;
subplot(1,3,1);
hold on;
for i = 1:length(Context),
    plot(Distances, Means(i,:), [Colours(i), 's-'], 'LineWidth', 1);
    plot(repmat(Distances(:), 1, 2)', Means_BootCIs{i}', Colours(i))
end
xlabel('Distance from female (cm)');
ylabel('# of INs');
axis auto;

subplot(1,3,2);
hold on;
for i = 1:length(Distances),
    Indices = find(Response_PredictorVar(:,2) == i);
    DistanceMeans(i) = mean(Response_PredictorVar(Indices,1));
    DistanceCIs(i,:) = bootci(10000, @mean, Response_PredictorVar(Indices,1));
end
plot(Distances, DistanceMeans, 'ko');
plot(repmat(Distances(:), 1, 2)', DistanceCIs', 'k');
xlabel('Distance from female (cm)');
ylabel('# of INs');
title('Case 9');
axis auto;

subplot(1,3,3);
hold on;
for i = 1:length(Context),
    Indices = find(Response_PredictorVar(:,3) == i);
    ContextMeans(i) = mean(Response_PredictorVar(Indices,1));
    ContextCIs(i,:) = bootci(10000, @mean, Response_PredictorVar(Indices,1));
end
plot(Context, ContextMeans, 'ko');
plot(repmat(Context(:), 1, 2)', ContextCIs', 'k');
xlabel('Context');
ylabel('# of INs');
axis auto;

NumINs = Response_PredictorVar(:,1);
Distance = Distances(Response_PredictorVar(:,2));
Context_Values = Response_PredictorVar(:,3);
disp('Case 9: # of INs decreases by different rates for each successive distance in context 1, and means are different clat distance 1 across context');
mdl = fitglm([Distances(Response_PredictorVar(:,2))' Response_PredictorVar(:,3)], Response_PredictorVar(:,1), 'interactions', 'CategoricalVars', [2], 'VarNames', {'Distance', 'Context_Values', 'NumINs'}, 'Distribution', 'poisson')

disp('Finished');

