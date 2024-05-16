%% Plot by questionnaire
conditions = {'LiveFast', 'LiveSlow', 'RecordedFast', 'RecordedSlow'};
%for every participant: 3 values (condition, questionnaire score, PLV)
allQs = [avgLiveFast' avgLiveSlow' avgRecFast' avgRecSlow']';
allPLVs = [squeeze(aggregatePLVs(:,1,16))' squeeze(aggregatePLVs(:,2,16))' squeeze(aggregatePLVs(:,3,16))' squeeze(aggregatePLVs(:,4,16))']' 
allconditions = [ones(1,21) 2*ones(1,21) 3*ones(1,21) 4*ones(1,21)]'

scatter3(allconditions,allQs,allPLVs);
xticks([1 2 3 4])
xticklabels(conditions);
xlabel("Condition")
yticks([1 2 3 4 5])
ylabel("Avg questionairre score")
zlabel("PLV at 3.2 Hz")
%%
qNum = 5;
allQs = [aggregate(:,1,qNum)' aggregate(:,2,qNum)' aggregate(:,3,qNum)' aggregate(:,4,qNum)']';
allPLVs = [squeeze(aggregatePLVs(:,1,16))' squeeze(aggregatePLVs(:,2,16))' squeeze(aggregatePLVs(:,3,16))' squeeze(aggregatePLVs(:,4,16))']' 
allconditions = [ones(1,21) 2*ones(1,21) 3*ones(1,21) 4*ones(1,21)]'

scatter3(allconditions,allQs,allPLVs);
xticks([1 2 3 4])
xticklabels(conditions);
xlabel("Condition")
yticks([1 2 3 4 5])
ylabel("Avg questionairre score")
zlabel("PLV at 3.2 Hz")

figure;
scatter(allQs,allPLVs)
lsline
%%
[P,S] = polyfit(allQs,allPLVs,1)
hold on;
prediction = P(1)*allQs+P(2);
Rsq = 1 - sum((allPLVs - prediction).^2)/sum((allPLVs - mean(allPLVs)).^2)
%% Format data in table
participant = nan(84,1);
performanceCondition = repmat({nan},84,1);
tempoCondition = repmat({nan},84,1);
score = nan(84,1);
qNum = 1:7;
for i = 1:21
    participant((i*4-3):i*4) = [i i i i];
    performanceCondition((i*4-3):i*4) = {'Live' 'Live' 'Recorded' 'Recorded'};
    tempoCondition((i*4-3):i*4) = {'Fast' 'Slow' 'Fast' 'Slow'};
    score((i*4-3):i*4) = [squeeze(mean(aggregate(i,1,qNum),3)) squeeze(mean(aggregate(i,2,qNum),3)) squeeze(mean(aggregate(i,3,qNum),3)) squeeze(mean(aggregate(i,4,qNum),3))];
end

data = table(participant, performanceCondition, tempoCondition,score);

lme = fitlme(data,'score ~ tempoCondition+performanceCondition+tempoCondition:performanceCondition+(1|participant)');
% can also consider 'score ~ tempoCondition*performanceCondition+(tempoCondition*performanceCondition|participant)'