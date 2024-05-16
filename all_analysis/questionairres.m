clear variables;
close all;

%%% FilePath Parameters
% addpath('./EEGlab'); % add eeglab to path
study_path = '/Users/arun/Library/Mobile Documents/com~apple~CloudDocs/Documents/MINDLab/performance_study_final'; %% Only thing to update!

%%%
participant_table = readtable('Participants.xlsx'); % have it in current directory


%%
% loadSurveyResponsesParticipants(participant_table, study_path, 1);
% calculates for all participants

%% Load the data
aggregate = aggregateResponse(participant_table, study_path, 1);

%% invert scoring of question 7
aggregate(:,:,7) = 6-aggregate(:,:,7);
newaggr = aggregate;
%% remove question 3
% newaggr = cat(3, aggregate(:,:,1:2), aggregate(:,:,4:8));
% % newaggr = aggregate
%% Plot avg data
figure;
tiledlayout('flow');
avgLiveFast = squeeze(mean(newaggr(:,1,:),3));
avgLiveSlow = squeeze(mean(newaggr(:,2,:),3));
avgRecFast = squeeze(mean(newaggr(:,3,:),3));
avgRecSlow = squeeze(mean(newaggr(:,4,:),3));
nexttile;
plotConditionAverages(avgLiveFast, avgLiveSlow,avgRecFast,avgRecSlow);
nexttile;
plotTwoCondition((avgLiveFast+avgLiveSlow)/2, (avgRecFast+avgRecSlow)/2)
set(gcf,"Color","white");
fontsize(20, 'points')


% exportgraphics(gcf,strcat(study_path, '/FIGURES/Q_scores.png'),'Resolution',600,'BackgroundColor','none')
%%
LiveFast_q = squeeze(newaggr(:,1,:));
LiveSlow_q = squeeze(newaggr(:,2,:));
RecFast_q = squeeze(newaggr(:,3,:));
RecSlow_q = squeeze(newaggr(:,4,:));


%% PLOT

%% Anova
%{
Anova: 
row factor is Tempo (fast/slow)
column factor is Performance Medium (Live/Rec)
A = 
[
[LiveFast_part1 RecordedFast_part1],
[LiveFast_part2 RecordedFast_part2],
...
...
...
[LiveFast_part21 RecordedFast_part21],
[LiveSlow_part1 RecordedSlow_part1]
...
...
...
[LiveSlow_part21 RecordedSlow_part21]
]
%}
A = zeros(42,2); %% row factor is tempo f/s, column factor is live/rec
A(1:21,1) = avgLiveFast;
A(1:21,2) = avgRecFast;
A(22:42,1) = avgLiveSlow;
A(22:42,2) = avgRecSlow;

[~,~,stats] = anova2(A,21); figure;
c1 = multcompare(stats); figure;
c2 = multcompare(stats,"Estimate","row");
% tbl1 = array2table(c1,"VariableNames", ...
%     ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
%%
questionA = 1;
questionB = 5;

collection1 = [squeeze(aggregate(:,1,questionA))', squeeze(aggregate(:,2,questionA))', squeeze(aggregate(:,3,questionA))', squeeze(aggregate(:,4,questionA))'];
collection2 = [squeeze(aggregate(:,1,questionB))', squeeze(aggregate(:,2,questionB))', squeeze(aggregate(:,3,questionB))', squeeze(aggregate(:,4,questionB))'];
collection3 = [squeeze(mean(aggregate(:,1,:),3))', squeeze(mean(aggregate(:,2,:),3))', squeeze(mean(aggregate(:,3,:),3))', squeeze(mean(aggregate(:,4,:),3))'];
figure;scatter(collection1,collection3);



%% Plot each question
fig = figure();
plotAllQuestions(LiveFast_q, LiveSlow_q,RecFast_q,RecSlow_q)
fontsize(fig, 20, "points")
set(gcf,"Color","white");


%% Plot all questions, only Live v Rec
fig = figure();
plotAllTwoCondition((LiveFast_q+LiveSlow_q)/2,(RecFast_q+RecSlow_q)/2)
fontsize(fig, 20, "points")
set(gcf,"Color","white");
% exportgraphics(gcf,strcat(study_path, '/FIGURES/All_Qs_Live-Rec.png'),'Resolution',600,'BackgroundColor','none')

%% PLot all q single graph
Live_q = (LiveFast_q+LiveSlow_q)/2;
Live_q = [Live_q mean(Live_q,2)];
Rec_q = (RecFast_q+RecSlow_q)/2;
Rec_q = [Rec_q mean(Rec_q,2)];
eLive = std(Live_q)/sqrt(size(Live_q,1));
eRec = std(Rec_q)/sqrt(size(Rec_q,1));

data = transpose(reshape([mean(Live_q, 1); mean(Rec_q,1)], [2, 9 ]));
errs = transpose(reshape([eLive; eRec], [2, 9 ]));
questions = {'Engagement','Enjoyment','Familiarity','Pleasure','Focus','Investment','Distraction','Spontaneity','Aggregate'};

figure;
colororder([.7 .1 .1; 0 .3 .7]) % first color is live, second color is recorded
% % Example data 
% model_series = [10 40 50 60; 20 50 60 70; 30 60 80 90]; 
% model_error = [1 4 8 6; 2 5 9 12; 3 6 10 13]; 
b = bar(data, 1);
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(data);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
    if i == nbars
        % b(i).FaceColor = [0 0 0]
    end
end
b(1).FaceColor = 'flat';
b(1).CData(9,:) = [.9 0 0];
b(2).FaceColor = 'flat';
b(2).CData(9,:) = [0 0 .9];
legend({'Live', 'Recorded'})
% Plot the errorbars
errorbar(x',data,errs, 'k','linestyle','none','LineWidth', 1.5, 'CapSize', 2, 'HandleVisibility','off');
xticklabels(questions)
xtickangle(-20)
yticks(1:5)
ylim([1 5])
ylabel("Rating")
title('Questionnaire Responses--Live vs. Recorded Excerpts');

hold off
set(gcf, 'color', 'white')
fontsize(20, 'points')

set(get(gca, 'XAxis'), 'FontWeight', 'bold');


 %% Corr between questions
responses_by_question = zeros(84,9); % 84 trials for 8 questions + 1 aggregate
for i = 1:8
    responses_by_question(:,i)=reshape(aggregate(:,:,i),84,1);
end
responses_by_question(:,9) = reshape(mean(aggregate,3),84,1);

[q_corr, q_pvals] = corrcoef(responses_by_question);

% replaces reduntant entries with nans
for i = 1:9
    for j = 1:9
        if i < j
            q_corr(i,j) = nan;
            q_pvals(i,j) = nan;
        end
    end
end

questions = {'Engagement','Enjoyment','Familiarity','Pleasure','Focus','Investment','Distraction (Reverse-Scored)','Spontaneity','Aggregate'};
figure;
t = tiledlayout('flow');
nexttile;
h=heatmap(questions,questions,q_corr,'MissingDataColor',[1 1 1],'ColorbarVisible','off','GridVisible','off');
h.CellLabelColor = 'none';
grid off;
h.Title = 'Question Correlations';
set(gcf,'color','white')
nexttile;
h = heatmap(questions,questions,q_pvals,'MissingDataColor',[1 1 1],'ColorbarVisible','off','GridVisible','off');
h.CellLabelColor = 'none';
% h.ColorLimits = [0, .2];
h.Title = 'Question Pvals';
colormap(h,bone)
set(gcf,'color','white')
fontsize(16, 'points')

%%
%loadSurveyResponsesParticipants: used for manually inputting responses from questionnaires
function loadSurveyResponsesParticipants(participant_table, study_path, startAt)
    for i = startAt: size(participant_table,1)
        participantID = char(participant_table{i,1}); 
        disp("Working on participant: " + participantID);
        for j = 1:4
            excerptNum = j;
            disp("---Excerpt " + num2str(excerptNum));
            % preprocessed_data_path = strcat(study_path,'/Data/Preprocessed_Data/',participantID, '/Excerpt_',num2str(excerptNum), '/');
            condition = char(participant_table{i,j+1});
            switch condition
                case 'LF'
                    condition = 'LiveFast';
                case 'LS'
                    condition = 'LiveSlow';
                case 'RF'
                    condition = 'RecordedFast';
                case 'RS'
                    condition = 'RecordedSlow';
            end
            save_file_path = strcat(study_path,'/Analysis/Participants/',participantID,'/',condition);
            if ~ isfolder(save_file_path)
                mkdir(save_file_path);
            end
            q_responses = zeros(8,1);
            for k = 1:8
                response = input(strcat("Response for question ", num2str(k), ": "));
                if isnumeric(response)
                    q_responses(k) = response;
                else
                    error("is not a number");
                end
            end
            save(strcat(save_file_path,'/survey'),"q_responses");
            writematrix(q_responses,strcat(save_file_path,'/survey.csv'));
           
        end
        disp("Finished participant: " + participantID);
    end
end


function aggregate = aggregateResponse(participant_table, study_path, startAt)
    study_path = strcat(study_path,'/Analysis');
    conditions = {'LiveFast', 'LiveSlow', 'RecordedFast', 'RecordedSlow'};
    totalnumParticipants = size(participant_table,1);
    aggregate = zeros(totalnumParticipants-startAt+1,4,8);
    for i = startAt:totalnumParticipants
        participantID = char(participant_table{i,1}); 
        disp("Working on participant: " + participantID);
        for j = 1:size(conditions,2)
            condition = conditions{j};
            disp (strcat('processing condition-', condition));
            load(strcat(study_path,'/Participants/', participantID, '/', condition,'/survey.mat'),"q_responses");
            aggregate(i,j,:) = q_responses;  
        end
        disp("Finished participant: " + participantID);
    end
end

function plotConditionAverages(LiveFast, LiveSlow,RecFast,RecSlow)
    eLiveFast = std(LiveFast)/sqrt(size(LiveFast,1));
    eLiveSlow = std(LiveSlow)/sqrt(size(LiveSlow,1));
    eRecFast = std(RecFast)/sqrt(size(RecFast,1));
    eRecSlow = std(RecSlow)/sqrt(size(RecSlow,1));
    x = ["LiveFast", "LiveSlow", "RecFast", "RecSlow"];
    data = [mean(LiveFast), mean(LiveSlow), mean(RecFast), mean(RecSlow)]';
    errhigh = [eLiveFast, eLiveSlow, eRecFast, eRecSlow];
    errlow = [eLiveFast, eLiveSlow, eRecFast, eRecSlow];
    
    bar(x,data)                
    ylim([1 5]);
    hold on
    
    er = errorbar([1 2 3 4],data,errlow, errhigh, 'LineWidth', 1.5);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    title('Mean Questionnaire Scores--All Conditions');
    hold off
end

function plotTwoCondition(Live, Rec)
    eLive = std(Live)/sqrt(size(Live,1));
    eRec = std(Rec)/sqrt(size(Rec,1));
    
    x = ["Live", "Rec"];
    data = [mean(Live), mean(Rec)]';
    errhigh = [eLive, eRec];
    errlow = errhigh;
    
    bar(x,data)                
    ylim([1 5]);
    hold on
    
    er = errorbar([1 2],data,errlow, errhigh, 'LineWidth', 1.5); %'Capsize', 20
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    title('Mean Questionnaire Scores--Live vs. Recorded');
    hold off
end

function plotAllTwoCondition(Live, Rec)
    t =tiledlayout('flow');
    t.Title.String = 'Reponses By Question--Live vs. Recorded';
    % title(t, "Ratings by Question Live vs. Recorded")
    nexttile;
    plotTwoCondition(Live(:,1), Rec(:,1));
    % title("How engaged were you during the musical excerpt?");
    title("Engagement");
    ylim([1 5]);
    yticks([1 5])
    
    nexttile;
    plotTwoCondition(Live(:,2), Rec(:,2));
    % title( "How much did you enjoy listening to the musical excerpt?");
    title("Enjoyment");
    ylim([1 5]);
    yticks([1 5]);
    
    nexttile;
    plotTwoCondition(Live(:,3), Rec(:,3));
    % title("How familiar was this excerpt?");
    title("Familiarity");
    ylim([1 5]);
    yticks([1 5])
    
    nexttile;
    plotTwoCondition(Live(:,4), Rec(:,4));
    % title("How pleasing was this excerpt?");
    title("Pleasure");
    ylim([1 5]);
    yticks([1 5])

    nexttile;
    plotTwoCondition(Live(:,5), Rec(:,5));
    % title("How often were you focused on the music?");
    title("Focus");
    ylim([1 5]);
    yticks([1 5])

    nexttile;
    plotTwoCondition(Live(:,6), Rec(:,6));
    % title("How emotionally invested were you while listening?");
    title("Investment");
    ylim([1 5]);
    yticks([1 5])

    nexttile;
    plotTwoCondition(Live(:,7), Rec(:,7));
    % title("How often were you distracted by musical mistakes? (Reverse Scored)");
    title("Distraction (Reversed)");
    ylim([1 5]);
    yticks([1 5])

    nexttile;
    plotTwoCondition(Live(:,8), Rec(:,8));
    % title("How spontaneous did this performance feel?");
    title("Spontaneity");
    ylim([1 5]);
    yticks([1 5])
end

function plotAllQuestions(LiveFast, LiveSlow,RecFast,RecSlow)
    tiledlayout('flow');

    nexttile;
    plotConditionAverages(LiveFast(:,1), LiveSlow(:,1),RecFast(:,1),RecSlow(:,1));
    % title("How engaged were you during the musical excerpt?");
    title("Engagement");
    ylim([1 5]);
    
    nexttile;
    plotConditionAverages(LiveFast(:,2), LiveSlow(:,2),RecFast(:,2),RecSlow(:,2));
    % title( "How much did you enjoy listening to the musical excerpt?");
    title("Enjoyment")
    ylim([1 5]);

    nexttile;
    plotConditionAverages(LiveFast(:,3), LiveSlow(:,3),RecFast(:,3),RecSlow(:,3));
    % title("How familiar was this excerpt?");
    title("Familiarity")
    ylim([1 5]);

    nexttile;
    plotConditionAverages(LiveFast(:,4), LiveSlow(:,4),RecFast(:,4),RecSlow(:,4));
    % title("How pleasing was this excerpt?");
    title("Pleasure");
    ylim([1 5]);

    nexttile;
    plotConditionAverages(LiveFast(:,5), LiveSlow(:,5),RecFast(:,5),RecSlow(:,5));
    % title("How often were you focused on the music?");
    title("Focus");
    ylim([1 5]);

    nexttile;
    plotConditionAverages(LiveFast(:,6), LiveSlow(:,6),RecFast(:,6),RecSlow(:,6));
    % title("How emotionally invested were you while listening?");
    title("Investment")
    ylim([1 5]);

    nexttile;
    plotConditionAverages(LiveFast(:,7), LiveSlow(:,7),RecFast(:,7),RecSlow(:,7));
    title("How often were you distracted by musical mistakes? (Reverse Scored)");
    title("Distraction (Reverse Scored)")
    ylim([1 5]);

    nexttile;
    plotConditionAverages(LiveFast(:,8), LiveSlow(:,8),RecFast(:,8),RecSlow(:,8));

    % title("How spontaneous did this performance feel?");
    title("Spontaneity")
    ylim([1 5]);

end

function allpower = avgFreqPower(study_path,participant_table, freqRange)
    study_path = strcat(study_path,'/Analysis');
    conditions = {'LiveFast', 'LiveSlow', 'RecordedFast', 'RecordedSlow'};
    allpower = zeros(size(participant_table,1), 4);
    for i = 1: size(participant_table,1)
        participantID = char(participant_table{i,1}); 
        disp(strcat("starting participant: ", participantID));
        for j = 1:size(conditions,2)
            condition = conditions{j};
            disp (strcat('processing condition-', condition));
            load(strcat(study_path,'/Participants/', participantID, '/', condition,'/Pwelch.mat'));
            freqs = matrix(:,1);
            avgpsd = matrix(:,2);
            % allpsds = matrix(:,3:end);
            fRange = freqs > freqRange(1) & freqs < freqRange(2);
            avgPower = mean(avgpsd.*fRange);
            allpower(i,j) = avgPower;
        end
    end
end