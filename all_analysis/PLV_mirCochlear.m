%% Gammatone PLV cochlear
clear variables;
study_path = '/Users/arun/Library/Mobile Documents/com~apple~CloudDocs/Documents/MINDLab/performance_study_final'; %% Only thing to update!
% addpath = strcat(study_path, '/Analysis/');
addpath(strcat(study_path, '/Analysis/EEGlab'))
participant_table = readtable("Participants.xlsx");
eeglab; close;
%%
numBins = 101; % number of frequencies (frequency resolution)
lowFreq = 0.2;
highFreq = 20.2;
fois = lowFreq:(highFreq-lowFreq)/(numBins-1):highFreq; % frequency values at which we calculate PLV
cycles = 5; % determines resolution of wavelet
aggregateChannelPLVs = zeros(size(participant_table,1),4,29,numBins);
aggregateChannelCoherences = zeros(size(participant_table,1),4,29,numBins);
%%
for i = 1: size(participant_table,1)
    for j = 1:4
        PLVs = zeros(numBins,1);
        Coherences = zeros(numBins,1);
        participantID = char(participant_table{i,1}); 
        excerptNum = j;
        disp("---Excerpt " + num2str(excerptNum));
        preprocessed_data_path = strcat(study_path,'/Data/Preprocessed_Data/',participantID, '/Excerpt_',num2str(excerptNum), '/');
        condition = char(participant_table{i,j+1});
        switch condition
            case 'LF'
                condition = 'LiveFast'; conditionNum = 1;
            case 'LS'
                condition = 'LiveSlow'; conditionNum = 2;
            case 'RF'
                condition = 'RecordedFast'; conditionNum = 3;
            case 'RS'
                condition = 'RecordedSlow'; conditionNum = 4;
        end
        analysis_file_path = strcat(study_path,'/Analysis/Participants/',participantID,'/',condition);
        
        % Load data files audio + EEG
        EEG = pop_loadset('filename','finished(ica).set','filepath', preprocessed_data_path);
        audio_signal = load(strcat(analysis_file_path,'/gammatone_cochlear_envelope.mat'),"filterbank_audio_data").filterbank_audio_data;
        audio_fs = 500;
        % Trim

        for k = 1:29 %iterate through every channel
            eegdata = EEG.data(k,:)';
            if (size(eegdata,1) > size(audio_signal,1))
                eegdata = eegdata(1:size(audio_signal));
            else
                audio_signal = audio_signal(1:size(eegdata));
            end

    
        % Calc PLV
    
            for l = 1:numBins
                cmw = designMorlet(fois(l), cycles, EEG.srate); % build a wavelet with N cycles
                [EEG_phase, EEG_amplitude, EEG_filtered, EEG_power] = waveletConvolution(eegdata,cmw); % convolve with EEG
                [Audio_phase, Audio_amplitude, Audio_filtered, Audio_power] = waveletConvolution(audio_signal,cmw); % convolve with Audio
                PLVs(l) = phaseCoherence(EEG_phase, Audio_phase);
                Coherences(l) = Coherence(Audio_phase, Audio_amplitude, EEG_phase, EEG_amplitude);
            end
            aggregateChannelPLVs(i,conditionNum,k,:) = PLVs;
            aggregateChannelCoherences(i,conditionNum,k,:) = Coherences;
        end
    end
end
%% Topo
topoSlice = squeeze(mean(aggregateChannelPLVs(1,:,:,5),[1 2]));
figure;
topoplot(topoSlice, EEG.chanlocs,  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo); 
%%
avgPLVs = squeeze(mean(aggregateChannelPLVs(:,:,:,40:48),[1 4]));
topoSlice_LiveFast = squeeze(avgPLVs(1,:,:));
topoSlice_LiveSlow = squeeze(avgPLVs(2,:,:));
topoSlice_RecFast = squeeze(avgPLVs(3,:,:));
topoSlice_RecSlow = squeeze(avgPLVs(4,:,:));
%% Plot raw topos
% eeglab redraw;
color_lim = [0 .10];
figure;
t = tiledlayout('flow');
t.Title.String = 'PLV Topos: 8-9.6 Hz';
t.Title.FontWeight = 'bold';
%frontal cluster channels: [1 2 3 4 6 7 8 9 12 24 25 26 27 28 29]

nexttile;
topoplot((topoSlice_LiveFast), EEG.chanlocs,'plotchans',[],'style','map','electrodes','ptslabels','chaninfo', EEG.chaninfo);  %, 'whitebk', 'on'
title('LiveFast');
clim(color_lim);
nexttile;
topoplot((topoSlice_LiveSlow), EEG.chanlocs,  'plotchans',[],'style','map','electrodes','ptslabels', 'chaninfo', EEG.chaninfo); 
title('LiveSlow');
clim(color_lim);
nexttile;
topoplot((topoSlice_RecFast), EEG.chanlocs, 'plotchans',[],'style','map','electrodes','ptslabels', 'chaninfo', EEG.chaninfo); 
title('RecordedFast');
clim(color_lim);
nexttile;
topoplot((topoSlice_RecSlow), EEG.chanlocs,  'plotchans',[],'style','map','electrodes','ptslabels', 'chaninfo', EEG.chaninfo); 
title('RecordedSlow');
clim(color_lim);

cb = colorbar;
cb.Layout.Tile = "east";
set(gcf,"color","white")
%% Plot raw topos with only some electrodes labeled
% eeglab redraw;
electrodes_to_keep = ["Fz", "F3", "F4"];
electrodes_to_keep_nums = [4, 26, 8];

EEG_copy = EEG;
for chan_num=1:29
    matchFound = false;
    for electrode=1:length(electrodes_to_keep)
        if convertCharsToStrings(EEG_copy.chanlocs(chan_num).labels) == convertCharsToStrings(electrodes_to_keep(electrode))
            matchFound = true;
        end
    end
    if ~matchFound
        EEG_copy.chanlocs(chan_num).labels = ' ';
    end
end

color_lim = [0 .10];
figure;
t = tiledlayout('flow');
t.Title.String = 'PLV Topos: 8-9.6 Hz';
t.Title.FontWeight = 'bold';
%frontal cluster channels: [1 2 3 4 6 7 8 9 12 24 25 26 27 28 29]

nexttile;
topoplot((topoSlice_LiveFast), EEG_copy.chanlocs,'plotchans',[],'style','both','electrodes','ptslabels','emarker2',{electrodes_to_keep_nums,'.','k', 11},'emarker', {'.','k',[],1}, 'chaninfo',  EEG.chaninfo);  %, 'whitebk', 'on'
title('LiveFast');
clim(color_lim);
nexttile;
topoplot((topoSlice_LiveSlow), EEG_copy.chanlocs,  'plotchans',[],'style','both','electrodes','ptslabels', 'emarker2',{electrodes_to_keep_nums,'.','k', 11},'emarker', {'.','k',[],1}, 'chaninfo',  EEG.chaninfo);
title('LiveSlow');
clim(color_lim);
nexttile;
topoplot((topoSlice_RecFast), EEG_copy.chanlocs, 'plotchans',[],'style','both','electrodes','ptslabels', 'emarker2',{electrodes_to_keep_nums,'.','k', 11},'emarker', {'.','k',[],1}, 'chaninfo',  EEG.chaninfo);
title('RecordedFast');
clim(color_lim);
nexttile;
topoplot((topoSlice_RecSlow), EEG_copy.chanlocs,  'plotchans',[],'style','both','electrodes','ptslabels','emarker2',{electrodes_to_keep_nums,'.','k', 11},'emarker', {'.','k',[],1}, 'chaninfo',  EEG.chaninfo);
title('RecordedSlow');
clim(color_lim);

fontsize(20, "points")

cb = colorbar;
cb.Layout.Tile = "east";
set(gcf,"color","white")
%% Plot raw topos Live v Rec
% eeglab redraw;
color_lim = [0 .1];
figure;
t = tiledlayout('flow');
t.Title.String = 'Phase Locking at Frequencies of Interest: 3-3.4 Hz';
t.Title.FontWeight = 'bold';
%frontal cluster channels: [1 2 3 4 6 7 8 9 12 24 25 26 27 28 29]

nexttile;
topoplot((topoSlice_LiveFast+topoSlice_LiveSlow)/2, EEG.chanlocs,'plotchans',[],'style','map','electrodes','ptslabels','chaninfo', EEG.chaninfo,'whitebk', 'on');  %, 'whitebk', 'on'
title('Live');
clim(color_lim);
nexttile;
topoplot((topoSlice_RecFast+topoSlice_RecSlow)/2, EEG.chanlocs, 'plotchans',[],'style','map','electrodes','ptslabels', 'chaninfo', EEG.chaninfo,'whitebk', 'on'); 
title('Recorded');
clim(color_lim);
set(gcf,"color","white")
cb = colorbar;
cb.Layout.Tile = "east";
%% Plot raw topos Fast v Slow
% eeglab redraw;
color_lim = [0 .1];
figure;
t = tiledlayout('flow');
t.Title.String = 'Fast Vs. Slow PLV at 3.2 Hz';
t.Title.FontWeight = 'bold';
%frontal cluster channels: [1 2 3 4 6 7 8 9 12 24 25 26 27 28 29]

nexttile;
topoplot((topoSlice_LiveFast+topoSlice_RecFast)/2, EEG.chanlocs,'plotchans',[],'style','map','electrodes','ptslabels','chaninfo', EEG.chaninfo);  %, 'whitebk', 'on'
title('Fast');
clim(color_lim);
nexttile;
topoplot((topoSlice_LiveSlow+topoSlice_RecSlow)/2, EEG.chanlocs, 'plotchans',[],'style','map','electrodes','ptslabels', 'chaninfo', EEG.chaninfo); 
title('Slow');
clim(color_lim);
set(gcf,"color","white")
cb = colorbar;
cb.Layout.Tile = "east";
%% Plot normalized topos
% eeglab redraw;
figure;
t = tiledlayout('flow');
t.Title.String = 'Phase locking NORMALIZED Topos from .2 - 10.2 Hz';
t.Title.FontWeight = 'bold';

nexttile;
topoplot(normalize(topoSlice_LiveFast), EEG.chanlocs,  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo); 
title('LiveFast');
clim([-1 1]);

nexttile;
topoplot(normalize(topoSlice_LiveSlow), EEG.chanlocs,  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo); 
title('LiveSlow');
clim([-1 1]);

nexttile;
topoplot(normalize(topoSlice_RecFast), EEG.chanlocs,  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo); 
title('RecordedFast');
clim([-1 1]);

nexttile;
topoplot(normalize(topoSlice_RecSlow), EEG.chanlocs,  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo); 
title('RecordedSlow');
clim([-1 1]);

cb = colorbar;
cb.Layout.Tile = "east";

%% Plot live rec at Fz
% avg [1 2 3 4 6 7 8 9 12 24 25 26 27 28 29]
% 4 and 12

aggregatePLVs = (squeeze(mean(aggregateChannelPLVs(:,:,[26 4 8],:),3)));
% aggregatePLVs = (squeeze(mean(aggregateChannelPLVs(:,:,1:9,:),3)) + squeeze((aggregateChannelPLVs(:,:,12,:))) + squeeze(mean(aggregateChannelPLVs(:,:,24:29,:),3)))/3;
% aggregatePLVs = (squeeze(mean(aggregateChannelPLVs(:,:,1:9,:),3)) + squeeze((aggregateChannelPLVs(:,:,12,:))) + squeeze(mean(aggregateChannelPLVs(:,:,24:29,:),3)))/3;
PLV_Live = squeeze(mean(aggregatePLVs(:,1:2,:),[1 2]));
PLV_Rec = squeeze(mean(aggregatePLVs(:,3:4,:),[1 2]));

figure;
hold on;
plot(fois, PLV_Live, 'LineWidth',2, 'Color',[1 0 0]);
plot(fois, PLV_Rec, 'LineWidth',2, 'Color',[0 0 1]);
legend('live', 'rec')
hold off;
%% Plot Fast Slow at Fz -> refer to other script. it should just run.
%% Plot by diff Q score Live/Rec
Q_data = squeeze(mean(aggregate(:,conditions,:),3));
participantDiffscores = [1:size(Q_data,1); zeros(1,size(Q_data,1))];
for i = 1:size(Q_data, 1)
    participantDiffscores(2,i) = mean(Q_data(i,1:2), 2) - mean(Q_data(i,3:4),2);   %Get live minus recorded
end
participant_ranking = sortrows(participantDiffscores',2)
figure;
bar(participant_ranking(:,2))
title('Difference In Ratings Between Live and Recorded')
xticks(1:21)
xticklabels(participant_ranking(:,1))
xlabel('Participant #')
ylabel('Live-Rec Aggregate Rating Difference')
fontsize(20,'points')
set(gcf, 'color', 'white')
%% continue Plot by diff Q score All Condition
electrodes = [26 4 8];
figure;
t = tiledlayout(2,1);
t.Title.String = 'PLVs Grouped by Questionnaire Difference (Live-Rec)';
% hold on


for condition = [1 3 2 4]
    switch condition
        case 1
            % condition = [1 2]
            nexttile; hold on;
            colorHigh = [0.8, .2, 0];
            colorLow = [0.6, .2, 0];
            label = "Live,";
            title("Fast")
            xlim([0.2,16.2])
            ylim([0, .14])
            yticks(0:.02:.14)
            % xticks(sort([xticks 8.4])); xline([8.4],"HandleVisibility","off","LineWidth",1);
            xlabel("Frequency [Hz]"); ylabel("PLV");
        case 3
            colorHigh = [0, .2, 0.8];
            colorLow = [0, .2, 0.4];
            label = "Recorded,";
        case 2
            % condition = [3 4]
            hold off; nexttile; hold on;
            colorHigh = [0.8, .2, 0];
            colorLow = [0.6, .2, 0];
            label = "Live,";
            title("Slow")
            xlim([0.2,6.2])
            ylim([0, .14])
            yticks(0:.02:.14)
            % xticks(sort([xticks .8])); xline([.8],"HandleVisibility","off","LineWidth",1);
            xlabel("Frequency [Hz]"); ylabel("PLV");
    
        case 4
            colorHigh = [0, .2, 0.8];
            colorLow = [0, .2, 0.4];
            label = "Recorded,";
    end
    disp(condition)
    lowPLV = squeeze(mean(aggregateChannelPLVs(participant_ranking(1:10,1),condition,electrodes,:),[2 3]));
    highPLV = squeeze(mean(aggregateChannelPLVs(participant_ranking(12:21,1),condition,electrodes,:),[2 3]));
    
    errlow = std(lowPLV,0,1)/sqrt(size(lowPLV,1));
    errhigh = std(highPLV,0,1)/sqrt(size(highPLV,1));
    
    hp = patch([8 8 9.6 9.6],[0 .14 .14 0],'k',...
    'facecolor',[.5 .5 .5],'facealpha', 0.2,'edgecolor','none', 'HandleVisibility', 'off') ;
    hp2 = patch([0.6 0.6 1 1],[0 .14 .14 0],'k',...
    'facecolor',[.5 .5 .5],'facealpha', 0.2,'edgecolor','none', 'HandleVisibility', 'off') ;
    xline(8,'-',{'8 Hz'}, 'LineWidth',1, 'HandleVisibility', 'off', 'color', [.8, 0, .2], 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center','LabelOrientation', 'align');
    xline(9.6,'-',{'9.6 Hz'}, 'LineWidth',1, 'HandleVisibility', 'off', 'color', [.8, 0, .2], 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center','LabelOrientation', 'align');
    xline(0.6,'-',{'0.6 Hz'}, 'LineWidth',1, 'HandleVisibility', 'off', 'color', [.8, 0, .2], 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center','LabelOrientation', 'align');
    xline(1,'-',{'1 Hz'}, 'LineWidth',1,'HandleVisibility', 'off', 'color', [.8, 0, .2], 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center','LabelOrientation', 'align');

    shadedErrorBar(fois,mean(highPLV),errhigh,'lineprops',{'Color',colorHigh,'LineWidth',2,'DisplayName',label + " High Rating Diff"})
    shadedErrorBar(fois,mean(lowPLV),errlow,'lineprops',{'Color',colorLow,'LineWidth',2,'DisplayName',label + " Low Rating Diff"})
    legend;
end


hold off;
xlabel("Frequency [Hz]")
ylabel("PLV")

% legend;
fontsize(20,'points')
set(gcf,'color','white')
%% test maxes within participant %% Check out participant 16??
figure
plot(fois,squeeze(mean(aggregateChannelPLVs(18,1,electrodes,:),3)))

frequency_indeces =  fois >= 1.8 & fois<=2.2; %fois >= 8 & fois<=9.6; %

max_PLVs = zeros(21, 4);


for i = 1:size(max_PLVs, 1)
    max_PLVs(i,:) = max(squeeze(mean(aggregateChannelPLVs(i,:,electrodes,frequency_indeces), 3)), [], 2);
end


% 
% livefast = max_PLVs(participant_ranking(12:21,1),1);
% recfast = max_PLVs(participant_ranking(12:21,1),3); % change for conditions

livefast = max_PLVs(:,1);
recfast = max_PLVs(:,3);

[h, p, ci, stats] = ttest(livefast, recfast)

disp(num2str(p))

figure
hold on
xlim([.5, 2.5])
ylim([0, .37])
title("max plv in range")
pl = [recfast; livefast];
pl_f = [repmat({'RecFast'},length(recfast),1); repmat({'LiveFast'},length(livefast),1)];
boxplot(pl,pl_f)
hold off
%%  LME Model

% PLV ~ Liveness*Tempo*Q_Diff + 1|Participant + 1

frequency_indeces =  fois >= 8 & fois<=9.6; %fois >= 0.6 & fois<=1; 
pruned_aggregatePLV = squeeze(mean(aggregatePLVs(:,:,frequency_indeces), 3));

pruned_aggregatePLV = max_PLVs;

% pruned_aggregatePLV = squeeze(mean(zScorePLV(:,:,frequency_indeces), 3));


% 1 is live, 2 is recorded
performanceMediums = zeros(84, 1); 
% 1 is fast, 2 is slow
tempi = zeros(84, 1); 

qdiffs = zeros(84, 1); 

participant_Nums = zeros(84, 1);
% participant_Nums = cell(84);

q_datas = zeros(84, 1);
% q_hilo = zeros(84, 1);
all_PLVs = zeros(84, 1); 

diff_category = zeros(84,1);

for participant=participantDiffscores(1,:)
    qdiffs(participant*4-3:participant*4) = participantDiffscores(2,participant);
    if participantDiffscores(2,participant) > median(participantDiffscores(2,:))
        diff_category(participant*4-3:participant*4) = 0;
    elseif participantDiffscores(2,participant) < median(participantDiffscores(2,:))
        diff_category(participant*4-3:participant*4) = 1;
    end
    participant_Nums(participant*4-3:participant*4) = participant;
    performanceMediums(participant*4-3:participant*4-2) = 0; performanceMediums(participant*4-1:participant*4) = 1;
    tempi([participant*4-3, participant*4-1]) = 0; tempi([participant*4-2, participant*4]) = 1;
    all_PLVs(participant*4-3:participant*4) = squeeze(pruned_aggregatePLV(participant,:));
    q_datas(participant*4-3:participant*4) = squeeze(Q_data(participant,:));
end
participant_Nums = categorical(participant_Nums,1:21,string(transpose(1:21)));
performanceMediums = categorical(performanceMediums,[0 1],{'Live' 'Recorded'});
tempi = categorical(tempi,[0 1],{'Fast' 'Slow'});
diff_category = categorical(diff_category, [0 1], {'High', 'Low'});
tbl = table(all_PLVs, performanceMediums, tempi, qdiffs, diff_category, q_datas, participant_Nums,'VariableNames',{'PLV','Liveness','Tempo', 'Q_diff', 'Diff_Split','qData','ParticipantNum'});
% lme_function = 'PLV~1+Liveness*Tempo + (1|Q_diff:ParticipantNum)';

% lme_function = 'PLV~1+Q_diff *Tempo*Liveness + (1|ParticipantNum)';
% lme_function = 'PLV~1+ Tempo*Liveness + (1|ParticipantNum)'; 
% lme_function = 'PLV~1+Diff_Split *Tempo*Liveness + (1|ParticipantNum)'; % test questionnaire effect on low peak

% lme_function = 'PLV~1+Tempo*Liveness + (1|Q_diff:ParticipantNum) + (1|Q_diff) + (1|ParticipantNum)';
lme_function = 'PLV~1+qData*Tempo*Liveness + (1|ParticipantNum)';
% lme_function = 'qData~1+Liveness*Tempo + (1|ParticipantNum)';
% lme_function = 'qData~1+Liveness*Tempo*PLV +(1|ParticipantNum)';
% lme_function = ' PLV~1+Liveness*Tempo +(';
lme = fitlme(tbl,lme_function, 'FitMethod', 'REML');

lme;

%%  LME Model for just Fast conditions

% PLV ~ Liveness*Tempo*Q_Diff + 1|Participant + 1

frequency_indeces = fois >= 8 & fois<=9.6; %fois >= 1.8 & fois<=2.4;  %fois >= 0.6 & fois<=1; %

pruned_aggregatePLV = squeeze(mean(aggregatePLVs(:,[1, 3],frequency_indeces), 3));

% pruned_aggregatePLV = squeeze(mean(zScorePLV(:,[1, 3],frequency_indeces), 3));


% 1 is live, 2 is recorded
performanceMediums = zeros(42, 1); 

qdiffs = zeros(42, 1); 

participant_Nums = zeros(42, 1);
% participant_Nums = cell(84);

q_datas = zeros(42, 1);
% q_hilo = zeros(84, 1);
all_PLVs = zeros(42, 1); 

% diff_category = zeros(84,1);

for participant=participantDiffscores(1,:)
    qdiffs(participant*2-1:participant*2) = participantDiffscores(2,participant);
    % if participantDiffscores(2,participant) > median(participantDiffscores(2,:))
    %     diff_category(participant*4-3:participant*4) = 0;
    % elseif participantDiffscores(2,participant) < median(participantDiffscores(2,:))
    %     diff_category(participant*4-3:participant*4) = 1;
    % end
    participant_Nums(participant*2-1:participant*2) = participant;
    performanceMediums(participant*2-1) = 0; performanceMediums(participant*2) = 1;
   
    all_PLVs(participant*2-1:participant*2) = squeeze(pruned_aggregatePLV(participant,:));
    q_datas(participant*2-1:participant*2) = squeeze(Q_data(participant,[1 3]));
end
participant_Nums = categorical(participant_Nums,1:21,string(transpose(1:21)));
performanceMediums = categorical(performanceMediums,[0 1],{'Live' 'Recorded'});
% tempi = categorical(tempi,[0 1],{'Fast' 'Slow'});
% diff_category = categorical(diff_category, [0 1], {'High', 'Low'});
% tbl = table(all_PLVs, performanceMediums, tempi, qdiffs, diff_category, q_datas, participant_Nums,'VariableNames',{'PLV','Liveness','Tempo', 'Q_diff', 'Diff_Split','qData','ParticipantNum'});
tbl = table(all_PLVs, performanceMediums, qdiffs, q_datas, participant_Nums,'VariableNames',{'PLV','Liveness', 'Q_diff','qData','ParticipantNum'});

% lme_function = 'PLV~1+Liveness*Tempo + (1|Q_diff:ParticipantNum)';

lme_function = 'PLV~1+Q_diff *Liveness + (1|ParticipantNum)';
% lme_function = 'PLV~1+qData *Liveness + (1|ParticipantNum)';


% lme_function = 'PLV~1+Liveness + (1|Q_diff:ParticipantNum) + (1|Q_diff) + (1|ParticipantNum)';

lme = fitlme(tbl,lme_function, 'FitMethod', 'REML');
%% T-test high/low participant diff score live fast

frequency_indeces =   fois >= 8 & fois<=9.6; % fois >= 0.6 & fois<=1;fois == 2; %fois ==0.8; %fois >= 0.6 & fois<=1;%
pruned_aggregatePLV = squeeze(mean(aggregatePLVs(:,:,frequency_indeces), 3));

% low_participants = participant_ranking(1:10,1);
% low_participants = participant_ranking(1:8,1); % for median split
high_participants = participant_ranking(12:21,1);
% high_participants = participant_ranking(12:21,1); % for median split

liveFasts = pruned_aggregatePLV(high_participants, 1);  
recFasts = pruned_aggregatePLV(high_participants, 3); 
% liveFasts = pruned_aggregatePLV(:, 1);  
% recFasts = pruned_aggregatePLV(:, 3); 

[h, p, ci, stats] = ttest(liveFasts, recFasts, 'Tail', 'Right');
%%

diffScores = sort(participant_ranking, 1)
figure
scatter(diffScores(:, 2), pruned_aggregatePLV(:, 1)- pruned_aggregatePLV(:, 3))

y = pruned_aggregatePLV(:, 1)- pruned_aggregatePLV(:, 3)
p = polyfit(diffScores(:, 2),y,1)
yfit = polyval(p,diffScores(:, 2));

yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal
rsq_adj = 1 - SSresid/SStotal * (length(y)-1)/(length(y)-length(p))



%%
[h, p] = ttest(liveFasts, mean(pruned_aggregatePLV(:,1),1) ,'Tail','Right');
% %% Anova for high rated people
% frequency_indeces = fois >= 8 & fois<=9.6; % fois >= 0.6 & fois<=1;
% avgLiveFast = squeeze(mean(aggregatePLVs(high_participants,1,frequency_indeces), 3));
% avgLiveSlow = squeeze(mean(aggregatePLVs(high_participants,2,frequency_indeces), 3));
% avgRecFast = squeeze(mean(aggregatePLVs(high_participants,3,frequency_indeces), 3));
% avgRecSlow = squeeze(mean(aggregatePLVs(high_participants,4,frequency_indeces), 3));
% anova_size = length(high_participants);
% A = zeros(anova_size*2,2); %% row factor is tempo f/s, column factor is live/rec
% A(1:anova_size,1) = avgLiveFast;
% A(1:anova_size,2) = avgRecFast;
% A((anova_size+1):(anova_size*2),1) = avgLiveSlow;
% A((anova_size+1):(anova_size*2),2) = avgRecSlow;
% [~,~,stats] = anova2(A,anova_size); figure;
% c1 = multcompare(stats); figure;
% c2 = multcompare(stats,"Estimate","row");


%% PLOT PLV BY Q Score

questions = {'Engagement','Enjoyment','Familiarity','Pleasure','Focus','Investment','Distraction (Reverse-Scored)','Spontaneity'};
conditions = 1:4;
electrodes = [4, 8, 26];
figure;
% t = tiledlayout('flow');
% t.Title.String = "PLV by Questionnaire High/Low Rating";
% t.Title.FontWeight = 'Bold';


Q_data = squeeze(mean(aggregate(:,conditions,:),3));
q_low = find(Q_data<median(Q_data(:, 4),'all'));
q_high = find(Q_data>median(Q_data(:, 4),'all'));
n_low = length(q_low);
n_high = length(q_high);
[rlow, clow] = ind2sub([21,4],q_low);
[rhigh, chigh] = ind2sub([21,4],q_high);
low_PLV = zeros(length(q_low),size(aggregateChannelPLVs,4));
high_PLV = zeros(length(q_high),size(aggregateChannelPLVs,4));
numLow = 0;
numHigh = 0;
for i = 1:length(q_low)
    if (clow(i) == 4 ) % check condition is live or rec
        numLow = numLow + 1;
        low_PLV(i,:) = squeeze(mean(aggregateChannelPLVs(rlow(i),clow(i),electrodes,:),3));
    end
end
for i = 1:length(q_high)
    if (chigh(i) == 4)
        numHigh = numHigh + 1;
        high_PLV(i,:) = squeeze(mean(aggregateChannelPLVs(rhigh(i),chigh(i),electrodes,:),3));
    end
end
% test_PLV = low_PLV(all(low_PLV,2),:);
% high_TRF_flux = high_TRF_flux(all(high_TRF_flux,2),:);
% 
errlow = std(low_PLV,0,1)/sqrt(length(low_PLV));
errhigh = std(high_PLV,0,1)/sqrt(length(high_PLV));

hold on
shadedErrorBar(fois,mean(high_PLV),errhigh,'lineprops',{'Color',[0.8, .2, 0],'LineWidth',2,'DisplayName',"High Rating n=" + num2str(numHigh)})
shadedErrorBar(fois,mean(low_PLV),errlow,'lineprops',{'Color',[0, .2, 0.8],'LineWidth',2,'DisplayName',"Low Rating n=" + num2str(numLow)})
hold off;
title("Aggregate")
xlabel("Frequency [Hz]")
ylabel("PLV")
ylim([0, .1])
legend;
fontsize(15,'points')
set(gcf,'color','white')





%% Plot all at Fz

 aggregatePLVs = (squeeze(mean(aggregateChannelPLVs(:,:,:,:),3)));
 % aggregatePLVs = (squeeze(aggregateChannelPLVs(:,:,2,:)) + squeeze(aggregateChannelPLVs(:,:,4,:)) + squeeze(aggregateChannelPLVs(:,:,8,:)) + squeeze(aggregateChannelPLVs(:,:,26,:)) + squeeze(aggregateChannelPLVs(:,:,12,:)))/5;
% aggregatePLVs = (squeeze(mean(aggregateChannelPLVs(:,:,1:9,:),3)) + squeeze((aggregateChannelPLVs(:,:,12,:))) + squeeze(mean(aggregateChannelPLVs(:,:,24:29,:),3)))/3;
PLV_LiveFast = squeeze(mean(aggregatePLVs(:,1,:),[1 2]));
PLV_LiveSlow = squeeze(mean(aggregatePLVs(:,2,:),[1 2]));
PLV_RecFast = squeeze(mean(aggregatePLVs(:,3,:),[1 2]));
PLV_RecSlow = squeeze(mean(aggregatePLVs(:,4,:),[1 2]));

figure;
hold on;
plot(fois, PLV_LiveFast, 'LineWidth',2, 'Color',[1 0 0],'DisplayName','LiveFast');
plot(fois, PLV_LiveSlow, 'LineWidth',2, 'Color',[1 .5 0],'DisplayName','LiveSlow');
plot(fois, PLV_RecFast, 'LineWidth',2, 'Color',[0 0 1],'DisplayName','RecFast');
plot(fois, PLV_RecSlow, 'LineWidth',2, 'Color',[0 1 0],'DisplayName','RecSlow');
hold off;
legend;
%% Anova at 4.2 Hz for Fz
avgLiveFast = squeeze(aggregateChannelPLVs(:,1,4,20));
avgLiveSlow = squeeze(aggregateChannelPLVs(:,2,4,20));
avgRecFast = squeeze(aggregateChannelPLVs(:,3,4,20));
avgRecSlow = squeeze(aggregateChannelPLVs(:,4,4,20));
A = zeros(42,2); %% row factor is tempo f/s, column factor is live/rec
A(1:21,1) = avgLiveFast;
A(1:21,2) = avgRecFast;
A(22:42,1) = avgLiveSlow;
A(22:42,2) = avgRecSlow;
[~,~,stats] = anova2(A,21); figure;
c1 = multcompare(stats); figure;
c2 = multcompare(stats,"Estimate","row");
%% Anova at 4.2-6.2 Hz for Fz
avgLiveFast = squeeze(mean(aggregatePLVs(:,1,3:6),3));
avgLiveSlow = squeeze(mean(aggregatePLVs(:,2,3:6),3));
avgRecFast = squeeze(mean(aggregatePLVs(:,3,3:6),3));
avgRecSlow = squeeze(mean(aggregatePLVs(:,4,3:6),3));
A = zeros(42,2); %% row factor is tempo f/s, column factor is live/rec
A(1:21,1) = avgLiveFast;
A(1:21,2) = avgRecFast;
A(22:42,1) = avgLiveSlow;
A(22:42,2) = avgRecSlow;
[~,~,stats] = anova2(A,21); figure;
c1 = multcompare(stats); figure;
c2 = multcompare(stats,"Estimate","row");

 %% Plot significance at Frequencies not working
% sigsLiveRec = zeros(length(fois),1);
% sigsFastSlow = zeros(length(fois),1);
% for i = 1:length(fois)
%     avgLiveFast = squeeze(aggregatePLVs(:,1,i));
%     avgLiveSlow = squeeze(aggregatePLVs(:,2,i));
%     avgRecFast = squeeze(aggregatePLVs(:,3,i));
%     avgRecSlow = squeeze(aggregatePLVs(:,4,i));
%     A = zeros(42,2); %% row factor is tempo f/s, column factor is live/rec
%     A(1:21,1) = avgLiveFast;
%     A(1:21,2) = avgRecFast;
%     A(22:42,1) = avgLiveSlow;
%     A(22:42,2) = avgRecSlow;
%     [~,~,stats] = anova2(A,21); %figure;
%     c1 = multcompare(stats); %figure;
%     c2 = multcompare(stats,"Estimate","row");
%     sigsLiveRec(i) = c1(6);
%     sigsFastSlow(i) = c1(3);
% end
% 

%tempo c1(3) live c1(6)
%% barGraph PLV Fz 4.2 all cond
    eLiveFast = std(avgLiveFast)/sqrt(size(avgLiveFast,1));
    eLiveSlow = std(avgLiveSlow)/sqrt(size(avgLiveSlow,1));
    eRecFast = std(avgRecFast)/sqrt(size(avgRecFast,1));
    eRecSlow = std(avgRecSlow)/sqrt(size(avgRecSlow,1));
    x = ["LiveFast", "LiveSlow", "RecFast", "RecSlow"];
    data = [mean(avgLiveFast), mean(avgLiveSlow), mean(avgRecFast), mean(avgRecSlow)]';
    errhigh = [eLiveFast, eLiveSlow, eRecFast, eRecSlow];
    errlow = [eLiveFast, eLiveSlow, eRecFast, eRecSlow];
    figure;
    bar(x,data)                
    ylim([0 .2]);
    hold on
    
    er = errorbar([1 2 3 4],data,errlow, errhigh);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    title('PLV at Fz at 4.2 Hz');
    set(gcf,'color','w');
    hold off
%% barGraph PLV Fz 4.2 Live v Rec
    avgLive = (avgLiveFast + avgLiveSlow)/2;
    avgSlow = (avgRecFast + avgRecSlow)/2;
    eLive = std(avgLive)/sqrt(size(avgLive,1));
    eSlow = std(avgSlow)/sqrt(size(avgSlow,1));
   
    x = ["Live", "Recorded"];
    data = [mean(avgLive), mean(avgSlow)]';
    errhigh = [eLive, eSlow];
    errlow = errhigh;
    figure;
    bar(x,data)                
    ylim([0 .2]);
    hold on
    
    er = errorbar([1 2],data,errlow, errhigh);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    title('PLV at Fz at 4.2 Hz Live Vs. Rec');
    set(gcf,'color','w');
    hold off
%% All Part All Condition PLV Fz 4.2 Plot
all_table_fz = [avgLiveFast avgLiveSlow avgRecFast avgRecSlow]';
figure; plot(all_table_fz)
ylabel("PLV"); xlabel("Condition");
xticks([1 2 3 4])
xticklabels({'LiveFast', 'LiveSlow', 'RecFast', 'RecSlow'})
title("All Participants Live/Slow PLV Fz 4.2 Hz")
set(gcf,'color','w');
%% All Part Live Rec PLV Fz 4.2 Plot
all_table_fz = [(avgLiveFast+avgLiveSlow)/2 (avgRecFast+avgRecSlow)/2]';
figure; plot(all_table_fz)
ylabel("PLV"); xlabel("Condition");
xticks([1 2])
xticklabels({'Live', 'Rec'})
title("All Participants PLV Fz 4.2 Hz")
set(gcf,'color','w');
%% All Part Fast Slow PLV Fz 4.2 Plot
all_table_fz = [(avgLiveFast+avgRecFast)/2 (avgLiveSlow+avgRecSlow)/2]';
figure; plot(all_table_fz)
ylabel("PLV"); xlabel("Condition");
xticks([1 2])
xticklabels({'Fast', 'Slow'})
title("All Participants PLV Fz 4.2 Hz")
set(gcf,'color','w');
%%
figure;
hold on;
plot(fois, squeeze(mean(aggregatePLVs(:,1:2,:),2)), 'LineWidth',2, 'Color',[1 0 0 .5]);
plot(fois, squeeze(mean(aggregatePLVs(:,3:4,:),2)), 'LineWidth',2, 'Color',[0 0 1 .5]);
legend('live', 'rec')
hold off;
% [up, lo] = envelope(audioRS,600);
% env = (up+lo)/2;

%%
figure;
hold on;
plot(fois, squeeze(mean(aggregatePLVs(:,1,:),2)), 'LineWidth',1, 'Color',[1 0 0 .5]);
plot(fois, squeeze(mean(aggregatePLVs(:,3,:),2)), 'LineWidth',1, 'Color',[0 0 1 .5]);
% legend('liveFast', 'recFast')

% [up, lo] = envelope(audioRS,600);
% env = (up+lo)/2;

plot(fois, squeeze(mean(aggregatePLVs(:,2,:),2)), 'LineWidth',1, 'Color',[1 .5 0 .5]);
plot(fois, squeeze(mean(aggregatePLVs(:,4,:),2)), 'LineWidth',1, 'Color',[0 1 0 .5]);
% legend('liveSlow', 'recSlow')
hold off;
% [up, lo] = envelope(audioRS,600);
% env = (up+lo)/2;
%%
LiveFast = squeeze(aggregatePLVs(:,1,2)); % at 2nd foi
LiveSlow = squeeze(aggregatePLVs(:,2,2));
RecFast = squeeze(aggregatePLVs(:,3,2));
RecSlow = squeeze(aggregatePLVs(:,4,2));
    eLiveFast = std(LiveFast)/sqrt(size(LiveFast,1));
    eLiveSlow = std(LiveSlow)/sqrt(size(LiveSlow,1));
    eRecFast = std(RecFast)/sqrt(size(RecFast,1));
    eRecSlow = std(RecSlow)/sqrt(size(RecSlow,1));
    x = ["LiveFast", "LiveSlow", "RecFast", "RecSlow"];
    data = [mean(LiveFast), mean(LiveSlow), mean(RecFast), mean(RecSlow)];
    errhigh = [eLiveFast, eLiveSlow, eRecFast, eRecSlow];
    errlow = [eLiveFast, eLiveSlow, eRecFast, eRecSlow];
    figure;
    bar(x,data)                
    ylim([0 .2]);
    hold on
    
    er = errorbar([1 2 3 4],data,errlow, errhigh);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    title('Questionairre Scores for All Conditions');
    hold off
%%
function PLV = phaseCoherence(angles1, angles2) %%phase coherence for two angle time series
    PLV = abs(sum(exp(1i*(angles1-angles2))))/length(angles1);
end

function coherence = Coherence(angles1, amplitudes1, angles2, amplitudes2) %% Wollman Neural Entrainment melodic spectral complexity
    numerator = abs(sum(amplitudes1.*amplitudes2.*exp(1i*(angles1-angles2))));
    denomenator = sqrt(sum(amplitudes1.^2)*sum(amplitudes2.^2));
    coherence = numerator/denomenator;
end

% reference with sample code: https://mikexcohen.com/lectures.html 
function [cmw] = designMorlet(waveletFreq, waveletCycles, srate) %% more cycles is better frequency precision but worse time precision
    time  = -2:1/srate:2; % best practice is to have time=0 at the center of the wavelet
    % create complex sine wave
    sine_wave = exp( 1i*2*pi*waveletFreq.*time );
    % create Gaussian window
    s = waveletCycles / (2*pi*waveletFreq); % this is the standard deviation of the gaussian
    gaus_win  = exp( (-time.^2) ./ (2*s^2) );
    % now create Morlet wavelet
    cmw = gaus_win.*sine_wave;
end

function [phase_timeSeries, amplitude_timeSeries, filtered_timeSeries, power_timeSeries] = waveletConvolution(data, cmw)
    nData = length(data);
    nKern = length(cmw);
    nConv = nData + nKern - 1;
    %FFTs:
    % note that the "N" parameter is the length of convolution, NOT the length
    % of the original signals! Super-important!
    % FFT of wavelet, and amplitude-normalize in the frequency domain
    cmwX = fft(cmw,nConv);
    cmwX = cmwX ./ max(cmwX);
    % FFT of data
    dataX = fft(data,nConv)';
    conv_res = dataX.*cmwX; %% --> we have just performed a convolution in the frequency domain (potentially the coolest thing ever) 
    % now back to the time domain
    % cut 1/2 of the length of the wavelet from the beginning and from the end
    half_wav = floor( length(cmw)/2 )+1;
    % take inverse Fourier transform
    conv_res_timedomain = ifft(conv_res);
    conv_res_timedomain = conv_res_timedomain(half_wav-1:end-half_wav);
    conv_res_timedomain = conv_res_timedomain(100:end-100); %remove edge artifacts from convolution by trimming 100 samples from both ends
    phase_timeSeries = angle(conv_res_timedomain);
    amplitude_timeSeries = abs(conv_res_timedomain);
    filtered_timeSeries = real(conv_res_timedomain);
    power_timeSeries = abs(conv_res_timedomain).^2;
end

