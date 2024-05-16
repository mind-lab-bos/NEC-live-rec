clear variables;
% study_path = '/Users/arun/Library/Mobile Documents/com~apple~CloudDocs/Documents/performance study MINDLab'; %% Only thing to update!
%addpath(strcat(study_path, '/Analysis/')); %??
conditions = {'LiveFast', 'LiveSlow', 'RecordedFast', 'RecordedSlow'};
participant_table = readtable("Participants.xlsx");
set(0,'defaultfigurecolor',[1 1 1])
 % set(0,'defaulttextinterpreter','latex')
% set(groot,'defaultAxesFontName','Verdana')
% set(groot,'defaultFontName','Verdana')
set(0,'defaultAxesFontSize',20)
set(0, 'DefaultLineLineWidth', 2); %for styling
%eeglab;
%%
numBins = 46; % number of frequencies (frequency resolution)
lowFreq = .2;
highFreq = 45.2;
cycles = 7; % must be odd number
fois = lowFreq:(highFreq-lowFreq+1)/numBins:highFreq; % frequency values at which we calculate PLV
%%
aggregatePLVs = zeros(size(participant_table,1),4,numBins);
aggregateCoherence = zeros(size(participant_table,1),4,numBins);

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

        % Resample audio
        EEG = pop_loadset('filename','finished(ica).set','filepath', preprocessed_data_path);
        eegdata = squeeze(mean(EEG.data,1))';
        audio_signal = load(strcat(analysis_file_path,'/cochlear_envelope.mat'),"finalSignal").finalSignal;
        audio_fs = 500;

        % Trim
        if (size(eegdata,1) > size(audio_signal,1))
            eegdata = eegdata(1:size(audio_signal));
        else
            audio_signal = audio_signal(1:size(eegdata));
        end

        % Calc PLV and Coherence

        for k = 1:numBins
            cmw = designMorlet(fois(k), cycles, EEG.srate); % build a wavelet
            [EEG_phase, EEG_amplitude, EEG_filtered, EEG_power] = waveletConvolution(eegdata,cmw); % convolve with EEG
            [Audio_phase, Audio_amplitude, Audio_filtered, Audio_power] = waveletConvolution(audio_signal,cmw); % convolve with Audio
            PLVs(k) = phaseCoherence(EEG_phase, Audio_phase);
            Coherences(k) = Coherence(Audio_phase, Audio_amplitude, EEG_phase, EEG_amplitude);
        end
        aggregatePLVs(i,conditionNum,:) = PLVs;
        aggregateCoherence(i,conditionNum,:) = Coherences;
    end
end
%  %% Plot Examle Phase vs Amplitude vs Filter vs Raw
% figure;
% set(gcf, "color", "white"); 
% hold on;
% plot(eegdata(2000:2700)*2, "DisplayName","Raw","LineWidth",2)
% 
%  plot(EEG_filtered(2000:2700)*10, "DisplayName","Filtered at 4.2 Hz","LineWidth",2)
%  plot(EEG_amplitude(2000:2700)*10, "DisplayName","Amplitude at 4.2 Hz","LineWidth",2)
%  plot(EEG_phase(2000:2700)*2, "DisplayName","Phase at 4.2 Hz","LineWidth",2)
% legend;
% title("Results of Wavelet Convolution from Raw EEG Signal");
% yticks([])
% xlabel("Time")
% hold off
%% Plot
PLV_Live = squeeze(mean(aggregatePLVs(:,1:2,:),[1 2]));
PLV_Rec = squeeze(mean(aggregatePLVs(:,3:4,:),[1 2]));

figure;
hold on;
plot(fois, PLV_Live, 'LineWidth',2, 'Color',[1 0 0]);
plot(fois, PLV_Rec, 'LineWidth',2, 'Color',[0 0 1]);
legend('live', 'rec')
hold off;

%% Shaded errors ALL
PLV_LiveFast = squeeze(mean(aggregatePLVs(:,1,:),1));
PLV_LiveSlow = squeeze(mean(aggregatePLVs(:,2,:),1));
PLV_RecFast = squeeze(mean(aggregatePLVs(:,3,:),1));
PLV_RecSlow = squeeze(mean(aggregatePLVs(:,4,:),1));
eLiveFast = calcErrors(squeeze(aggregatePLVs(:,1,:)));
eLiveSlow = calcErrors(squeeze(aggregatePLVs(:,2,:)));
eRecFast = calcErrors(squeeze(aggregatePLVs(:,3,:)));
eRecSlow = calcErrors(squeeze(aggregatePLVs(:,4,:)));


% figure;
% tiledlayout('flow')
nexttile;
hold on;

hp = patch([8 8 9.6 9.6],[0 .14 .14 0],'k',...
    'facecolor',[.5 .5 .5],'facealpha', 0.2,'edgecolor','none', 'HandleVisibility', 'off') ;
hp2 = patch([0.6 0.6 1 1],[0 .14 .14 0],'k',...
    'facecolor',[.5 .5 .5],'facealpha', 0.2,'edgecolor','none', 'HandleVisibility', 'off') ;

shadedErrorBar(fois,PLV_LiveFast,eLiveFast,'lineprops',{'Color',[.8 0 0],'LineWidth',4,'DisplayName',"Live Fast"})
shadedErrorBar(fois,PLV_LiveSlow,eLiveSlow,'lineprops',{'Color',[.7 .4 0],'LineWidth',4,'DisplayName',"Live Slow"})
shadedErrorBar(fois,PLV_RecFast,eRecFast,'lineprops',{'Color',[0 .6 0.4],'LineWidth',4,'DisplayName',"Recorded Fast"})
shadedErrorBar(fois,PLV_RecSlow,eRecSlow,'lineprops',{'Color',[0 0 .8],'LineWidth',4,'DisplayName',"Recorded Slow"})

xline(8,'-',{'8 Hz'}, 'LineWidth',1, 'HandleVisibility', 'off', 'color', [.8, 0, .2], 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center','LabelOrientation', 'align');
xline(9.6,'-',{'9.6 Hz'}, 'LineWidth',1,'HandleVisibility', 'off', 'color', [.8, 0, .2], 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center','LabelOrientation', 'align');
xline(0.6,'-',{'0.6 Hz'},'LineWidth',1, 'HandleVisibility', 'off', 'color', [.8, 0, .2], 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center','LabelOrientation', 'align');
xline(1,'-',{'1 Hz'},'LineWidth',1, 'HandleVisibility', 'off', 'color', [.8, 0, .2], 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center','LabelOrientation', 'align');

legend;
hold off;
 title('PLVs By Condition')
xlabel('Frequency [Hz]')
ylabel('Phase Coherence')
xlim([.2 12.2]);
ylim([0 0.14])
set(gcf,'color','w');
fontsize(25, "points")


% xticks(sort([xticks .8 8.4])); xline([.8 8.4],"HandleVisibility","off","LineWidth",2);
%% Anova All at Peak
    frequency_indeces = 3:5;
runAnova(squeeze(mean(aggregatePLVs(:,1,frequency_indeces),3)), squeeze(mean(aggregatePLVs(:,2,frequency_indeces),3)), squeeze(mean(aggregatePLVs(:,3,frequency_indeces),3)), squeeze(mean(aggregatePLVs(:,4,frequency_indeces),3)))

%% Q score PLV for fast condition
questions = {'Engagement','Enjoyment','Familiarity','Pleasure','Focus','Investment','Distraction (Reverse-Scored)','Spontaneity'};
conditions = 1:4;
figure;
t = tiledlayout('flow');
title(t, 'PLVs (Fast Conditions) by Questionnaire Responses');
% xlabel(t, 'Frequency[Hz]');
% ylabel(t, 'PLV');
for Q_num = 1:8
    nexttile;
    Q_data = squeeze(aggregate(:,conditions,Q_num));
    q_low = find(Q_data<mean(Q_data,'all'));
    q_high = find(Q_data>mean(Q_data,'all')); % mean split
    n_low = length(q_low);
    n_high = length(q_high);
    [rlow, clow] = ind2sub([21,4],q_low);
    [rhigh, chigh] = ind2sub([21,4],q_high);
    lowPLV = zeros(n_low,size(aggregatePLVs,3));
    highPLV = zeros(n_high,size(aggregatePLVs,3));
    numLow = 0;
    numHigh = 0;
    for i = 1:length(q_low)
        if (clow(i) == 1 || clow(i) == 3)
            numLow = numLow + 1;
            lowPLV(i,:) = squeeze(aggregatePLVs(rlow(i),clow(i),:));
        end
    end
    for i = 1:length(q_high)
        if (chigh(i) == 1 || chigh(i) == 3)
            numHigh = numHigh + 1;
            highPLV(i,:) = squeeze(aggregatePLVs(rhigh(i),chigh(i),:));
        end
    end
    lowPLV = lowPLV(all(lowPLV,2),:);
    highPLV = highPLV(all(highPLV,2),:);
    errlow = std(lowPLV,0,1)/sqrt(size(lowPLV,1));
    errhigh = std(highPLV,0,1)/sqrt(size(highPLV,1));
    
    switch Q_num
        case 1
            low_Engagement = lowPLV;
            high_Engagement = highPLV;
        case 3
            low_Familiarity = lowPLV;
            high_Familiarity = highPLV;
        case 4
            low_Pleasure = lowPLV;
            high_Pleasure = highPLV;
        case 8
            low_Spontaneity = lowPLV;
            high_Spontaneity= highPLV;
    end

    hold on
    shadedErrorBar(fois,mean(highPLV),errhigh,'lineprops',{'Color',[0.8, .2, 0],'LineWidth',2,'DisplayName',"High Rating n=" + num2str(numHigh)})
    shadedErrorBar(fois,mean(lowPLV),errlow,'lineprops',{'Color',[0, .2, 0.8],'LineWidth',2,'DisplayName',"Low Rating n=" + num2str(numLow)})
    legend;
    hold off;
    title(questions(Q_num))
    xlabel("Frequency [Hz]")
    ylabel("PLV")
    xlim([.2 20.2]);
    legend;
    % xticks(sort([xticks .8 3.2])); xline([.8 3.2],"HandleVisibility","off","LineWidth",2);
    yticks(0:.02:.1)
    % if Q_num == 1
    %     frequency_indeces = fois >= 8 & fois<=9.6;
    %     sample1 = mean(lowPLV(:, frequency_indeces), 2)
    %     sample2 = mean(highPLV(:, frequency_indeces), 2)
    %     [h, p] = ttest2(sample1, sample2);
    % end
end



Q_data = squeeze(mean(aggregate(:,conditions,:),3));
q_low = find(Q_data<mean(Q_data,'all'));
q_high = find(Q_data>mean(Q_data,'all'));
n_low = length(q_low);
n_high = length(q_high);
[rlow, clow] = ind2sub([21,4],q_low);
[rhigh, chigh] = ind2sub([21,4],q_high);
lowPLV = zeros(n_low,size(aggregatePLVs,3));
highPLV = zeros(n_high,size(aggregatePLVs,3));
numLow = 0;
numHigh = 0;
for i = 1:length(q_low)
    if (clow(i) == 1 || clow(i) == 3)
        numLow = numLow + 1;
        lowPLV(i,:) = squeeze(aggregatePLVs(rlow(i),clow(i),:));
    end
end
for i = 1:length(q_high)
    if (chigh(i) == 1 || chigh(i) == 3)
        numHigh = numHigh + 1;
        highPLV(i,:) = squeeze(aggregatePLVs(rhigh(i),chigh(i),:));
    end
end
lowPLV = lowPLV(all(lowPLV,2),:);
highPLV = highPLV(all(highPLV,2),:);

disp(num2str(size(lowPLV,1)))

errlow = std(lowPLV,0,1)/sqrt(size(lowPLV,1));
errhigh = std(highPLV,0,1)/sqrt(size(highPLV,1));

nexttile;
hold on;
% plot(fois, PLV_Live, 'LineWidth',2, 'Color',[.5 0 0], 'DisplayName','Live',"LineWidth",4);
% errorbar(fois,PLV_Live,eLive,'Color',[.5 0 0 .5],'HandleVisibility','off',"LineWidth",2);
hold on
shadedErrorBar(fois,mean(highPLV),errhigh,'lineprops',{'Color',[0.8, .2, 0],'LineWidth',2,'DisplayName',"High Rating n=" + num2str(numHigh)})
shadedErrorBar(fois,mean(lowPLV),errlow,'lineprops',{'Color',[0, .2, 0.8],'LineWidth',2,'DisplayName',"Low Rating n=" + num2str(numLow)})
legend;
hold off;

title('Aggregate')
xlabel('Frequency [Hz]')
ylabel('PLV')
xlim([.2 20.2]);
set(gcf,'color','w');
fontsize(20, "points")
% xticks(sort([xticks .8 3.2])); xline([.8 3.2],"HandleVisibility","off","LineWidth",2);
yticks(0:.02:.1)
%% T test (Welch's>) for aggregate high vs low 
lowPLV = low_Spontaneity; highPLV = high_Spontaneity;
frequency_indeces =  fois == 8.6; %fois >= 8 & fois<=9.6;
[h, p] = ttest2(mean(lowPLV(:, frequency_indeces), 2), mean(highPLV(:, frequency_indeces), 2), 'Vartype', 'unequal');
%% Scatter/BoxPlot high v low

% lowPLV = low_Familiarity; highPLV = high_Familiarity;

lows = mean(lowPLV(:, frequency_indeces), 2);
% xlows = zeros(length(lows),1);

highs = mean(highPLV(:, frequency_indeces), 2);
% xhighs = ones(length(highs),1);
% figure;
% hold on
% scatter(xlows, lows);
% xlim([-1, 2])
% 
% scatter(xhighs, highs);
% 
% plot([0, 1], [mean(lows), mean(highs)])
% 
% hold off

figure
hold on
xlim([.5, 2.5])
ylim([0, .32])
title("pleasure")
pl = [lows; highs];
pl_f = [repmat({'Low'},length(lows),1); repmat({'High'},length(highs),1)];
boxplot(pl,pl_f)
hold off
%% Shaded errors Live v Rec % calculate with 7 cycle 101 bin avg Chan
PLV_Live = squeeze(mean((aggregatePLVs(:,1,:)+aggregatePLVs(:,2,:))/2,1));
PLV_Rec = squeeze(mean((aggregatePLVs(:,3,:)+aggregatePLVs(:,4,:))/2,1));
eLive = calcErrors(squeeze((aggregatePLVs(:,1,:)+aggregatePLVs(:,2,:))/2));
eRec = calcErrors(squeeze((aggregatePLVs(:,3,:)+aggregatePLVs(:,4,:))/2));
% figure;
% tiledlayout('flow')
nexttile;
hold on;
% plot(fois, PLV_Live, 'LineWidth',2, 'Color',[.5 0 0], 'DisplayName','Live',"LineWidth",4);
% errorbar(fois,PLV_Live,eLive,'Color',[.5 0 0 .5],'HandleVisibility','off',"LineWidth",2);
shadedErrorBar(fois,PLV_Live,eLive,'lineprops',{'Color',[.5 0 0 .5],'LineWidth',4,'DisplayName',"Live"})

% plot(fois, PLV_Rec, 'LineWidth',2, 'Color',[0 .5 0], 'DisplayName','Rec',"LineWidth",4);
shadedErrorBar(fois,PLV_Rec,eRec,'lineprops',{'Color',[0 .5 0],'LineWidth',4,'DisplayName',"Recorded"})
legend;
hold off;
 title('Cerebro-Acoustic Phase Locking Between Live and Recorded Conditions')
xlabel('Frequency [Hz]')
ylabel('Phase Coherence')
xlim([.2 4.2]);
set(gcf,'color','w');
fontsize(25, "points")
xticks(sort([xticks .8 8.4])); xline([.8 8.4],"HandleVisibility","off","LineWidth",2);
yticks(0:.02:.1)
%% Shaded errors Fast v Slow % calculate with 7 cycle 101 bin avg Chan
PLV_Fast = squeeze(mean((aggregatePLVs(:,1,:)+aggregatePLVs(:,3,:))/2,1));
PLV_Slow = squeeze(mean((aggregatePLVs(:,2,:)+aggregatePLVs(:,4,:))/2,1));
eFast = calcErrors(squeeze((aggregatePLVs(:,1,:)+aggregatePLVs(:,3,:))/2));
eSlow = calcErrors(squeeze((aggregatePLVs(:,2,:)+aggregatePLVs(:,4,:))/2));
nexttile;
hold on;
% plot(fois, PLV_Live, 'LineWidth',2, 'Color',[.5 0 0], 'DisplayName','Live',"LineWidth",4);
% errorbar(fois,PLV_Live,eLive,'Color',[.5 0 0 .5],'HandleVisibility','off',"LineWidth",2);
shadedErrorBar(fois,PLV_Fast,eFast,'lineprops',{'Color',[.5 0 .5],'LineWidth',4,'DisplayName',"Fast"})

% plot(fois, PLV_Rec, 'LineWidth',2, 'Color',[0 .5 0], 'DisplayName','Rec',"LineWidth",4);
shadedErrorBar(fois,PLV_Slow,eSlow,'lineprops',{'Color',[.5 .5 .1],'LineWidth',4,'DisplayName',"Slow"})
legend;
hold off;
 title(' Fast vs. Slow Stimuli')
xlabel('Frequency [Hz]')
ylabel('Phase Coherence')
xlim([.2 4.2]);
set(gcf,'color','w');
fontsize(25, "points")
xticks(sort([xticks, 0.8, 8.4])); xline([0.8, 8.4],"HandleVisibility","off","LineWidth",2);

%% Multiple Paired T-Test Live/Rec
% 
% pvals = zeros(length(fois),1);
% for i= 1: length(pvals)
%     [~, pvals(i)] = ttest(mean(aggregatePLVs(:,[1, 2],i),2),mean(aggregatePLVs(:,[3 4],i),2));
% end








%% Live v Rec
PLV_Live = squeeze(mean((aggregatePLVs(:,1,:)+aggregatePLVs(:,2,:))/2,1));
PLV_Rec = squeeze(mean((aggregatePLVs(:,3,:)+aggregatePLVs(:,4,:))/2,1));
eLive = calcErrors(squeeze((aggregatePLVs(:,1,:)+aggregatePLVs(:,2,:))/2));
eRec = calcErrors(squeeze((aggregatePLVs(:,3,:)+aggregatePLVs(:,4,:))/2));
f = figure;
hold on;
plot(fois, PLV_Live, 'LineWidth',2, 'Color',[.5 0 0], 'DisplayName','Live',"LineWidth",4);
errorbar(fois,PLV_Live,eLive,'Color',[.5 0 0 .5],'HandleVisibility','off',"LineWidth",2);

plot(fois, PLV_Rec, 'LineWidth',2, 'Color',[0 .5 0], 'DisplayName','Rec',"LineWidth",4);
errorbar(fois,PLV_Rec,eRec,'Color',[0 .5 0 .5],'HandleVisibility','off',"LineWidth",2);
legend;
hold off;
 title('PLVs PLVs from .2-20.2 Hz Live vs. Rec')
xlabel('Frequency [Hz]')
ylabel('Phase Coherence')
xlim([.2 20.2]);
set(gcf,'color','w');
% fontsize(f, 40, "points")
xticks(sort([xticks .8 3.2])); xline([.8 3.2],"HandleVisibility","off","LineWidth",3);
%% Fast vs. Slow PLV
PLV_Fast = squeeze(mean((aggregatePLVs(:,1,:)+aggregatePLVs(:,3,:))/2,1));
PLV_Slow = squeeze(mean((aggregatePLVs(:,2,:)+aggregatePLVs(:,4,:))/2,1));
eFast = calcErrors(squeeze((aggregatePLVs(:,1,:)+aggregatePLVs(:,3,:))/2));
eSlow = calcErrors(squeeze((aggregatePLVs(:,2,:)+aggregatePLVs(:,4,:))/2));
figure; 
hold on;
plot(fois, PLV_Fast, 'LineWidth',2, 'Color',[.5 0 .5], 'DisplayName','Fast');
errorbar(fois,PLV_Fast,eFast,'Color',[.5 0 .5 .5],'HandleVisibility','off');

plot(fois, PLV_Slow, 'LineWidth',2, 'Color',[.5 .5 1], 'DisplayName','Slow');
errorbar(fois,PLV_Slow,eSlow,'Color',[.5 .5 1 .5],'HandleVisibility','off');
legend;
hold off;
title('PLVs from .2-20.2 Hz Fast vs. Slow')
xlabel('Frequency [Hz]')
ylabel('Phase Coherence')
xlim([.2 20.2]);
xticks(sort([xticks .8 3.2])); xline(3.2,"HandleVisibility","off"); xline(.8,"HandleVisibility","off")
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


%% Errors across frequencies
LiveFast = squeeze(aggregatePLVs(:,1,:)); 
LiveSlow = squeeze(aggregatePLVs(:,2,:));
RecFast = squeeze(aggregatePLVs(:,3,:));
RecSlow = squeeze(aggregatePLVs(:,4,:));

eLiveFast = calcErrors(LiveFast); eLiveSlow = calcErrors(LiveSlow);
eRecFast = calcErrors(RecFast); eRecSlow = calcErrors(RecSlow);


%%
figure;
hold on;
plot(fois, squeeze(mean(aggregatePLVs(:,1,:),1)), 'LineWidth',1, 'Color',[1 0 0 ],'DisplayName','Live Fast','LineWidth',2);
errorbar(fois,squeeze(mean(aggregatePLVs(:,1,:),1)),eLiveFast,'Color',[1 0 0 .5],'HandleVisibility','off');

plot(fois, squeeze(mean(aggregatePLVs(:,2,:),1)), 'LineWidth',1, 'Color',[1 .5 0], 'DisplayName','Live Slow','LineWidth',2);
errorbar(fois,squeeze(mean(aggregatePLVs(:,2,:),1)),eLiveSlow,'Color',[1 .5 0 .5],'HandleVisibility','off');

plot(fois, squeeze(mean(aggregatePLVs(:,3,:),1)), 'LineWidth',1, 'Color',[0 0 1 ], 'DisplayName','Recorded Fast','LineWidth',2);
errorbar(fois,squeeze(mean(aggregatePLVs(:,3,:),1)),eRecFast,'Color',[0 0 1 .5],'HandleVisibility','off');

plot(fois, squeeze(mean(aggregatePLVs(:,4,:),1)), 'LineWidth',1, 'Color',[0 1 0],'DisplayName','Recorded Slow','LineWidth',2);
errorbar(fois,squeeze(mean(aggregatePLVs(:,4,:),1)),eRecSlow,'Color',[0 1 0 .5],'HandleVisibility','off');
title('PLVs from .2-20.2 Hz')
xlabel('Frequency [Hz]')
ylabel('Phase Coherence')
xlim([0.2 20.2]);
xticks(sort([xticks .8 3.2])); xline(3.2,"HandleVisibility","off"); xline(.8,"HandleVisibility","off")
fontsize('increase');
legend;
hold off;

%%
figure;
hold on;
plot(fois, squeeze(mean(aggregatePLVs(:,1,:),2)), 'LineWidth',1, 'Color',[1 0 0 .5]);
plot(fois, squeeze(mean(aggregatePLVs(:,2,:),2)), 'LineWidth',1, 'Color',[1 .5 0 .5]);
plot(fois, squeeze(mean(aggregatePLVs(:,3,:),2)), 'LineWidth',1, 'Color',[0 0 1 .5]);
% legend('liveFast', 'recFast')

% [up, lo] = envelope(audioRS,600);
% env = (up+lo)/2;


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
    title('PLVs at 1.2 Hz by Condition');
    hold off
%% PLOT COHERENCES
figure;
hold on;


plot(fois, squeeze(mean(aggregateCoherence(:,1,:),1)), 'LineWidth',1, 'Color',[1 0 0 ],'DisplayName','Live Fast');
errorbar(fois,squeeze(mean(aggregateCoherence(:,1,:),1)),eLiveFast,'Color',[1 0 0 .5],'HandleVisibility','off');


plot(fois, squeeze(mean(aggregateCoherence(:,2,:),1)), 'LineWidth',1, 'Color',[1 .5 0], 'DisplayName','Live Slow');
errorbar(fois,squeeze(mean(aggregateCoherence(:,2,:),1)),eLiveSlow,'Color',[1 .5 0 .5],'HandleVisibility','off');

plot(fois, squeeze(mean(aggregateCoherence(:,3,:),1)), 'LineWidth',1, 'Color',[0 0 1 ], 'DisplayName','Recorded Fast');
errorbar(fois,squeeze(mean(aggregateCoherence(:,3,:),1)),eRecFast,'Color',[0 0 1 .5],'HandleVisibility','off');

plot(fois, squeeze(mean(aggregateCoherence(:,4,:),1)), 'LineWidth',1, 'Color',[0 1 0],'DisplayName','Recorded Slow');
 errorbar(fois,squeeze(mean(aggregateCoherence(:,4,:),1)),eRecSlow,'Color',[0 1 0 .5],'HandleVisibility','off');
title('Coherence at FOIs from .2-45.2 Hz')
xlabel('Frequency')
ylabel('Coherence')
xlim([0 45.2]);
fontsize('increase');
legend;
hold off;

%% Fast vs. Slow
C_Fast = squeeze(mean((aggregateCoherence(:,1,:)+aggregateCoherence(:,3,:))/2,1));
C_Slow = squeeze(mean((aggregateCoherence(:,2,:)+aggregateCoherence(:,4,:))/2,1));
eCFast = calcErrors(squeeze((aggregateCoherence(:,1,:)+aggregateCoherence(:,3,:))/2));
eCSlow = calcErrors(squeeze((aggregateCoherence(:,2,:)+aggregateCoherence(:,4,:))/2));
figure;
hold on;
plot(fois, C_Fast, 'LineWidth',2, 'Color',[.5 0 .5], 'DisplayName','Fast');
errorbar(fois,C_Fast,eCFast,'Color',[.5 0 .5 .5],'HandleVisibility','off');

plot(fois, C_Slow, 'LineWidth',2, 'Color',[.5 .5 1], 'DisplayName','Slow');
errorbar(fois,C_Slow,eCSlow,'Color',[.5 .5 1 .5],'HandleVisibility','off');
legend;
hold off;
title('Coherence at FOIs from .2-45.2 Hz between Tempo Conditions')
xlabel('Frequency')
ylabel('Phase Coherence')
xlim([0 45.2]);
%% Live vs Rec
C_Live = squeeze(mean(aggregateCoherence(:,1:2,:),[1 2]));
C_Rec= squeeze(mean(aggregateCoherence(:,3:4,:),[1 2]));
eCLive = calcErrors(squeeze((aggregateCoherence(:,1,:)+aggregateCoherence(:,2,:))/2));
eCRec = calcErrors(squeeze((aggregateCoherence(:,3,:)+aggregateCoherence(:,4,:))/2));
figure;
hold on;
plot(fois, C_Live, 'LineWidth',2, 'Color',[.5 0 0 1], 'DisplayName','Live');
errorbar(fois,C_Live,eCLive,'Color',[.5 0 0 .5],'HandleVisibility','off');

plot(fois, C_Rec, 'LineWidth',2, 'Color',[0 .5 0 1], 'DisplayName','Rec');
errorbar(fois,C_Rec,eCRec,'Color',[0 .5 0 .5],'HandleVisibility','off');
legend;
hold off;
title('Coherence at FOIs from .2-45.2 Hz between Performance Mediums')
xlabel('Frequency')
ylabel('Phase Coherence')
xlim([0 45.2]);


%%
function PLV = phaseCoherence(angles1, angles2) %%phase coherence for two angle time series %% Harding Cortical Tracking music and speech
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
    power_timeSeries = amplitude_timeSeries.^2;
end


function errors = calcErrors(dataMatrix) %column is participats and row is data points
    errors = zeros(1, size(dataMatrix,2));
    for i = 1: size(dataMatrix,2)
        errors(i) = std(dataMatrix(:,i))/sqrt(size(dataMatrix,1));
    end
end

function runAnova(LiveFast, LiveSlow, RecFast, RecSlow)
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
    
    %8-9.6 Hz
    A = zeros(42,2); %% row factor is tempo f/s, column factor is live/rec
    A(1:21,1) = LiveFast;
    A(1:21,2) = RecFast;
    A(22:42,1) = LiveSlow;
    A(22:42,2) = RecSlow;
    
    [~,~,stats] = anova2(A,21); figure;
    c1 = multcompare(stats); figure;
    c2 = multcompare(stats,"Estimate","row");
    % tbl1 = array2table(c1,"VariableNames", ...
    %     ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

end

