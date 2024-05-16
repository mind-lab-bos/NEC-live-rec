% Toolboxes - statistics and machine learning, signal processing
 clear; %clc;
study_path = '/Users/arun/Library/Mobile Documents/com~apple~CloudDocs/Documents/performance study MINDLab'; %% Only thing to update!
addpath(strcat(study_path, '/Analysis/EEGlab')); eeglab;
%% 
participant_table = readtable("Participants.xlsx");
%%
% savePSD(EEG_ZacLive, save_file_path, 'Wavelet');
%%
% iterateParticipants(participant_table, study_path, 'Pwelch');
%%
[LiveFast, LF_all] = avgParticipants(participant_table, study_path, 'Pwelch', 'LiveFast');
[LiveSlow, LS_all] = avgParticipants(participant_table, study_path, 'Pwelch', 'LiveSlow');
[RecordedFast, RF_all] = avgParticipants(participant_table, study_path, 'Pwelch', 'RecordedFast');
[RecordedSlow, RS_all] = avgParticipants(participant_table, study_path, 'Pwelch', 'RecordedSlow');
%% PSD All Participants
figure;
tiledlayout('flow');
for i = 1:21
    nexttile;
    hold on;
    plot(LiveFast(:,1),log10(squeeze(LF_all(i,:,2))),'r','LineWidth',2)
    plot(LiveSlow(:,1),log10(squeeze(LS_all(i,:,2))),'g', 'LineWidth',2)
    plot(RecordedFast(:,1),log10(squeeze(RF_all(i,:,2))),'b', 'LineWidth',2)
    plot(RecordedSlow(:,1),log10(squeeze(RS_all(i,:,2))),'Color', [1 .4 0 1], 'LineWidth',2)
    xlim([0 58])
    title(strcat('Participant-',num2str(i)));
    hold off;
end
leg = legend('LiveFast','LiveSlow','RecordedFast','RecordedSlow');
leg.Layout.Tile = 'north';
set(gcf,"color","white")
%% Plot PSD PWELCH NO FOOOF
figure;
hold on;
title("Pwelch PSD for All Conditions and All Participants")
plot(LiveFast(:,1),log10(LiveFast(:,2)),'r','LineWidth',2)
plot(LiveSlow(:,1),log10(LiveSlow(:,2)),'g', 'LineWidth',2)
plot(RecordedFast(:,1),log10(RecordedFast(:,2)),'b', 'LineWidth',2)
plot(RecordedSlow(:,1),log10(RecordedSlow(:,2)),'Color', [1 .4 0 1], 'LineWidth',2)

plot(LiveFast(:,1),log10(squeeze(LF_all(:,:,2))),'Color',[1 0 0 .2],'LineWidth',1)
plot(LiveSlow(:,1),log10(squeeze(LS_all(:,:,2))),'Color',[0 1 0 .2], 'LineWidth',1)
plot(RecordedFast(:,1),log10(squeeze(RF_all(:,:,2))),'Color',[0 0 1 .2], 'LineWidth',1)
plot(RecordedSlow(:,1),log10(squeeze(RS_all(:,:,2))),'Color',[1 .4 0 .2],'LineWidth',1)

legend('LiveFast','LiveSlow','RecordedFast','RecordedSlow')
xlim([.2 58])
xlabel("Frequency [Hz]")
ylabel("Log(V^2/Hz)")
set(gcf, "color", "white"); fontsize("increase")


axes('Position',[0.293567251461988,0.684210526315789,0.278362573099415,0.195488721804513])
box on;
hold on;
xrange = find(LiveFast(:,1)<12 & LiveFast(:,1)>8);
plot(LiveFast(xrange,1), log10(LiveFast(xrange,2)),'r','LineWidth',1);
plot(LiveSlow(xrange,1),log10(LiveSlow(xrange,2)),'g', 'LineWidth',1)
plot(RecordedFast(xrange,1),log10(RecordedFast(xrange,2)),'b', 'LineWidth',1)
plot(RecordedSlow(xrange,1),log10(RecordedSlow(xrange,2)),'Color', [1 .4 0 1], 'LineWidth',1)
hold off;
axis tight;
box off;

axes('Position',[0.166081871345029,0.312781954887218,0.112280701754386,0.174436090225564])
box on;
hold on;
xrange = find(LiveFast(:,1)<8 & LiveFast(:,1)>4);
plot(LiveFast(xrange,1), log10(LiveFast(xrange,2)),'r','LineWidth',1);
plot(LiveSlow(xrange,1),log10(LiveSlow(xrange,2)),'g', 'LineWidth',1)
plot(RecordedFast(xrange,1),log10(RecordedFast(xrange,2)),'b', 'LineWidth',1)
plot(RecordedSlow(xrange,1),log10(RecordedSlow(xrange,2)),'Color', [1 .4 0 1], 'LineWidth',1)
hold off;
axis tight;
box off;

axes('Position',[0.421052631578947,0.485714285714286,0.109941520467836,0.11578947368421])
box on;
hold on;
xrange = find(LiveFast(:,1)<20 & LiveFast(:,1)>16);
plot(LiveFast(xrange,1), log10(LiveFast(xrange,2)),'r','LineWidth',1);
plot(LiveSlow(xrange,1),log10(LiveSlow(xrange,2)),'g', 'LineWidth',1)
plot(RecordedFast(xrange,1),log10(RecordedFast(xrange,2)),'b', 'LineWidth',1)
plot(RecordedSlow(xrange,1),log10(RecordedSlow(xrange,2)),'Color', [1 .4 0 1], 'LineWidth',1)
hold off;
axis tight;
box off;

hold off

figure;
title("Pwelch PSD for All Conditions")

hold on;
plot(LiveFast(:,1),log10(LiveFast(:,2)),'r','LineWidth',2)
plot(LiveSlow(:,1),log10(LiveSlow(:,2)),'g', 'LineWidth',2)
plot(RecordedFast(:,1),log10(RecordedFast(:,2)),'b', 'LineWidth',2)
plot(RecordedSlow(:,1),log10(RecordedSlow(:,2)),'Color', [1 .4 0 1], 'LineWidth',2)

legend('LiveFast','LiveSlow','RecordedFast','RecordedSlow')
xlim([0 58])
hold off;
%% Plot Normalized PSDs

figure;
hold on;
title("Normalized Pwelch PSD for All Conditions and All Participants")
plot(LiveFast(:,1),log10(normalize(LiveFast(:,2)','norm',1)'),'r','LineWidth',2)
plot(LiveSlow(:,1),log10(normalize(LiveSlow(:,2)','norm',1)'),'g', 'LineWidth',2)
plot(RecordedFast(:,1),log10(normalize(RecordedFast(:,2)','norm',1)'),'b', 'LineWidth',2)
plot(RecordedSlow(:,1),log10(normalize(RecordedSlow(:,2)','norm',1)'),'Color', [1 .4 0 1], 'LineWidth',2)

% plot(LiveFast(:,1),log10(normalize(squeeze(LF_all(:,:,2))','norm',1)'),'Color',[1 0 0 .2],'LineWidth',1)
% plot(LiveSlow(:,1),log10(normalize(squeeze(LS_all(:,:,2))','norm',1)'),'Color',[0 1 0 .2],'LineWidth',1)
% plot(RecordedFast(:,1),log10(normalize(squeeze(RF_all(:,:,2))','norm',1)'),'Color',[0 0 1 .2],'LineWidth',1)
% plot(RecordedSlow(:,1),log10(normalize(squeeze(RS_all(:,:,2))','norm',1)'),'Color',[1 .4 0 .2],'LineWidth',1)

legend('LiveFast','LiveSlow','RecordedFast','RecordedSlow')
xlim([0 58])
%% plot Topo
participantID = '230814_SSZO';
excerptNum = 1;
preprocessed_data_path = strcat(study_path,'/Data/Preprocessed_Data/',participantID, '/Excerpt_',num2str(excerptNum), '/');
EEG = pop_loadset('filename','finished(ica).set','filepath', preprocessed_data_path);

% Live - Rec
Live = (LiveFast+LiveSlow)/2;
Recorded = (RecordedFast+RecordedSlow)/2;
frequencies = Live(:,1);
avgPSD = squeeze(Live(:,2) - Recorded(:,2))';
pxx = (Live(:,3:end) - Recorded(:,3:end))';
freqRange = [8 12];

plotFreqTopo(freqRange, avgPSD, pxx, frequencies, EEG, 'avg')
title('Live - Rec at 8-12 Hz');

freqRange = [4 8];
plotFreqTopo(freqRange, avgPSD, pxx, frequencies, EEG, 'avg')
title('Live - Rec at 4-8 Hz');
%% Plot each condition by participant score
a = figure(); plotByScore(LF_all,sortPart_LF,[30 45], "LiveFast 30-45 Hz", EEG, [-2 2]);
cb = colorbar;
cb.Layout.Tile = "east";
savefig(a, 'LiveFast Gamma by survey');

b = figure(); plotByScore(LS_all,sortPart_LS,[30 45], "LiveSlow 30-45 Hz", EEG, [-2 2]);
cb = colorbar;
cb.Layout.Tile = "east";
savefig(b, 'LiveSlow Gamma by survey');

c = figure(); plotByScore(RF_all,sortPart_RF,[30 45], "RecordedFast 30-45 Hz", EEG,  [-2 2]);
cb = colorbar;
cb.Layout.Tile = "east";
savefig(c,'RecordedFast Gamma by survey');

d = figure(); plotByScore(RS_all,sortPart_RS,[30 45], "RecordedSlow 30-45 Hz", EEG,  [-2 2]);
cb = colorbar;
cb.Layout.Tile = "east";
savefig(d, 'RecordedSlow Gamma by survey');

LiveAvg = (LS_all + LF_all)/2;
RecAvg = (RS_all + RF_all)/2;

e = figure(); plotByScore(LiveAvg,sortPart_LS,[30 45], "Live 30-45 Hz", EEG, [-2 2]);
cb = colorbar;
cb.Layout.Tile = "east";
savefig(e,'Live Gamma by survey');

f = figure(); plotByScore(RecAvg,sortPart_RF,[30 45], "Recorded 30-45 Hz", EEG,  [-2 2]);
cb = colorbar;
cb.Layout.Tile = "east";
savefig(f,'Recorded Gamma by survey');

%% Plot by survey diff live/rec
freqRange = [8 12];
frequencies = squeeze(LiveAvg(1,:,1))';
g = figure();
t = tiledlayout(7,6);
% t.TileIndexing = "rowmajor";
for i = 1:21
    nexttile;
     avgPSDLive = squeeze(LiveAvg(sortPart_diff(i),:,2))';
     pxxLive = squeeze((LiveAvg(sortPart_diff(i),:,3:end)))';
    plotFreqTopo(freqRange, avgPSD, pxxLive, frequencies, EEG, 'avg')
    clim([-2 2]);
    title(strcat('Live, Participant - ', num2str(sortPart_diff(i))));

    nexttile;
     avgPSDRec = squeeze(RecAvg(sortPart_diff(i),:,2))';
     pxxRec = squeeze((RecAvg(sortPart_diff(i),:,3:end)))';
    plotFreqTopo(freqRange, avgPSD, pxxRec, frequencies, EEG, 'avg')
    clim([-2 2]);
    title(strcat('Rec, Participant - ', num2str(sortPart_diff(i))));

    % nexttile;
    % avgPSD = avgPSDLive - avgPSDRec;
    % pxx = pxxLive - pxxRec;
    % fRange = frequencies > freqRange(1) & frequencies < freqRange(2);
    % topoSlice = (mean(pxx(:,fRange),2));
    % topoplot(topoSlice, EEG.chanlocs,  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo); 
    % title(strcat('Diff, Participant - ', num2str(sortPart_diff(i))));
    % 
    % clim([-2 2]);
end 
cb = colorbar;
cb.Layout.Tile = "south";
t.Title.String = strcat('Live/Rec By Survey Diff Score (', num2str(freqRange(1)),'-',num2str(freqRange(2)), ' Hz)');
t.Title.FontWeight = 'bold';

savefig(g,'Live and Rec By Survey Diff Score 8-12 Hz');


%% Plot DIFF by survey diff live/rec
freqRange = [30 45];
frequencies = squeeze(LiveAvg(1,:,1))';
g = figure();
t = tiledlayout(7,3);
% t.TileIndexing = "rowmajor";
for i = 1:21
     avgPSDLive = squeeze(LiveAvg(sortPart_diff(i),:,2))';
     pxxLive = squeeze((LiveAvg(sortPart_diff(i),:,3:end)))';

     avgPSDRec = squeeze(RecAvg(sortPart_diff(i),:,2))';
     pxxRec = squeeze((RecAvg(sortPart_diff(i),:,3:end)))';

    nexttile;
    avgPSD = avgPSDLive - avgPSDRec;
    pxx = pxxLive - pxxRec;
    fRange = frequencies > freqRange(1) & frequencies < freqRange(2);
    topoSlice = (mean(pxx(:,fRange),2));
    topoplot(topoSlice, EEG.chanlocs,  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo); 
    title(strcat('Diff, Participant - ', num2str(sortPart_diff(i))));

    clim([-.5 .5]);
end 
cb = colorbar;
cb.Layout.Tile = "east";
t.Title.String = strcat('Live/Rec By Survey Diff Score (', num2str(freqRange(1)),'-',num2str(freqRange(2)), ' Hz)');
t.Title.FontWeight = 'bold';

savefig(g,'Live_Rec Diff By Survey Diff Score 30-45 Hz');



%  figure;
% tiledlayout('flow')
% hold on;
% for i = 1:21
%     nexttile;
%     avgPSD = squeeze(LF_all(sortPart_LF(i),:,2))';
%     pxx = squeeze((LF_all(sortPart_LF(i),:,3:end)))';
%     freqRange = [8 12];
%     plotFreqTopo(freqRange, avgPSD, pxx, frequencies, EEG, 'avg')
%     clim([-30 30]);
%     title(strcat('LiveFast - RecFast at 8-12 Hz Participant: ', num2str(i)));
% end
% 
%  figure;
% tiledlayout('flow')
% hold on;
% for i = 1:21
%     nexttile;
%     avgPSD = squeeze(LF_all(sortPart_LF(i),:,2))';
%     pxx = squeeze((LF_all(sortPart_LF(i),:,3:end)))';
%     freqRange = [8 12];
%     plotFreqTopo(freqRange, avgPSD, pxx, frequencies, EEG, 'avg')
%     clim([-30 30]);
%     title(strcat('LiveFast - RecFast at 8-12 Hz Participant: ', num2str(i)));
% end


%% 
% LiveFast - RecFast
    % % % frequencies = LiveFast(:,1);
    % % % avgPSD = squeeze(LiveFast(:,2) - RecordedFast(:,2))';
    % % % pxx = (LiveFast(:,3:end) - RecordedFast(:,3:end))';
    % % % freqRange = [8 12];
    % % % 
    % % % plotFreqTopo(freqRange, avgPSD, pxx, frequencies, EEG, 'avg')
    % % % title('LiveFast - RecFast at 8-12 Hz');
    % % % 
    % % % freqRange = [4 8];
    % % % plotFreqTopo(freqRange, avgPSD, pxx, frequencies, EEG, 'avg')
    % % % title('LiveFast - RecFast at 4-8 Hz');

%  figure;
% tiledlayout('flow')
% hold on;
% for i = 1:21
%     nexttile;
%     avgPSD = squeeze(LF_all(i,:,2) - RF_all(i,:,2))';
%     pxx = squeeze((LF_all(i,:,3:end) - RF_all(i,:,3:end)))';
%     freqRange = [8 12];
%     plotFreqTopo(freqRange, avgPSD, pxx, frequencies, EEG, 'avg')
%     title(strcat('LiveFast - RecFast at 8-12 Hz Participant: ', num2str(i)));
% 
% end

%%
% LiveSlow - RecSlow
frequencies = LiveSlow(:,1);
avgPSD = squeeze(LiveSlow(:,2) - RecordedSlow(:,2))';
pxx = (LiveSlow(:,3:end) - RecordedSlow(:,3:end))';
freqRange = [8 12];

plotFreqTopo(freqRange, avgPSD, pxx, frequencies, EEG, 'avg')
title('LiveSlow - RecSlow at 8-12 Hz');

freqRange = [4 8];
plotFreqTopo(freqRange, avgPSD, pxx, frequencies, EEG, 'avg')
title('LiveSlow - RecSlow at 4-8 Hz');

figure;
tiledlayout('flow')
hold on;
for i = 1:21
    nexttile;
    avgPSD = squeeze(LS_all(i,:,2) - RS_all(i,:,2))';
    pxx = squeeze((LS_all(i,:,3:end) - RS_all(i,:,3:end)))';
    freqRange = [8 12];
    plotFreqTopo(freqRange, avgPSD, pxx, frequencies, EEG, 'avg')
    title(strcat('LiveFast - RecFast at 8-12 Hz Participant: ', num2str(i)));

end

%%
frequencies = Live(:,1);
avgPSD = squeeze(Live(:,2))';
pxx = (Live(:,3:end))';
freqRange = [16 20];
plotFreqTopo(freqRange, avgPSD, pxx, frequencies, EEG, 'avg')
title('Live 16-20 Hz');
clim([-2.5 2.5])
frequencies = Recorded(:,1);
avgPSD = squeeze(Recorded(:,2))';
pxx = (Recorded(:,3:end))';
freqRange = [16 20];
plotFreqTopo(freqRange, avgPSD, pxx, frequencies, EEG, 'avg')
title('Recorded 16-20 Hz');
clim([-2.5 2.5])

%%
frequencies = Live(:,1);
avgPSD = squeeze(LiveSlow(:,2) - RecordedSlow(:,2))';
pxx = (LiveSlow(:,3:end) - RecordedSlow(:,3:end))';
freqRange = [16 20];
plotFreqTopo(freqRange, avgPSD, pxx, frequencies, EEG, 'avg')
title('Live-Rec 16-20 Hz');
clim([-.5 .5])
%%

participants = readlines('Participants.txt');
% writematrix(xmatrix);
eeglab;
EEG_ZacLive = pop_loadset('filename','finished(ica).set','filepath','/Users/arun/Documents/eeglab2023.0/preprocessed/Zac/live');
EEG_ZacRec = pop_loadset('filename','finished(ica).set','filepath','/Users/arun/Documents/eeglab2023.0/preprocessed/Zac/rec');
EEG_ZacRec2 = pop_loadset('filename','finished(ica).set','filepath','/Users/arun/Documents/eeglab2023.0/preprocessed/Zac/recNew');
EEG_AnjLive = pop_loadset('filename','live-ica.set','filepath','/Users/arun/Documents/eeglab2023.0/preprocessed/Anja/Live');
EEG_AnjRec = pop_loadset('filename','rec-ica.set','filepath','/Users/arun/Documents/eeglab2023.0/preprocessed/Anja/Recorded');
EEG_JakobLive = pop_loadset('filename','live-postICA.set','filepath','/Users/arun/Documents/eeglab2023.0/preprocessed/Jakob/Live');
EEG_JakobRec = pop_loadset('filename','rec_finalPreprocess.set','filepath','/Users/arun/Documents/eeglab2023.0/preprocessed/Jakob/Recorded');
EEG_EvaLive = pop_loadset('filename','live-ica.set','filepath','/Users/arun/Documents/eeglab2023.0/preprocessed/Eva/Live');
EEG_EvaRec = pop_loadset('filename','rec-ICA.set','filepath','/Users/arun/Documents/eeglab2023.0/preprocessed/Eva/Recorded');
close;


%%
clear; clc; close all;
eeglab;
EEG_LaineLiveFast = pop_loadset('filename','finished(ica).set','filepath','/Users/arun/Documents/eeglab2023.0/preprocess/Laine/Live_BachC');
EEG_LaineRecFast = pop_loadset('filename','finished(ica).set','filepath','/Users/arun/Documents/eeglab2023.0/preprocess/Laine/Rec_BachE');
EEG_LaineLiveSlow = pop_loadset('filename','finished(ica).set','filepath','/Users/arun/Documents/eeglab2023.0/preprocess/Laine/Live_BachA');
EEG_LaineRecSlow = pop_loadset('filename','finished(ica).set','filepath','/Users/arun/Documents/eeglab2023.0/preprocess/Laine/Rec_BachG');
%%

plotLiveRecorded(EEG_LaineLiveFast,EEG_LaineRecFast, 'Laine Fast');
plotLiveRecorded(EEG_LaineLiveSlow,EEG_LaineRecSlow, 'Laine Slow');
%%
[LaineavgLiveFast,LaineavgRecFast,LainepxxLiveFast,LainepxxRecFast,LainefLiveFast, LainefRecFast] = plotWaveletPSD(EEG_LaineLiveFast,EEG_LaineRecFast, 'Laine Fast');
%%
[LaineavgLiveSlow,LaineavgRecSlow,LainepxxLiveSlow,LainepxxRecSlow,LainefLiveSlow, LainefRecSlow] = plotWaveletPSD(EEG_LaineLiveSlow,EEG_LaineRecSlow, 'Laine Slow');
%%
plotLiveRecorded(EEG_ZacLive,EEG_ZacRec, 'Zac');
%%
hold on;
[rec2avg, rec2pxx, rec2f] = PSD(EEG_ZacRec2);
plotPSD(rec2f, rec2avg, rec2pxx, 'blue', [1 1 1 0]);
hold off;
%%
plotLiveRecorded(EEG_AnjLive,EEG_AnjRec, 'Anj');
%%
plotLiveRecorded(EEG_JakobLive,EEG_JakobRec, 'Jakob');
%%
plotLiveRecorded(EEG_EvaLive,EEG_EvaRec, 'Eva');
%%
[ZacavgLive,ZacavgRec,ZacpxxLive,ZacpxxRec,ZacfLive, ZacfRec] = plotWaveletPSD(EEG_ZacLive,EEG_ZacRec, 'Zac');
%%    
cwt(mean(EEG_ZacLive.data,1), "amor", EEG_ZacLive.srate); 
%%
cwt(mean(EEG_EvaLive.data,1), "amor", EEG_EvaLive.srate); 
%%
[AnjavLive,AnjavgRec,AnjpxxLive,AnjpxxRec,AnjfLive, AnjfRec] = plotWaveletPSD(EEG_AnjLive,EEG_AnjRec, 'Anj');
%%
[EvaavgLive,EvaavgRec,EvapxxLive,EvapxxRec,EvafLive, EvafRec] = plotWaveletPSD(EEG_EvaLive,EEG_EvaRec, 'Eva');
%%
[JakavgLive,JakavgRec,JakpxxLive,JakpxxRec,JakfLive, JakfRec] = plotWaveletPSD(EEG_JakobLive,EEG_JakobRec, 'Jakob');

%%
livemaxes = islocalmax(ZacavgLive);
recmaxes = islocalmax(ZacavgRec);
fRange = f > .5 & f < 3.5;
livefreqPeaks = fRange & livemaxes';
recfreqPeaks = fRange & recmaxes';

hold on;
plot(f(livemaxes)', log10(ZacavgLive(livemaxes))', 'r*', 'Color', 'red','LineWidth',4);
plot(f(recmaxes)', log10(ZacavgRec(recmaxes))', 'r*', 'Color', 'green','LineWidth',4);
hold off;

ZacliveTopo = mean(ZacpxxLive(:,freqPeaks),2);

figure; topoplot(ZacliveTopo, EEG_ZacLive.chanlocs,  'electrodes', 'labelpoint', 'chaninfo', EEG_ZacLive.chaninfo);

%% Plot topos at note frequencies
ZacrectopoSlice = plotFreqTopo([.9 2],ZacavgRec,ZacpxxRec,ZacfRec,EEG_ZacRec,"avg");
ZaclivetopoSlice = plotFreqTopo([.9 2],ZacavgLive,ZacpxxLive,ZacfLive,EEG_ZacLive,"avg");
figure; topoplot(ZaclivetopoSlice-ZacrectopoSlice, EEG_ZacLive.chanlocs,  'electrodes', 'labelpoint', 'chaninfo', EEG_ZacLive.chaninfo);  
colorbar;
clim([-6 6]);

%% Plot topos at alpha avg
ZacrectopoSlice = plotFreqTopo([8 12],ZacavgRec,ZacpxxRec,ZacfRec,EEG_ZacRec, "avg");
ZaclivetopoSlice = plotFreqTopo([8 12],ZacavgLive,ZacpxxLive,ZacfLive,EEG_ZacLive, "avg");
figure; topoplot(ZaclivetopoSlice-ZacrectopoSlice, EEG_ZacLive.chanlocs,  'electrodes', 'labelpoint', 'chaninfo', EEG_ZacLive.chaninfo);  



%% Plot topos at Gamma avg
ZacrectopoSlice = plotFreqTopo([30 44],ZacavgRec,ZacpxxRec,ZacfRec,EEG_ZacRec, "avg");
ZaclivetopoSlice = plotFreqTopo([30 44],ZacavgLive,ZacpxxLive,ZacfLive,EEG_ZacLive, "avg");
figure; topoplot(ZaclivetopoSlice-ZacrectopoSlice, EEG_ZacLive.chanlocs,  'electrodes', 'labelpoint', 'chaninfo', EEG_ZacLive.chaninfo);  


%PSDs
function [avgPSD, pxx, f] = PSD(EEG) % returns power spectral estimate V^2/Hz
    save('EEG.mat', 'EEG');
    [pxx,f] = pwelch(EEG.data',500,300,1000,EEG.srate);
    avgPSD = mean(pxx,2);
end

function [] = plotPSD(f,avg,pxx, colorAvg, colorAll)
    plot(f, log10(avg), 'Color',colorAvg,'LineWidth',4);
    plot(f, log10(pxx), 'Color',colorAll);
end

function [] = plotLiveRecorded(EEGlive, EEGrec, name)
   figure;
   hold on;
   [avgLive, pxxLive, fLive] = PSD(EEGlive);
   [avgRec, pxxRec, fRec] = PSD(EEGrec);
   plotPSD(fLive,avgLive,pxxLive,'red', [1 1 0 .3]);
   plotPSD(fRec,avgRec,pxxRec,'green', [0 0 1 .3]);
   xlim([0 58]);
   ylim([-1.5 2.5]);
   title(name);
   hold off;
end

%Wavelet
function [waveletTransform, f] = wavelet(EEG)
    [cft, f] = cwt(EEG.data(1,:), "amor", EEG.srate);
    waveletTransform = zeros([size(EEG.data,1), size(cft)]);
    for i = 1:size(EEG.data,1)
        row = EEG.data(i,:);
        cft = cwt(row, "amor", EEG.srate);
        waveletTransform(i,:,:) = cft;
    end
    waveletTransform = abs(waveletTransform);
end

function [waveavgPSD, wavePSD, f] = waveletPSD(EEG)
    [waveletTransform, f] = wavelet(EEG);
    wavePSD = squeeze(mean(waveletTransform, 3))';
    waveavgPSD = squeeze(mean(wavePSD, 2));
end

function [avgLive,avgRec,pxxLive,pxxRec,fLive,fRec] = plotWaveletPSD(EEGlive, EEGrec, name)
   figure;
   hold on;
   [avgLive, pxxLive, fLive] = waveletPSD(EEGlive);
   [avgRec, pxxRec, fRec] = waveletPSD(EEGrec);
   plotPSD(fLive,avgLive,pxxLive,'red', [1 1 0 .3]);
   plotPSD(fRec,avgRec,pxxRec,'green', [0 0 1 .3]);
   xlim([0 58]);
   ylim([0 1.5]);
   title(name);
    %plot local maxima
   livemaxes = islocalmax(avgLive);
   recmaxes = islocalmax(avgRec);
   plot(fLive(livemaxes)', log10(avgLive(livemaxes))', 'r*', 'Color', 'red','LineWidth',4);
   plot(fRec(recmaxes)', log10(avgRec(recmaxes))', 'r*', 'Color', 'green','LineWidth',4);
   hold off;
end

function [topoSlice] = plotFreqTopo(freqRange, avgPSD, pxx, frequencies, EEG, method)
    fRange = frequencies > freqRange(1) & frequencies < freqRange(2);
    if method == "maxes"
        maxes = islocalmax(avgPSD);
        freqPeaks = fRange & maxes';
    elseif method == "avg"
        freqPeaks = fRange;
    end
    topoSlice = log10(mean(pxx(:,freqPeaks),2));
    topoplot(topoSlice, EEG.chanlocs,  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo); 
    title(strcat("at Freq ", string(round(freqRange(1))),"-",string(round(freqRange(2))), " Hz"));
    % colorbar;
    % clim([-10 10]);
end

function savePSD(EEG, save_file_path, fn)
    switch fn
        case 'Pwelch'
            [avgPSD, pxx, f] = PSD(EEG);
        case 'Wavelet'
            [avgPSD, pxx, f] = waveletPSD(EEG);
    end
    matrix = [f avgPSD pxx];
    writematrix(matrix,strcat(save_file_path,'/',fn,'.csv'));
    save(strcat(save_file_path,'/',fn,'.mat'),"matrix");
end

function [result, result_all_participants] = avgParticipants(participant_table, study_path, method, condition)
    participantID = char(participant_table{1,1});
    load(strcat(study_path,'/Analysis/Participants/',participantID,'/',condition, '/',method,'.mat'),'matrix');
    result_all_participants = zeros(size(participant_table,1), size(matrix,1), size(matrix,2));
    for i = 1: size(participant_table,1)
        participantID = char(participant_table{i,1});
        load(strcat(study_path,'/Analysis/Participants/',participantID,'/',condition, '/',method,'.mat'),'matrix');
        result_all_participants(i,:,:) = matrix;
    end
    result = squeeze(mean(result_all_participants,1));
end

function iterateParticipants(participant_table, study_path, fn)
    preprocessed_data_filename = 'finished(ica).set'; 
    for i = 1: size(participant_table,1)
        participantID = char(participant_table{i,1}); 
        disp("Working on participant: " + participantID);
        for j = 1:4
            excerptNum = j;
            disp("---Excerpt " + num2str(excerptNum));
            preprocessed_data_path = strcat(study_path,'/Data/Preprocessed_Data/',participantID, '/Excerpt_',num2str(excerptNum), '/');
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
            EEG = pop_loadset('filename',preprocessed_data_filename,'filepath', preprocessed_data_path);
            switch fn
                case {'Pwelch', 'Wavelet'}
                    savePSD(EEG, save_file_path, fn);
            end 
        end
        disp("Finished participant: " + participantID);
    end
end

function plotByScore(conditionPSD, sortPart, freqRange, my_title, EEG, clim_param)
    t = tiledlayout('flow');
    frequencies = squeeze(conditionPSD(1,:,1))';
    for i = 1:21
        nexttile;
        avgPSD = squeeze(conditionPSD(sortPart(i),:,2))';
        pxx = squeeze((conditionPSD(sortPart(i),:,3:end)))';
        plotFreqTopo(freqRange, avgPSD, pxx, frequencies, EEG, 'avg')
        clim(clim_param);
        title(strcat('Participant - ', num2str(sortPart(i))));
    end
    
    t.Title.String = my_title;
    t.Title.FontWeight = 'bold';
end