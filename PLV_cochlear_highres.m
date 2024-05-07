%% High resolution PLV cochlear
clear variables;
study_path = '/Users/arun/Library/Mobile Documents/com~apple~CloudDocs/Documents/performance study MINDLab'; %% Only thing to update!
addpath = strcat(study_path, '/Analysis/');
participant_table = readtable("Participants.xlsx");
%%
numBins = 101; % number of frequencies (frequency resolution)
lowFreq = 0.2;
highFreq = 20.2;
fois = lowFreq:(highFreq-lowFreq)/(numBins-1):highFreq; % frequency values at which we calculate PLV
cycles = 11; % determines resolution of wavelet
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
        
        % Resample audio
        EEG = pop_loadset('filename','finished(ica).set','filepath', preprocessed_data_path);
        audio_signal = load(strcat(analysis_file_path,'/cochlear_envelope.mat'),"finalSignal").finalSignal;
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
avgPLVs = squeeze(mean(aggregateChannelPLVs(:,:,:,15:17),[1 4]));
topoSlice_LiveFast = squeeze(avgPLVs(1,:,:));
topoSlice_LiveSlow = squeeze(avgPLVs(2,:,:));
topoSlice_RecFast = squeeze(avgPLVs(3,:,:));
topoSlice_RecSlow = squeeze(avgPLVs(4,:,:));
%% Plot raw topos
% eeglab redraw;
color_lim = [0 .08];
figure;
t = tiledlayout('flow');
t.Title.String = 'Phase locking Topos at 3-3.4 Hz';
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

%% Plot raw topos Live v Rec
% eeglab redraw;
color_lim = [0 .08];
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
color_lim = [0 .08];
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

aggregatePLVs = (squeeze(mean(aggregateChannelPLVs(:,:,4,:),3)));
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

%% Plot all at Fz

 aggregatePLVs = (squeeze(mean(aggregateChannelPLVs(:,:,4,:),3)));
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

