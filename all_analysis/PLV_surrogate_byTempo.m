%% Gammatone PLV cochlear
clear variables;
study_path = '/Users/arun/Library/Mobile Documents/com~apple~CloudDocs/Documents/MINDLab/performance_study_final'; %% Only thing to update!
% addpath = strcat(study_path, '/Analysis/');
addpath(strcat(study_path, '/Analysis/EEGlab'))
participant_table = readtable("Participants.xlsx");
eeglab; close;
%% define bootstrap random picks
% conditions = {'LiveFast', 'LiveSlow', 'RecordedFast', 'RecordedSlow'};

surrogate_size = 100; %number of random pairings

randomSeed = rng(0); %if you want reproduceable results, seed with 0


random_participantNumsFast = randi([1 21], surrogate_size, 1);
random_audio_excerptNumsFast = randi([1, 2], surrogate_size,1);
random_EEG_excerptNumsFast = randi([1, 2], surrogate_size,1);
random_EEG_channelsFast = randi([1, 29], surrogate_size, 1);


random_participantNumsSlow = randi([1 21], surrogate_size, 1);
random_audio_excerptNumsSlow = randi([1, 2], surrogate_size,1);
random_EEG_excerptNumsSlow = randi([1, 2], surrogate_size,1);
random_EEG_channelsSlow = randi([1, 29], surrogate_size, 1);

%% Parameters for PLV calculation
numBins = 101; % number of frequencies (frequency resolution)
lowFreq = 0.2;
highFreq = 20.2;
fois = lowFreq:(highFreq-lowFreq)/(numBins-1):highFreq; % frequency values at which we calculate PLV
cycles = 5; % determines resolution of wavelet

surrogatePLVs = zeros(surrogate_size, 2, numBins); % two is for Fast and Slow conditions

%% calc all random PLVs
conditions = {'LiveFast', 'LiveSlow', 'RecordedFast', 'RecordedSlow'};
fastConditions = {"LF", "RF"};
slowConditions = {"LS", "RS"};


for pair_index = 1:surrogate_size
    disp("working on pair index: " + num2str(pair_index))
    for condition = [1, 2]
        switch condition
            case 1
                participantNum = random_participantNumsFast(pair_index);
                participantID = char(participant_table{participantNum,1});
                participantExcerpts = {participant_table.First{participantNum}, participant_table.Second{participantNum}, participant_table.Third{participantNum}, participant_table.Fourth{participantNum}};
                    
                EEG_excerptNum = sum(times(1:4, (participantExcerpts == fastConditions{random_EEG_excerptNumsFast(pair_index)}))); % get the excerpt (arg == excerptCondition)
          
                mapConditions = [1, 3];
                audio_condition = conditions{mapConditions(random_audio_excerptNumsFast(pair_index))};
                eegChannel = random_EEG_channelsFast(pair_index);
           
            case 2
                participantNum = random_participantNumsSlow(pair_index);                
                participantID = char(participant_table{random_participantNumsSlow(pair_index),1});
                participantExcerpts = {participant_table.First{participantNum}, participant_table.Second{participantNum}, participant_table.Third{participantNum}, participant_table.Fourth{participantNum}};
                EEG_excerptNum = sum(times(1:4, (participantExcerpts == slowConditions{random_EEG_excerptNumsSlow(pair_index)}))); % get the excerpt (arg == excerptCondition)
                mapConditions = [2, 4];
                audio_condition = conditions{mapConditions(random_audio_excerptNumsSlow(pair_index))};
                eegChannel = random_EEG_channelsSlow(pair_index);
        end
    
        PLVs = zeros(numBins,1); % initialize PLV placeholder 
    
        % participantID = char(participant_table{random_participantNums(pair_index),1}); 
        
        disp("Processing participant: " + participantID)
    
        % load EEG data
        % EEG_excerptNum = random_EEG_excerptNums(pair_index);
       
        preprocessed_data_path = strcat(study_path,'/Data/Preprocessed_Data/',participantID, '/Excerpt_',num2str(EEG_excerptNum), '/');
        EEG = pop_loadset('filename','finished(ica).set','filepath', preprocessed_data_path);
    
        disp("working on audio excerpt: " + audio_condition)
        % load Audio data
        % audio_condition = conditions{random_audio_excerptNums(pair_index)};
        analysis_file_path = strcat(study_path,'/Analysis/Participants/',participantID,'/',audio_condition);
        audio_signal = load(strcat(analysis_file_path,'/gammatone_cochlear_envelope.mat'),"filterbank_audio_data").filterbank_audio_data;
        audio_fs = 500;
    
        % get EEG data and trim
        % eegChannel = random_EEG_channels(pair_index);
        eegdata = EEG.data(eegChannel,:)';
        if (size(eegdata,1) > size(audio_signal,1))
            eegdata = eegdata(1:size(audio_signal));
        else
            audio_signal = audio_signal(1:size(eegdata, 1));
        end
    
        % Calc PLV
        for l = 1:numBins
            cmw = designMorlet(fois(l), cycles, EEG.srate); % build a wavelet with N cycles
            [EEG_phase, EEG_amplitude, EEG_filtered, EEG_power] = waveletConvolution(eegdata,cmw); % convolve with EEG
            [Audio_phase, Audio_amplitude, Audio_filtered, Audio_power] = waveletConvolution(audio_signal,cmw); % convolve with Audio
            PLVs(l) = phaseCoherence(EEG_phase, Audio_phase);
        end
        surrogatePLVs(pair_index,condition,:) = PLVs;
    end

end
%% Calculate mean and deviation

averageFast = mean(squeeze(surrogatePLVs(:,1,:)),1);
deviationFast = std(squeeze(surrogatePLVs(:,1,:)),1);

averageSlow = mean(squeeze(surrogatePLVs(:,2,:)),1);
deviationSlow = std(squeeze(surrogatePLVs(:,2,:)),1);
%% Z-score aggregate-plv

zScorePLV = zeros(size(aggregatePLVs));
all_zScorePLV = zeros(size(aggregateChannelPLVs));

for binNum = 1:numBins
    % Calc Fast
    zScorePLV(:,[1,3],binNum) = (aggregatePLVs(:,[1,3],binNum) - averageFast(binNum))/deviationFast(binNum);
    all_zScorePLV(:,[1,3],:,binNum) = (aggregateChannelPLVs(:,[1,3],:,binNum) - averageFast(binNum))/deviationFast(binNum);
    % Calc Slow
    zScorePLV(:,[2,4],binNum) = (aggregatePLVs(:,[2,4],binNum) - averageSlow(binNum))/deviationSlow(binNum);
    all_zScorePLV(:,[2,4],:,binNum) = (aggregateChannelPLVs(:,[2,4],:,binNum) - averageSlow(binNum))/deviationSlow(binNum);

end


%% PLOTTING


%% Shaded errors ALL
PLV_LiveFast = squeeze(mean(zScorePLV(:,1,:),1));
PLV_LiveSlow = squeeze(mean(zScorePLV(:,2,:),1));
PLV_RecFast = squeeze(mean(zScorePLV(:,3,:),1));
PLV_RecSlow = squeeze(mean(zScorePLV(:,4,:),1));
eLiveFast = calcErrors(squeeze(zScorePLV(:,1,:)));
eLiveSlow = calcErrors(squeeze(zScorePLV(:,2,:)));
eRecFast = calcErrors(squeeze(zScorePLV(:,3,:)));
eRecSlow = calcErrors(squeeze(zScorePLV(:,4,:)));


figure;
% tiledlayout('flow')
nexttile;
hold on;

hp = patch([8 8 9.6 9.6],[0 .14 .14 0],'k',...
    'facecolor',[.5 .5 .5],'facealpha', 0.2,'edgecolor','none', 'HandleVisibility', 'off') ;
hp2 = patch([0.6 0.6 1 1],[0 .14 .14 0],'k',...
    'facecolor',[.5 .5 .5],'facealpha', 0.2,'edgecolor','none', 'HandleVisibility', 'off') ;

shadedErrorBar(fois,PLV_LiveFast,eLiveFast,'lineprops',{'Color',[.8 0 0],'LineWidth',4,'DisplayName',"Live Fast"})
% shadedErrorBar(fois,PLV_LiveSlow,eLiveSlow,'lineprops',{'Color',[.7 .4 0],'LineWidth',4,'DisplayName',"Live Slow"})
shadedErrorBar(fois,PLV_RecFast,eRecFast,'lineprops',{'Color',[0 .6 0.4],'LineWidth',4,'DisplayName',"Recorded Fast"})
% shadedErrorBar(fois,PLV_RecSlow,eRecSlow,'lineprops',{'Color',[0 0 .8],'LineWidth',4,'DisplayName',"Recorded Slow"})

xline(8,'-',{'8 Hz'}, 'LineWidth',1, 'HandleVisibility', 'off', 'color', [.8, 0, .2], 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center','LabelOrientation', 'align');
xline(9.6,'-',{'9.6 Hz'}, 'LineWidth',1,'HandleVisibility', 'off', 'color', [.8, 0, .2], 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center','LabelOrientation', 'align');
xline(0.6,'-',{'0.6 Hz'},'LineWidth',1, 'HandleVisibility', 'off', 'color', [.8, 0, .2], 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center','LabelOrientation', 'align');
xline(1,'-',{'1 Hz'},'LineWidth',1, 'HandleVisibility', 'off', 'color', [.8, 0, .2], 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center','LabelOrientation', 'align');

legend;
hold off;
 title('PLVs By Condition')
xlabel('Frequency [Hz]')
ylabel('Phase Coherence')
xlim([.2 20.2]);
% ylim([0 0.14])
set(gcf,'color','w');
fontsize(25, "points")

%% BY DIFF Q
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
            ylim([-1, 2])
            % yticks(0:.02:.14)
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
            ylim([-1, 2])
            % yticks(0:.02:.14)
            % xticks(sort([xticks .8])); xline([.8],"HandleVisibility","off","LineWidth",1);
            xlabel("Frequency [Hz]"); ylabel("PLV");
    
        case 4
            colorHigh = [0, .2, 0.8];
            colorLow = [0, .2, 0.4];
            label = "Recorded,";
    end
    disp(condition)
    lowPLV = squeeze(mean(all_zScorePLV(participant_ranking(1:10,1),condition,electrodes,:),[2 3]));
    highPLV = squeeze(mean(all_zScorePLV(participant_ranking(12:21,1),condition,electrodes,:),[2 3]));
    
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
    % shadedErrorBar(fois,mean(lowPLV),errlow,'lineprops',{'Color',colorLow,'LineWidth',2,'DisplayName',label + " Low Rating Diff"})
    legend;
end


hold off;
xlabel("Frequency [Hz]")
ylabel("PLV")

% legend;
fontsize(20,'points')
set(gcf,'color','white')









%% Topo
topoSlice = squeeze(mean(aggregateChannelPLVs(1,:,:,5),[1 2]));
figure;
topoplot(topoSlice, EEG.chanlocs,  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo); 
%%
avgPLVs = squeeze(mean(all_zScorePLV(:,:,:,40:48),[1 4]));
topoSlice_LiveFast = squeeze(avgPLVs(1,:,:));
topoSlice_LiveSlow = squeeze(avgPLVs(2,:,:));
topoSlice_RecFast = squeeze(avgPLVs(3,:,:));
topoSlice_RecSlow = squeeze(avgPLVs(4,:,:));
%% Plot raw topos
% eeglab redraw;
color_lim = [-2 2];
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



%%

function errors = calcErrors(dataMatrix) %column is participats and row is data points
    errors = zeros(1, size(dataMatrix,2));
    for i = 1: size(dataMatrix,2)
        errors(i) = std(dataMatrix(:,i))/sqrt(size(dataMatrix,1));
    end
end

%%

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


