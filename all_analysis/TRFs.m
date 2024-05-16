%% MIR Gammatone TRFs
clear variables;
study_path = '/Users/arun/Library/Mobile Documents/com~apple~CloudDocs/Documents/MINDLab/performance_study_final'; %% Only thing to update!
addpath = strcat(study_path, '/Analysis/');
participant_table = readtable("Participants.xlsx");
%% UPDATED VERSION: Compute all audio-brain TRFs for all electrodes and all participants
% IF already computed, load flux_and_envelope_TRF.mat
all_envelope_TRFs = nan(21,4,29,65);
all_flux_TRFs = nan(21,4,29,65);
for i = 1: size(participant_table,1)
    for j = 1:4
        participantID = char(participant_table{i,1});
        disp(participantID);
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
        
        % Load EEG struct
        EEG = pop_loadset('filename','finished(ica).set','filepath', preprocessed_data_path);
        eegdata = EEG.data';
        % load amplitude env signal
        amplitude_env = load(strcat(analysis_file_path,'/gammatone_cochlear_envelope.mat'),"filterbank_audio_data").filterbank_audio_data;
        amplitude_envSrate = 500;
        % load spectral flux env signal
        fileID = fopen(strcat(analysis_file_path,'/mir_onset_env.txt'));
        flux_env = fscanf(fileID,'%f');
        fclose(fileID);
        flux_envSrate = 196;
        % Z-score audio feature data
        flux_env = normalize(flux_env);
        amplitude_env = normalize(amplitude_env);
        % Downsample all signals
        downsampleFreq = 128;
        flux_env = resample(flux_env,downsampleFreq,flux_envSrate);
        amplitude_env = resample(amplitude_env,downsampleFreq,amplitude_envSrate);
        eegdata = resample(eegdata,downsampleFreq,EEG.srate);
        
        signallengths = [length(eegdata), length(amplitude_env), length(flux_env)];
        minsignallength = min(signallengths);
        % error if signal lengths are more than 4 frames apart
        if any((signallengths - minsignallength) > 4)
            error("signal sizes don't match")
        end
        % Trim data if needed
        if (length(eegdata) > minsignallength)
            eegdata = eegdata(1:minsignallength,:);
        end
        if (length(amplitude_env) > minsignallength)
            amplitude_env = amplitude_env(1:minsignallength);
        end
        if (length(flux_env) > minsignallength)
            flux_env = flux_env(1:minsignallength);
        end
       % Calc TRFs for all electrodes
       [all_envelope_TRFs(i,conditionNum,:,:), time_lags] = calcTRF(amplitude_env,eegdata, downsampleFreq);
       all_flux_TRFs(i,conditionNum,:,:) = calcTRF(flux_env,eegdata, downsampleFreq);
    end
end
save("mir_flux_and_envelope_TRF","all_envelope_TRFs","all_flux_TRFs","time_lags")


%% Plot Stuff
% Flux v Env corr across time
electrodes = [4 8 26];

envLiveFast_corrlags = squeeze(mean(all_envelope_TRFs(:,1,electrodes,:),[1 3]));
envLiveSlow_corrlags = squeeze(mean(all_envelope_TRFs(:,2,electrodes,:),[1 3]));
envRecFast_corrlags = squeeze(mean(all_envelope_TRFs(:,3,electrodes,:),[1 3]));
envRecSlow_corrlags = squeeze(mean(all_envelope_TRFs(:,4,electrodes,:),[1 3]));

fluxLiveFast_corrlags = squeeze(mean(all_flux_TRFs(:,1,electrodes,:),[1 3]));
fluxLiveSlow_corrlags = squeeze(mean(all_flux_TRFs(:,2,electrodes,:),[1 3]));
fluxRecFast_corrlags = squeeze(mean(all_flux_TRFs(:,3,electrodes,:),[1 3]));
fluxRecSlow_corrlags = squeeze(mean(all_flux_TRFs(:,4,electrodes,:),[1 3]));

% errors
err_envLiveFast_corrlags = std(squeeze(mean(all_envelope_TRFs(:,1,electrodes,:),3)))/sqrt(size(all_envelope_TRFs,1));
err_envLiveSlow_corrlags = std(squeeze(mean(all_envelope_TRFs(:,2,electrodes,:),3)))/sqrt(size(all_envelope_TRFs,1));
err_envRecFast_corrlags = std(squeeze(mean(all_envelope_TRFs(:,3,electrodes,:),3)))/sqrt(size(all_envelope_TRFs,1));
err_envRecSlow_corrlags = std(squeeze(mean(all_envelope_TRFs(:,4,electrodes,:),3)))/sqrt(size(all_envelope_TRFs,1));

err_fluxLiveFast_corrlags = std(squeeze(mean(all_flux_TRFs(:,1,electrodes,:),3)))/sqrt(size(all_flux_TRFs,1));
err_fluxLiveSlow_corrlags = std(squeeze(mean(all_flux_TRFs(:,2,electrodes,:),3)))/sqrt(size(all_flux_TRFs,1));
err_fluxRecFast_corrlags = std(squeeze(mean(all_flux_TRFs(:,3,electrodes,:),3)))/sqrt(size(all_flux_TRFs,1));
err_fluxRecSlow_corrlags = std(squeeze(mean(all_flux_TRFs(:,4,electrodes,:),3)))/sqrt(size(all_flux_TRFs,1));

% plot envelope
figure;
% tiledlayout('flow')
% nexttile;
hold on;
plot(time_lags,envLiveFast_corrlags,'DisplayName','LiveFast','Color',[1 0 0],'LineWidth',3)
errorbar(time_lags,envLiveFast_corrlags,err_envLiveFast_corrlags,'HandleVisibility','off','Color',[1 0 0 .5],'CapSize',0)

plot(time_lags,envLiveSlow_corrlags,'DisplayName','LiveSlow','Color',[.6 .1 .3],'LineWidth',3)
errorbar(time_lags,envLiveSlow_corrlags,err_envLiveSlow_corrlags,'HandleVisibility','off','Color',[.6 .1 .3 .5],'CapSize',0)

plot(time_lags,envRecFast_corrlags,'DisplayName','RecordedFast','Color',[0 0 1],'LineWidth',3)
errorbar(time_lags,envRecFast_corrlags,err_envRecFast_corrlags,'HandleVisibility','off','Color',[0 0 1 .5],'CapSize',0)

plot(time_lags,envRecSlow_corrlags,'DisplayName','RecordedSlow','Color',[0 .2 .6],'LineWidth',3)
errorbar(time_lags,envRecSlow_corrlags,err_envRecSlow_corrlags,'HandleVisibility','off','Color',[0 .2 .6 .5],'CapSize',0)
yline(0,'--','LineWidth',3,'Color','black','HandleVisibility','off')
ylim([-.02 .04])

xlabel('Time lag (ms)')
ylabel('TRF corr')
title("Envelope TRF Correlations")
legend
set(gcf,'color','white')
fontsize(30,'points')
hold off;

% plot Flux
figure;
hold on;
plot(time_lags,fluxLiveFast_corrlags,'DisplayName','LiveFast','Color',[1 0 0],'LineWidth',3)
errorbar(time_lags,fluxLiveFast_corrlags,err_fluxLiveFast_corrlags,'HandleVisibility','off','Color',[1 0 0 .5],'CapSize',0)

plot(time_lags,fluxLiveSlow_corrlags,'DisplayName','LiveSlow','Color',[.6 .1 .3],'LineWidth',3)
errorbar(time_lags,fluxLiveSlow_corrlags,err_fluxLiveSlow_corrlags,'HandleVisibility','off','Color',[.6 .1 .3 .5],'CapSize',0)

plot(time_lags,fluxRecFast_corrlags,'DisplayName','RecordedFast','Color',[0 0 1],'LineWidth',3)
errorbar(time_lags,fluxRecFast_corrlags,err_fluxRecFast_corrlags,'HandleVisibility','off','Color',[0 0 1 .5],'CapSize',0)

plot(time_lags,fluxRecSlow_corrlags,'DisplayName','RecordedSlow','Color',[0 .2 .6],'LineWidth',3)
errorbar(time_lags,fluxRecSlow_corrlags,err_fluxRecSlow_corrlags,'HandleVisibility','off','Color',[0 .2 .6 .5],'CapSize',0)
yline(0,'--','LineWidth',3,'Color','black','HandleVisibility','off')
ylim([-.02 .04])

xlabel('Time lag (ms)')
ylabel('TRF corr')
title("Flux TRF Correlations")
legend
set(gcf,'color','white')
fontsize(30,'points')
hold off;
%% Live Vs Rec Fast Vs Slow
% flux
fluxLive_corrlags = squeeze(mean(all_flux_TRFs(:,[1 2],electrodes,:),[1,  2, 3]));
err_fluxLive_corrlags = std(squeeze(mean(all_flux_TRFs(:,[1 2],electrodes,:),[2, 3])))/sqrt(size(all_flux_TRFs,1));
fluxRec_corrlags = squeeze(mean(all_flux_TRFs(:,[3 4],electrodes,:),[1,  2, 3]));
err_fluxRec_corrlags = std(squeeze(mean(all_flux_TRFs(:,[3 4],electrodes,:),[2, 3])))/sqrt(size(all_flux_TRFs,1));

fluxFast_corrlags = squeeze(mean(all_flux_TRFs(:,[1 3],electrodes,:),[1,  2, 3]));
err_fluxFast_corrlags = std(squeeze(mean(all_flux_TRFs(:,[1 3],electrodes,:),[2, 3])))/sqrt(size(all_flux_TRFs,1));
fluxSlow_corrlags = squeeze(mean(all_flux_TRFs(:,[2 4],electrodes,:),[1,  2, 3]));
err_fluxSlow_corrlags = std(squeeze(mean(all_flux_TRFs(:,[2 4],electrodes,:),[2, 3])))/sqrt(size(all_flux_TRFs,1));

% plot Flux Live/Rec

figure;
tiledlayout('flow')
nexttile;
hold on;
plot(time_lags,fluxLive_corrlags,'DisplayName','Live','Color',[1 .4 0],'LineWidth',3)
errorbar(time_lags,fluxLive_corrlags,err_fluxLive_corrlags,'HandleVisibility','off','Color',[1 .4 0 .5],'CapSize',0)
yline(mean(fluxLive_corrlags),'--','LineWidth',1,'Color',[1 .4 0 .5],'HandleVisibility','off')

plot(time_lags,fluxRec_corrlags,'DisplayName','Recorded','Color',[0 .4 1],'LineWidth',3)
errorbar(time_lags,fluxRec_corrlags,err_fluxRec_corrlags,'HandleVisibility','off','Color',[0 .4 1 .5],'CapSize',0)
yline(mean(fluxRec_corrlags),'--','LineWidth',1,'Color',[0 .4 1 0.5],'HandleVisibility','off')

yline(0,'--','LineWidth',3,'Color','black','HandleVisibility','off')
ylim([-.02 .04])
xlabel('Time lag (ms)')
ylabel('TRF corr')
title("Flux Live/Rec TRF Correlations")
legend
set(gcf,'color','white')
fontsize(30,'points')
hold off;

% plot Flux Fast/Slow
% figure;
nexttile;
hold on;
plot(time_lags,fluxFast_corrlags,'DisplayName','Fast','Color',[.4 .4 0],'LineWidth',3)
errorbar(time_lags,fluxFast_corrlags,err_fluxFast_corrlags,'HandleVisibility','off','Color',[.4 .4 0 .5],'CapSize',0)
yline(mean(fluxFast_corrlags),'--','LineWidth',1,'Color',[.4 .4 0 .5],'HandleVisibility','off')

plot(time_lags,fluxSlow_corrlags,'DisplayName','Slow','Color',[0 .4 .4],'LineWidth',3)
errorbar(time_lags,fluxSlow_corrlags,err_fluxSlow_corrlags,'HandleVisibility','off','Color',[0 .4 .4 .5],'CapSize',0)
yline(mean(fluxSlow_corrlags),'--','LineWidth',1,'Color',[0 .4 .4 .5],'HandleVisibility','off')

yline(0,'--','LineWidth',3,'Color','black','HandleVisibility','off')
ylim([-.02 .04])
xlabel('Time lag (ms)')
ylabel('TRF corr')
title("Flux Live/Rec TRF Correlations")
legend
set(gcf,'color','white')
fontsize(30,'points')
hold off;
%% Flux vs envelope
flux =  squeeze(mean(all_flux_TRFs,[2 3]));
errflux = std(flux)/sqrt(size(flux,1));
envelope = squeeze(mean(all_envelope_TRFs,[2 3]));
errenv = std(envelope)/sqrt(size(envelope,1));

% plot Flux vs envelope
figure;
hold on;
% plot(time_lags,mean(flux,1),'DisplayName','Flux','Color',[1 0 .6],'LineWidth',3)
% errorbar(time_lags,mean(flux,1),errflux,'HandleVisibility','off','Color',[1 0 .6 .5],'CapSize',0)

shadedErrorBar(time_lags,mean(flux,1),errflux,'lineprops',{'Color',[1 0 .6],'LineWidth',3,'DisplayName',"Spectral Flux"})
% plot(time_lags,mean(envelope,1),'DisplayName','Envelope','Color',[.6 .8 0],'LineWidth',3)
% errorbar(time_lags,mean(envelope,1),errenv,'HandleVisibility','off','Color',[.6 .8 0 .5],'CapSize',0)

shadedErrorBar(time_lags,mean(envelope,1),errenv,'lineprops',{'Color',[0.6 0.8 0],'LineWidth',3,'DisplayName',"Amplitude Envelope"})


yline(0,'--','LineWidth',3,'Color','black','LineWidth', 1,'HandleVisibility','off')
xline(164.0625,'label','t = 164.06 ms','LineWidth', 1,'color',[0.8, 0.1, 0.2],'HandleVisibility','off')
ylim([-.02 .04])
xlabel('Time lag (ms)')
ylabel('TRF corr')
title("Spectral Flux vs. Amplitude Envelope TRF Correlations")
legend
set(gcf,'color','white')
fontsize(30,'points')
hold off;

%% T test for above plot
time_lag_index = 22;
[h, p, ci, stats] = ttest(flux(:,time_lag_index), envelope(:,time_lag_index));


%% TOPOS
%% plot Flux Topos by condition
% times = find(time_lags < 250); 
times = 22; %22 for peak

liveFast_topo = squeeze(mean(all_flux_TRFs(:,1,:,times),[1 4]));
liveSlow_topo = squeeze(mean(all_flux_TRFs(:,2,:,times),[1 4]));
recFast_topo = squeeze(mean(all_flux_TRFs(:,3,:,times),[1 4]));
recSlow_topo = squeeze(mean(all_flux_TRFs(:,4,:,times),[1 4]));
figure;
t = tiledlayout('flow');
nexttile;
topoplot(liveFast_topo, EEG.chanlocs,  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo); 
clim([-.05 .05])
title("Live Fast");
nexttile;
topoplot(liveSlow_topo, EEG.chanlocs,  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo); 
clim([-.05 .05])
title("Live Slow");
nexttile;
topoplot(recFast_topo, EEG.chanlocs,  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo); 
clim([-.05 .05])
title("Recorded Fast");
nexttile;
topoplot(recSlow_topo, EEG.chanlocs,  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo); 
clim([-.05 .05])
title("Recorded Slow");
c = colorbar;
c.Layout.Tile = 'East';
t.Title.String = 'Spectral Flux TRF Correlations at Optimal Time Lag';
t.Title.FontWeight = 'bold';
set(gcf, 'color', 'white')
 fontsize(20,'points')
%%
%% plot Env Topos by condition
liveFast_topo = squeeze(mean(all_envelope_TRFs(:,1,:,:),[1 4]));
liveSlow_topo = squeeze(mean(all_envelope_TRFs(:,2,:,:),[1 4]));
recFast_topo = squeeze(mean(all_envelope_TRFs(:,3,:,:),[1 4]));
recSlow_topo = squeeze(mean(all_envelope_TRFs(:,4,:,:),[1 4]));
figure;
t = tiledlayout('flow');
nexttile;
topoplot(liveFast_topo, EEG.chanlocs,  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo); 
clim([-.03 .03])
title("Live Fast");
nexttile;
topoplot(liveSlow_topo, EEG.chanlocs,  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo); 
clim([-.03 .03])
title("Live Slow");
nexttile;
topoplot(recFast_topo, EEG.chanlocs,  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo); 
clim([-.03 .03])
title("Recorded Fast");
nexttile;
topoplot(recSlow_topo, EEG.chanlocs,  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo); 
clim([-.03 .03])
title("Recorded Slow TRF corr");
c = colorbar;
c.Layout.Tile = 'East';
t.Title.String = 'Envelope TRF Correlation Topos by Condition';
t.Title.FontWeight = 'bold';
fontsize(20,'points')
%% plot means of conditions for Flux
electrodes = [8 4 26];
times = 22;
liveFast = squeeze(mean(all_flux_TRFs(:,1,electrodes,times),[3, 4]));
liveSlow = squeeze(mean(all_flux_TRFs(:,2,electrodes,times),[3, 4]));
recFast = squeeze(mean(all_flux_TRFs(:,3,electrodes,times),[3, 4]));
recSlow = squeeze(mean(all_flux_TRFs(:,4,electrodes,times),[3, 4]));
eLiveFast = std(liveFast)/sqrt(length(liveFast));
eLiveSlow = std(liveSlow)/sqrt(length(liveSlow));
eRecFast = std(recFast)/sqrt(length(recFast));
eRecSlow = std(recSlow)/sqrt(length(recSlow));
figure;
t = tiledlayout('flow');
t.Title.String = 'Mean TRF Corr by Condition'
t.Title.FontWeight = 'bold'
nexttile;
hold on;
x = ["LiveFast", "LiveSlow", "RecFast", "RecSlow"];
data = [mean(liveFast), mean(liveSlow), mean(recFast), mean(recSlow)]';
bar(x,data);
ylim([-0.03,0.05])
ylabel("TRF Corr")
er = errorbar([1 2 3 4],data,[eLiveFast, eLiveSlow, eRecFast, eRecSlow]);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off;

nexttile;
hold on;
live = squeeze(mean(all_flux_TRFs(:,[1 2], electrodes, times),[2,3,4]));
recorded = squeeze(mean(all_flux_TRFs(:,[3 4], electrodes, times),[2,3,4]));
elive = std(live)/sqrt(length(live));
erecorded = std(recorded)/sqrt(length(recorded));
x = ["Live", "Recorded"];
data = [mean(live), mean(recorded)]';
b = bar(x,data,'FaceColor','flat');
b.CData(1,:) = [1 .4 0];
b.CData(2,:) = [0 .4 1];
ylim([-0.03,0.05])
ylabel("TRF Corr")
er = errorbar([1 2],data,[elive, erecorded]);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off;

nexttile;
hold on;
fast = squeeze(mean(all_flux_TRFs(:,[1 3], electrodes, times),[2,3,4]));
slow = squeeze(mean(all_flux_TRFs(:,[2 4], electrodes, times),[2,3,4]));
efast = std(fast)/sqrt(length(fast));
eslow = std(slow)/sqrt(length(slow));
x = ["Fast", "Slow"];
data = [mean(fast), mean(slow)]';
b = bar(x,data,'FaceColor','flat');
b.CData(1,:) = [.4 .4 0];
b.CData(2,:) = [0 .4 .4];
ylim([-0.03,0.05])
ylabel("TRF Corr")
er = errorbar([1 2],data,[elive, erecorded]);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off;

fontsize(30,'points')
%%
% Inputs: audio_stim is time varying audio feature, EEGdata is time varying eeg
% signal at all electrodes (matrix: time points by electrode)
function [correlationByElectrode, t] = calcTRF(audio_stim, EEGdata, fs)
    stim = audio_stim; 
    %Calc TRFs for every electrode
    correlationByElectrode = zeros(29,65);
    for k = 1:29 
        resp = EEGdata(:,k);
        % Z-score EEG data
        resp = normalize(resp);
        % ---Cross-validation---
        % Generate training/test sets
        nfold = 6;
        testfold = 1;
        [stimtrain,resptrain,stimtest,resptest] = mTRFpartition(stim,resp,nfold,...
            testfold);
        % Model hyperparameters
        tmin = 0;
        tmax = 500;
        lambdas = 10.^(-6:2:6);
        
        % Run fast cross-validation (to get optimal ridge parameter) with
        % FORWARD MODEL using all-lags to optimize
        [cv, t] = mTRFcrossval(stimtrain,resptrain,fs,1,tmin,tmax,lambdas,...
            'zeropad',0,'fast',1,'type','multi');
        
        % ---Model training---
        
        % Get optimal hyperparameters
        [~,idx] = max(mean(cv.r));
        lambda = lambdas(idx);
        % nlambda = length(lambdas);
        
        % Train FORWARD model for every 'single' time lag
        model = mTRFtrain(stimtrain,resptrain,fs,1,tmin,tmax,lambda,'zeropad',0, 'type', 'single');
        
        % ---Model testing---
        
        % Test model
        [~,test] = mTRFpredict(stimtest,resptest,model,'zeropad',0);
        correlationByElectrode(k,:) = squeeze(test.r);
    end
end
    
