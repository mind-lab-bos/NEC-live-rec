%% FOOOF FOOOF FOOOF



%% Load Data
study_path = '/Users/arun/Library/Mobile Documents/com~apple~CloudDocs/Documents/performance study MINDLab/Analysis'; %% Only thing to update!
participant_table = readtable(strcat(study_path,'/Participants.xlsx'));

conditions = {'LiveFast', 'LiveSlow', 'RecordedFast', 'RecordedSlow'};
%%
%% Iterate participants and calculate fooof

for i = 1: size(participant_table,1)
        participantID = char(participant_table{i,1}); 
        disp(strcat("starting participant: ", participantID));
            for j = 1:size(conditions,2)
                condition = conditions{j};
                disp (strcat('processing condition-', condition));
                load(strcat(study_path,'/Participants/', participantID, '/', condition,'/Pwelch.mat'));
                freqs = matrix(:,1);
                avgpsd = matrix(:,2);
                allpsds = matrix(:,3:end);
                save(strcat(study_path, '/FOOOF_implement/power_spectra'), 'freqs', 'avgpsd');
                % cd( strcat(study_path, '/FOOOF_implement/'))
                % status = system('python3 fooof_imp.py');
                % disp(status);
                disp("run fooof_imp.py"); % execute: python3 fooof_imp.py
                yes = input("Did you run fooof_imp.py? (y) for yes: ", "s");
                if ~ (yes == 'y')
                    error('Please run the python fooof script');
                end
                % pyrunfile(strcat(study_path, '/FOOOF_implement/fooof_imp.py'));
                % fooof_results = load_fooof_results('fooof_results.mat');
                load('fooof_results.mat', 'aperiodic_params');
                ap = aperiodic_params(1) - log10(freqs.^aperiodic_params(2));
                corrected_avg = log10(avgpsd) - ap;
                corrected_all = log10(allpsds) - ap;
                save(strcat(study_path,'/Participants/', participantID, '/', condition,'/FOOOFED'), "corrected_all","corrected_avg");
                disp (strcat('finished condition-', condition));
            end
            disp(strcat("finished participant: ", participantID));
end


%% combine fooof data -- 1 avg fooof column and 29 channel fooof columns per participant per condition
all_fooofs = zeros(size(participant_table,1), size(conditions,2), size(freqs,1), 30); 
for i = 1: size(participant_table,1)
        participantID = char(participant_table{i,1}); 
        for j = 1:size(conditions,2)
            condition = conditions{j};
            load(strcat(study_path,'/Participants/', participantID, '/', condition,'/FOOOFED.mat'));
            all_fooofs(i,j, :, :) = [corrected_avg corrected_all];
        end
end

% keeps all participants
avglive_all = squeeze(mean(all_fooofs(:,1:2,:,:),2));
avgrec_all = squeeze(mean(all_fooofs(:,3:4,:,:),2));
avgfast_all = squeeze(mean(cat(2,all_fooofs(:,1,:,:),all_fooofs(:,3,:,:)), 2));
avgslow_all = squeeze(mean(cat(2,all_fooofs(:,2,:,:),all_fooofs(:,4,:,:)), 2));

livefast_all = squeeze(all_fooofs(:,1,:,:));
liveslow_all = squeeze(all_fooofs(:,2,:,:));
recfast_all = squeeze(all_fooofs(:,3,:,:));
recslow_all = squeeze(all_fooofs(:,4,:,:));



avglive = squeeze(mean(avglive_all,1));
avgrec = squeeze(mean(avgrec_all,1));
avgfast = squeeze(mean(avgfast_all,1));
avgslow = squeeze(mean(avgslow_all,1));

livefast = squeeze(mean(livefast_all,1));
liveslow = squeeze(mean(liveslow_all,1));
recfast = squeeze(mean(recfast_all,1));
recslow = squeeze(mean(recslow_all,1));

%% NEW PLOT ALL AVG FOOOFS All Condition

figure;
tiledlayout('flow');

nexttile;
title("Avg FOOOF all conditions")
hold on;
plot(freqs, squeeze(mean(all_fooofs(:,1,:,1),1)), 'r','LineWidth',3,'DisplayName','LiveFast')
plot(freqs, squeeze(mean(all_fooofs(:,2,:,1),1)), 'b','LineWidth',3,'DisplayName','LiveSlow')
plot(freqs, squeeze(mean(all_fooofs(:,3,:,1),1)), 'm','LineWidth',3,'DisplayName','RecordedFast')
plot(freqs, squeeze(mean(all_fooofs(:,4,:,1),1)), 'g','LineWidth',3,'DisplayName','RecordedSlow')
xlim([.2 40])
ylim([-.7 1.5])
xlabel("Frequency [Hz]")
ylabel("Relative activity [no meaningful unit]")
legend;
hold off;
set(gcf,"color","white"); fontsize("increase");



%% ALL FOOOFS


figure;
tiledlayout('flow');
nexttile;
title("Avg FOOOF all conditions")
hold on;
plot(freqs, squeeze(all_fooofs(:,1,:,1)), 'r','LineWidth',1,'DisplayName','LiveFast')
plot(freqs, squeeze(all_fooofs(:,2,:,1)), 'b','LineWidth',1,'DisplayName','LiveSlow')
plot(freqs, squeeze(all_fooofs(:,3,:,1)), 'm','LineWidth',1,'DisplayName','RecordedFast')
plot(freqs, squeeze(all_fooofs(:,4,:,1)), 'g','LineWidth',1,'DisplayName','RecordedSlow')
xlabel("Frequency [Hz]")
ylabel("Relative activity [no meaningful unit]")
xlim([.2 40])
ylim([-.7 2])
% legend;
hold off;
set(gcf,"color","white"); fontsize("increase");

%% PXX FOR TOPOS
pxx_live = avglive(:,2:end)';
pxx_rec = avgrec(:,2:end)';
pxx_fast = avgfast(:,2:end)';
pxx_slow = avgslow(:,2:end)';

pxx_livefast = livefast(:,2:end)';
pxx_liveslow = liveslow(:,2:end)';
pxx_recfast = recfast(:,2:end)';
pxx_recslow = recslow(:,2:end)';
%%
figure;
plotAllTopos(pxx_livefast,freqs,EEG, 'LiveFast');
set(gcf,"color","white"); fontsize("increase");
figure;
plotAllTopos(pxx_liveslow,freqs,EEG, 'LiveSlow');
set(gcf,"color","white"); fontsize("increase");

figure;
plotAllTopos(pxx_recfast,freqs,EEG, 'RecordedFast');
set(gcf,"color","white"); fontsize("increase");

figure;
plotAllTopos(pxx_recslow,freqs,EEG, 'RecordedSlow');
set(gcf,"color","white"); fontsize("increase");


%% Plot TOPOS live-rec and fast-slow
live_rec = pxx_live - pxx_rec;
fast_slow = pxx_fast - pxx_slow;
figure;
plotAllTopos(live_rec,freqs,EEG, 'Live-Rec');
set(gcf,"color","white"); fontsize("increase");
figure;
plotAllTopos(fast_slow,freqs,EEG, 'Fast-Slow');
set(gcf,"color","white"); fontsize("increase");

%%
figure;
plotAllTopos(pxx_fast,freqs,EEG, 'Fast');
figure;
plotAllTopos(pxx_slow,freqs,EEG, 'Slow');
%%
figure;
plotAllTopos(pxx_live,freqs,EEG, 'Live');
figure;
plotAllTopos(pxx_rec,freqs,EEG, 'Rec');


%% PSDs
% LiveFast, LiveSlow, RecFast, RecSlow
live_all = squeeze(mean(avglive_all(:,:,2:end),3));
rec_all = squeeze(mean(avgrec_all(:,:,2:end),3));
fast_all = squeeze(mean(avgfast_all(:,:,2:end),3));
slow_all = squeeze(mean(avgslow_all(:,:,2:end),3));

figure;
tiledlayout('flow');

nexttile;
title("Live(red) Rec(blue)")
hold on;
 plot(freqs, live_all,'r') 
plot(freqs, mean(pxx_live,1), 'r','LineWidth',3)


plot(freqs, rec_all,'b') 
plot(freqs, mean(pxx_rec,1), 'b','LineWidth',3)
 xlim([0 40])
ylim([-2 3])
hold off;
nexttile;

title("Fast(magenta) Slow(green)")
hold on;
plot(freqs, fast_all,'m') 
plot(freqs, mean(pxx_fast,1), 'm','LineWidth',3)
 plot(freqs, slow_all,'g') 
plot(freqs, mean(pxx_slow,1), 'g','LineWidth',3)
xlim([0 40])
ylim([-2 3])

set(gcf,"color","white")

%% PSD All condition
livefast_avgChan = squeeze(mean(livefast_all(:,:,2:end),3));
liveslow_avgChan = squeeze(mean(liveslow_all(:,:,2:end),3));
recfast_avgChan = squeeze(mean(recfast_all(:,:,2:end),3));
recslow_avgChan = squeeze(mean(recfast_all(:,:,2:end),3));

figure;
tiledlayout('flow');

nexttile;
title("All Conditions: LiveFast(red) LiveSlow(blue), RecFast(magenta) RecSlow(green)")
hold on;
plot(freqs, mean(pxx_livefast,1), 'r','LineWidth',3)
plot(freqs, mean(pxx_liveslow,1), 'b','LineWidth',3)
plot(freqs, mean(pxx_recfast,1), 'm','LineWidth',3)
plot(freqs, mean(pxx_recslow,1), 'g','LineWidth',3)
xlim([0 40])
ylim([-2 3])
hold off;

nexttile;
title("All Conditions: LiveFast(red) LiveSlow(blue), RecFast(magenta) RecSlow(green)")
hold on;
plot(freqs, livefast_avgChan,'r') 
plot(freqs, mean(pxx_livefast,1), 'r','LineWidth',3)
plot(freqs, liveslow_avgChan,'b') 
plot(freqs, mean(pxx_liveslow,1), 'b','LineWidth',3)
plot(freqs, recfast_avgChan,'m') 
plot(freqs, mean(pxx_recfast,1), 'm','LineWidth',3)
plot(freqs, recslow_avgChan,'g') 
plot(freqs, mean(pxx_recslow,1), 'g','LineWidth',3)
xlim([0 40])
ylim([-2 3])
hold off;
%% PSD PRUNE CHANNEL
channels_to_keep = [1 2 3 4 7 8 25 26 27 29];%[4 8 23];

pxx_live_pruned =pxx_live(channels_to_keep,:);
pxx_rec_pruned =pxx_rec(channels_to_keep,:);
pxx_fast_pruned =pxx_fast(channels_to_keep,:);
pxx_slow_pruned =pxx_slow(channels_to_keep,:);


pxx_livefast_pruned = pxx_livefast(channels_to_keep,:);
pxx_liveslow_pruned = pxx_liveslow(channels_to_keep,:);
pxx_recfast_pruned = pxx_recfast(channels_to_keep,:);
pxx_recslow_pruned = pxx_recslow(channels_to_keep,:);

live_all = squeeze(mean(avglive_all(:,:,channels_to_keep+1),3));
rec_all = squeeze(mean(avgrec_all(:,:,channels_to_keep+1),3));
fast_all = squeeze(mean(avgfast_all(:,:,channels_to_keep+1),3));
slow_all = squeeze(mean(avgslow_all(:,:,channels_to_keep+1),3));

figure;
tiledlayout('flow');

nexttile;
title("Live(red) Rec(blue)")
hold on;
plot(freqs, live_all,'r') 
plot(freqs, mean(pxx_live_pruned,1), 'r','LineWidth',3);

plot(freqs, rec_all,'b') 
plot(freqs, mean(pxx_rec_pruned,1), 'b','LineWidth',3);
 xlim([0 40])
ylim([-2 3]);

hold off;

nexttile;

title("Fast(magenta) Slow(green)")
hold on;
plot(freqs, fast_all,'m') 
plot(freqs, mean(pxx_fast_pruned,1), 'm','LineWidth',3)
 plot(freqs, slow_all,'g') 
plot(freqs, mean(pxx_slow_pruned,1), 'g','LineWidth',3)
xlim([0 40])
ylim([-2 3])


%% Plot diff live-rec all

%avg:
figure;
hold on;
title("difference between live and recorded FOOOFED PSDs")
plot(freqs,diff_live_rec_all_spectrum,'m') 
plot(freqs, mean(diff_live_rec_all_spectrum,1), 'm','LineWidth',3)
hold off;
xlim([0 40])
%%

function [topoSlice] = plotFreqTopo(freqRange, pxx, frequencies, EEG)
    fRange = frequencies > freqRange(1) & frequencies < freqRange(2);
    topoSlice = mean(pxx(:,fRange),2);
    topoplot(topoSlice, EEG.chanlocs,  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo); 
    title(strcat("at Freq ", string(round(freqRange(1))),"-",string(round(freqRange(2))), " Hz"));
    colorbar;
    clim([-5 5]);
end

function plotAllTopos(pxx, freqs, EEG, my_title)
    t = tiledlayout('flow');

    nexttile;
    freqRange = [8 12];
    plotFreqTopo(freqRange, pxx, freqs, EEG)
    title('8-12 Hz (Alpha)');
    clim([-2,2]);
    
    nexttile;
    freqRange = [4 8];
    plotFreqTopo(freqRange, pxx, freqs, EEG)
    title('4-8 Hz (Theta)');
    clim([-1,1]);
    
    nexttile;
    freqRange = [13 30];
    plotFreqTopo(freqRange, pxx, freqs, EEG)
    title('13-30 Hz (Beta)');
    clim([-1,1]);
    
    nexttile;
    freqRange = [30 45];
    plotFreqTopo(freqRange, pxx, freqs, EEG)
    title('30-45 Hz (Gamma)');
    clim([-1,1]);

    t.Title.String = my_title;
    t.Title.FontWeight = 'bold';
end

function new_pxx = pruneChannelsTopos(pxx, channels)
    
    zeros(size(channels,2),size(pxx,2));

end