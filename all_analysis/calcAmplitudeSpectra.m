%% Spectra Gammatone Cochlear Files
clear variables;
study_path = '/Users/arun/Library/Mobile Documents/com~apple~CloudDocs/Documents/MINDLab/performance_study_final'; %% Only thing to update!
% addpath = strcat(study_path, '/Analysis/');
addpath(strcat(study_path, '/Analysis/EEGlab'))
participant_table = readtable("Participants.xlsx");
%%
spectra = zeros(4,45000,21);

%%
for i = 1:size(participant_table,1)
    for j = 1:4
        participantID = char(participant_table{i,1}); 
        disp(participantID)
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
        audio_signal = audio_signal(1:90*audio_fs); % get first 90 sec
        spectra(conditionNum,:,i) = abs(fft(audio_signal));

    
    end
end

%%
avg_spectra = squeeze(mean(spectra, 3));
%% All condition
figure;
% plot(audio_fs/length(audio_signal)*(0:length(audio_signal)-1),abs(spec),"LineWidth",1, "Color", [0.3, 0.8, 0.4], 'DisplayName', 'LiveFast');
hold on;
plot(audio_fs/length(audio_signal)*(0:length(audio_signal)-1), avg_spectra(1,:),"LineWidth",2, "Color", [0.3, 0.8, 0.4], 'DisplayName', 'LiveFast');
% plot(audio_fs/length(audio_signal)*(0:length(audio_signal)-1), avg_spectra(2,:),"LineWidth",1, "Color", [.5, 0.2, 0.2], 'DisplayName', 'LiveSlow');
plot(audio_fs/length(audio_signal)*(0:length(audio_signal)-1), avg_spectra(3,:),"LineWidth",2, "Color", [0, 0.2, 0.6], 'DisplayName', 'RecordedFast');
% plot(audio_fs/length(audio_signal)*(0:length(audio_signal)-1), avg_spectra(4,:),"LineWidth",1, "Color", [0.3, 0.2, 0.6], 'DisplayName', 'RecordedSlow');
legend;
hold off;
ylim([0, 50])
xlim([0, 50])

%% Fast v Slow

figure;
t = tiledlayout('flow')
nexttile
% plot(audio_fs/length(audio_signal)*(0:length(audio_signal)-1),abs(spec),"LineWidth",1, "Color", [0.3, 0.8, 0.4], 'DisplayName', 'LiveFast');
hold on;

% get frequencies around peaks (8.8 Hz and .84 Hz) bandwidth depends on
% frequency a bit but is approximate


hp = patch([8 8 9.6 9.6],[0 80 80 0],'k',...
    'facecolor',[.5 .5 .5],'facealpha', 0.2,'edgecolor','none', 'HandleVisibility', 'off') ;
hp2 = patch([0.6 0.6 1 1],[0 80 80 0],'k',...
    'facecolor',[.5 .5 .5],'facealpha', 0.2,'edgecolor','none', 'HandleVisibility', 'off') ;


% TRY PLOTTING LOG(AMPLITUDE)
plot(audio_fs/length(audio_signal)*(0:length(audio_signal)-1), mean(avg_spectra([1,3],:),1),"LineWidth",2, "Color", [0.5, 0.9, 0.1], 'DisplayName', 'Fast Excerpts');
plot(audio_fs/length(audio_signal)*(0:length(audio_signal)-1), mean(avg_spectra([2, 4],:),1),"LineWidth",2, "Color", [0, 0.6, 0.2], 'DisplayName', 'Slow Excerpts');
l = legend;
% l.Location = 'Southeast'
xticks(0:3:20)
% xticks(sort([xticks 0.44 0.84 2.1 4.3 8.8])); xline([0.44 0.84 2.1 4.3 7.8 8.8 9.8 ],"HandleVisibility","off","LineWidth",0.7);
% xticks(sort([xticks 7.92 9.68]))
xline(8,'-',{'8 Hz'}, 'LineWidth',1,'HandleVisibility', 'off', 'color', [.8, 0, .2], 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center','LabelOrientation', 'align');
xline(9.6,'-',{'9.6 Hz'}, 'LineWidth',1,'HandleVisibility', 'off', 'color', [.8, 0, .2], 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center','LabelOrientation', 'align');
xline(0.6,'-',{'0.6 Hz'},'LineWidth',1, 'HandleVisibility', 'off', 'color', [.8, 0, .2], 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center','LabelOrientation', 'align');
xline(1,'-',{'1 Hz'},'LineWidth',1, 'HandleVisibility', 'off', 'color', [.8, 0, .2], 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center','LabelOrientation', 'align');

fontsize(20, "points")
xlabel('Frequency [Hz]')
ylabel('Amplitude')
hold off;
ylim([0, 80])
xlim([0.2, 12.2])
set(gcf, 'color','white')
title('Amplitude Spectra of Cochlear Filtered Audio')
% exportgraphics(gcf,strcat(study_path, '/FIGURES/Spectra.png'),'Resolution',600,'BackgroundColor','none')
