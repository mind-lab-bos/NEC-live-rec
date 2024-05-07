%% Updated Preprocessing
%%
clear;
close all;
%%% Update these for each participant
participantID = '230822_CLUE';
excerptNum = '4';
audio_filename = 'Excerpt4 - CLUE.mov';
%%%

%%% EEG Processing Parameters
low_threshhold = .2;
high_threshhold = 58;
refChannel = 'A2';
chanRange = 1:30;
%%%

%%% FilePath Parameters
study_path = '/Users/arun/Library/Mobile Documents/com~apple~CloudDocs/Documents/MINDLab/performance_study_final'; %% Only thing to update!
cd(study_path)
addpath('./EEGlab'); % add eeglab to path
audio_file_path = strcat(study_path,'/Data/Raw_Data/',participantID,'/Recordings/',audio_filename);
raw_data_path = strcat(study_path,'/Data/Raw_Data/',participantID,'/EEG/',participantID,'/');
raw_data_filename = strcat(participantID,'_',excerptNum,'.vhdr'); %.vhdr file
save_dataset_path = strcat(study_path,'/Data/Preprocessed_Data/',participantID,'/Excerpt_',excerptNum,'/');
mkdir(save_dataset_path);
chanlocs_filepath = strcat(study_path,'/Analysis/EEGlab/plugins/dipfit/standard_BEM/elec/standard_1005.elc');
%%%

checkFiles(audio_file_path,strcat(raw_data_path,raw_data_filename),chanlocs_filepath); %% checks that files are valid

%% Audio processing
[audio, audio_fs] = audioread(audio_file_path);
%% Inspection of Beginning of Audio
[scope, param] = inspectAudio(audio, audio_fs, "beginning");
scope(param); % you should zoom in and scroll within the scope for accuracy

% DETERMINE ONSET OF 2ND TO LAST BEEP AND VERIFY THAT ONSETS OF LAST TWO BEEPS ARE 3 SECONDS APART
% DETERMINE START OF VIOLIN PLAYING
%% WRITE 2nd to last BEEP and MUSIC START IN SECONDS
beepstartTime = 16.616; % INPUT THE TIME OF 2ND TO LAST BEEP -- BE EXTREMELY PRECISE
startOfExcerpt = 34.955; % INPUT THE TIME OF MUSIC STARTING
%% -- test beep sound --
clear sound;
playSound(audio,audio_fs,beepstartTime,"beginning"); %YOU SHOULD HEAR LAST TWO BEEPS THEN SILENCE 
%% -- test audio start --
clear sound;
playSound(audio,audio_fs,startOfExcerpt,"beginning"); %YOU SHOULD HEAR MUSIC START RIGHT AWAY
%% Inspection of End of Audio
clear sound;
[scope, param] = inspectAudio(audio, audio_fs, "end");
scope(param);
%% WRITE END OF EXCERPT IN SECONDS
endOfExcerpt = 238.08; % INPUT THE END TIME OF EXCERPT
%% -- test audio end --
clear sound;
playSound(audio,audio_fs,endOfExcerpt,"end"); %YOU SHOULD HEAR THEN END OF EXCERPT WHEN PLAYBACK STOPS
%% Save Trimmed Audio
[final_audio, audio_fs, lengthOfaudio, timeBeforeMusic] = saveTrimmedAudio(audio,audio_fs, beepstartTime, startOfExcerpt, endOfExcerpt, save_dataset_path);
%{
GREAT! You've completed the audio pre-processing for the excerpt. 
Now follow the next steps for EEG preprocessing.
%}

%% Pre-ICA Processing
clear scope;
eeglab;
[EEG, triggerDiffs] = pre_trim(raw_data_path,raw_data_filename);
close;
openvar('triggerDiffs');
%% WRITE FIRST BEEP TRIGGER EVENT NUMBER
% DETERMINE EVENT NUMBER OF 1ST BEEP's TRIGGER FROM 'triggerDiffs'
% LATENCY DIFFS SHOULD FOLLOW THIS PATTERN FROM THE FIRST TRIGGER: 2000 1500 1000 1500
        % ---- not needed but here is length of beeps and silences in
        % -----seconds: 2 —- 2 - 1 - 1 —- 2 (longer dashes are two seconds)
firstBeepTriggerNumber = 1;  % will be an odd number

%Calculates start time of music by adding: start time of first BEEP trigger, length of beep track, length of wait before music
musicStartTime = EEG.event(firstBeepTriggerNumber).latency/EEG.srate + 14 + timeBeforeMusic; 
%% Pre-ICA Processing
[EEG, preRejEEG] = process_pre_ICA(EEG, musicStartTime, lengthOfaudio, chanlocs_filepath, chanRange, low_threshhold, high_threshhold, refChannel, save_dataset_path);
%% plotting
pop_eegplot(EEG, 1, 1, 1);
% pop_eegplot(preRejEEG, 1, 1, 1);
figure; topoplot([],EEG.chanlocs, 'style', 'blanEEk',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
% plot channel spectra and maps for psd rejection
% eeglab redraw; % if additional channels need to be removed
%% Interpolate More Channels
more_channels_to_remove = {'PO7', 'FC5'}; % 'AF8', 'F8', 'T8', 'C4', 'FC6', 'F7', 'AF7'
EEG = remove_more_chans(EEG, more_channels_to_remove);
%% ICA
close all;
EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on');
EEG = pop_iclabel(EEG, 'default');
% eeglab redraw;
pop_viewprops(EEG, 0, 1:28); %% you may need to use the eeglab gui -- tools/Classify Components using ICLabel/Label Components
%% deleting components
deleteComponent = [1 2 3 10 11 13 14 18 20 22 23 24];
EEG = pop_subcomp( EEG, deleteComponent, 0);

%save final preprocessing
EEG = pop_saveset( EEG, 'filename','finished(ica).set','filepath',save_dataset_path);
save(strcat(save_dataset_path,'PreprocessedWorkspace'), '-regexp', '^(?!(param|audio|final_audio)$).') % saves all variables except param, audio and final_audio as a workspace file
pop_eegplot(EEG, 1, 1, 1);

% pop_topoplot(EEG, 0, 1:29 ,EEG.setname,[5 6] ,0,'electrodes','on');
%%
function [EEG, preRejEEG] = process_pre_ICA(EEG, sampleStart, sampleLength, chanlocs_filepath, chanRange, low_threshhold, high_threshhold, refChannel, save_dataset_path)
    time_range = [sampleStart sampleStart+sampleLength].*EEG.srate;
    fprintf('time range of points in seconds: %f\n',time_range./EEG.srate);
    EEG = pop_select(EEG,'point', time_range);
    disp('Completed EEG time trimming');
    EEG = pop_chanedit(EEG, 'lookup', chanlocs_filepath); 
    disp('Completed EEG channel location mapping');
    EEG = pop_select( EEG, 'channel',chanRange);
    EEG = pop_saveset( EEG, 'filename','trimmed.set','filepath',save_dataset_path);

    disp('Completed EEG trimming');
    EEG = pop_reref(EEG, refChannel);
    disp('Completed EEG rereferencing');
    EEG = pop_saveset( EEG, 'filename','reref.set','filepath',save_dataset_path);

    EEG = pop_eegfiltnew(EEG, 'locutoff',low_threshhold,'hicutoff',high_threshhold);
    disp('Completed EEG filtering');
    EEG = pop_saveset( EEG, 'filename','filter.set','filepath',save_dataset_path);

    preRejEEG = EEG;
    EEG = clean_artifacts(EEG,'ChannelCriterion',0.5,'LineNoiseCriterion',4,'BurstCriterion',10,'WindowCriterion','off','BurstRejection','off', 'ChannelCriterionMaxBadTime', .5);
    disp('Completed EEG artifact cleaning');
    EEG = pop_saveset( EEG, 'filename','cleaned.set','filepath',save_dataset_path);

    EEG = pop_interp(EEG, preRejEEG.chanlocs, 'spherical');
    disp('Completed EEG interpolation');
    EEG = pop_saveset( EEG, 'filename','interp.set','filepath',save_dataset_path);
end

function EEG = remove_more_chans(EEG, chansToRemove)
    preRejEEG = EEG;
    EEG = pop_select( EEG, 'rmchannel',chansToRemove);
    EEG = pop_interp(EEG, preRejEEG.chanlocs, 'spherical');
    disp('Removed and interpolated back channels'); disp(chansToRemove);
end

function [EEG, triggerDiffs] = pre_trim(raw_data_path,raw_data_filename)
    EEG = pop_loadbv(raw_data_path, raw_data_filename);
    triggerDiffs = processTriggers(EEG);
end

function triggerDiffs = processTriggers(EEG)
    triggerDiffs(1).num = 0;
    triggerDiffs(1).latency = 0;
    for i = 1:size(EEG.event,2)
        if isequal(EEG.event(i).type,'2048') %EEG.event(i).type,'0') || 
            triggerDiffs(end+1).num = EEG.event(i).bvmknum;
            triggerDiffs(end).latency = EEG.event(i).latency;
            if (i ~= 1)
                triggerDiffs(end).diff = triggerDiffs(end).latency - triggerDiffs(end-1).latency;
            end
        end
    end
    triggerDiffs = triggerDiffs(1, 2:size(triggerDiffs,2));
end

function [scope, param] = inspectAudio(audio, audio_fs, location)
    scope = timescope(SampleRate= audio_fs, YLimits = [-1 1], TimeSpanSource = "auto");
    scope.PeakFinder.Enabled = true;
    scope.PeakFinder.LabelPeaks = true;
    scope.CursorMeasurements.Enabled = true;
    if location == "beginning"
        length = 50; %% shows first 50 seconds of the audio recording
        param = audio(1:audio_fs*length,:); 
    elseif location == "end"
        param = audio; %% shows all of the audio recording
    end
end

function playSound(audio, audio_fs, time, location)
    soundLength = 10;
    if location == "beginning"
        sound(audio(time*audio_fs:(time+soundLength)*audio_fs,:),audio_fs);
    elseif location == "end"
        sound(audio((time-soundLength)*audio_fs:time*audio_fs,:),audio_fs);
    end
end

function[final_audio, audio_fs, lengthOfaudio, timeBeforeMusic] = saveTrimmedAudio(audio,audio_fs, beepstartTime, startOfExcerpt, endOfExcerpt, save_dataset_path)
    endOfBeeps = 5 + beepstartTime; % end of beeps is 5 seconds after the onset of the second to last beep
    timeBeforeMusic = startOfExcerpt - endOfBeeps;
    final_audio = audio(startOfExcerpt*audio_fs:endOfExcerpt*audio_fs,:);
    lengthOfaudio = endOfExcerpt - startOfExcerpt;
    audiowrite(strcat(save_dataset_path,'trimmed_audio.wav'),final_audio,audio_fs);
    disp("Trimmed audio file saved!");
end

function checkFiles(varargin)
    for i = 1:nargin
        if ~isfile(varargin{i})
            error(varargin{i} + " does not locate a valid file");
        end
    end
    disp("Yay all files and paths valid and loaded!");
end

