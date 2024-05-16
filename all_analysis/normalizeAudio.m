study_path = '/Users/arun/Library/Mobile Documents/com~apple~CloudDocs/Documents/MINDLab/performance_study_final'; %% Only thing to update!
addpath = strcat(study_path, '/Analysis/');
conditions = {'LiveFast', 'LiveSlow', 'RecordedFast', 'RecordedSlow'};
participant_table = readtable("Participants.xlsx");
%%
normalizeaudio(participant_table,study_path);


%%
function normalizeaudio(participant_table, study_path)
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
            [audio, fs] = audioread(strcat(preprocessed_data_path,'trimmed_audio.wav'));
            [loudness, ~] = integratedLoudness(audio,fs); %% from MathWorks audio normalization help
            target = -23; % Normalize to -23 LUFS
            gaindB = target - loudness;
            gain = 10^(gaindB/20);
            xn = audio.*gain;
            audiowrite(strcat(save_file_path,'/normalized_audio.wav'),xn,fs);
        end
        disp("Finished participant: " + participantID);
    end
end