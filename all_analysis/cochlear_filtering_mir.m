study_path = '/Users/arun/Library/Mobile Documents/com~apple~CloudDocs/Documents/MINDLab/performance_study_final'; %% Only thing to update!
addpath = strcat(study_path, '/Analysis/');
participant_table = readtable("Participants.xlsx");

%{
miraudio, mirfilterbank, mirenvelope, mirsum -> cochlear filterbank

miraudio, mirspectrum, mirflux
%}
%%
for i = 1: size(participant_table,1)
    participantID = char(participant_table{i,1}); 
    disp("Processing participant " + participantID);
    for j = 1:4
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
        audio = miraudio(strcat(analysis_file_path,'/normalized_audio.wav'));
        filt = mirfilterbank(audio, 'Gammatone');
        env = mirenvelope(filt, 'Sampling', 500);
        sum_bands = mirsum(env);
        filterbank_audio_data = get(sum_bands,'Data'); filterbank_audio_data = filterbank_audio_data{1}{1};
        filterbank_sr = get(sum_bands, 'Sampling'); filterbank_sr = filterbank_sr{1};
        % filterbank_audio_data_RS = resample(filterbank_audio_data,500,filterbank_sr); % resample audio to 500 hz to match EEG
        save(strcat(analysis_file_path,'/gammatone_cochlear_envelope'), "filterbank_audio_data", "filterbank_sr")
        
    end
end
