study_path = '/Users/arun/Library/Mobile Documents/com~apple~CloudDocs/Documents/MINDLab/performance_study_final'; %% Only thing to update!
addpath = strcat(study_path, '/Analysis/');
participant_table = readtable("Participants.xlsx");

%{
miraudio, mirfilterbank, mirenvelope, mirsum -> cochlear filterbank

miraudio, mirspectrum, mirflux
%}

%%
f_0 = 1;
f_e = 50;
numBands = 10;
x_0 = log10(.88+ (f_0/165.4) )/.06; 
x_e = log10(.88+ (f_e/165.4) )/.06;
b = 1:numBands;
bandFreqs = 165.4.*(10.^(.06.*(x_0 + (b-1)*(x_e-x_0)/(numBands-1))) - .88); % Greenwood cochlear function 
% bandWidth = bandFreqs(2) - bandFreqs(1);
bandWidth = 1;
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
        [audio, fs] = audioread(strcat(analysis_file_path,'/normalized_audio.wav'));
        audio = mean(audio,2);
        audioRS = resample(audio,500,fs); % resample audio to 500 hz
        bandpassed = zeros(size(audioRS,1),numBands);
        for k = 1:numBands
            freq = bandFreqs(k);
            band = [freq - bandWidth/2, freq + bandWidth/2];
            bandpassed(:,k) = bandpass(audioRS',band,500)';
        end
        narrowBandEnvelopes = zeros(size(audioRS,1),numBands);
        for k = 1:numBands
            freq = bandFreqs(k); 
            narrowBandEnvelopes(:,k) = abs(hilbert(bandpassed(:,k)));
        end
        finalFS = 500;
        finalSignal = sum(narrowBandEnvelopes,2);
        save(strcat(analysis_file_path,'/cochlear_envelope'),"finalSignal", "finalFS", "narrowBandEnvelopes","bandpassed");
    end
end

 %{
  F(x) = 165.4(10^(.06x) - .88) -> gives the frequency at distance x (in
  millimeters) for human cochlea 
  
 let us have a lower frequency of f_0 and upper frequency of f_e with n
 bands of equal cochlear distance
 
 x_0 = log(.88+ (f_0/165.4) )/.06
 x_e = log(.88+ (f_e/165.4) )/.06

 the frequency for band #b (in range 1 to n) is given by ->
 F(x_0 + (b-1)*(x_e-x_0)/(n-1)) where F(x) is the first equation
%}
