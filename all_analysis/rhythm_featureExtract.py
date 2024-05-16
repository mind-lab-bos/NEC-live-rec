import librosa
import soundfile as sf
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

filename = '/Users/Arun/Downloads/tomorrow.wav'
datapath = '/Users/arun/Library/Mobile Documents/com~apple~CloudDocs/Documents/performance study MINDLab/Analysis/Participants/'
participantfile = '/Users/arun/Library/Mobile Documents/com~apple~CloudDocs/Documents/performance study MINDLab/Analysis/participants.xlsx'
condition = 'RecordedSlow'
writeFileName = 'saveOnset.txt'
print(filename)

df = pd.read_excel(participantfile)
participantArray = df['Participant'].values

participantID = participantArray[1]
soundFileName = datapath+participantID+'/' + condition + '/normalized_audio.wav'
y, sr = librosa.load(soundFileName)
sampleClickLength = 250 # length of click in samples
time_units = np.array(list(range(sampleClickLength)));

click = np.sin(2*np.pi*440/sr*time_units)

y = np.append(y, np.zeros(sampleClickLength*2)) ## zero padding of the audio at the end -- length is twice the click length
onsets = librosa.onset.onset_detect(y=y, sr=sr, normalize = 'True', units = 'time')
onset_env = librosa.onset.onset_strength(y=y, sr=sr)
pulse = librosa.beat.plp(onset_envelope=onset_env, sr=sr) #local pulse
beats_plp = np.flatnonzero(librosa.util.localmax(pulse)) #local beats
times = librosa.times_like(onset_env, sr=sr) # gives us the time value at each frame
beatTimes = times[beats_plp]

print(len(times), len(onset_env))
insertedOnsetSignal = y
for onsetTime in onsets:
	if onset_env[np.where(times==onsetTime)[0][0]] > 1.5:
		np.put(insertedOnsetSignal,list(range(int(sr*onsetTime), int(sr*onsetTime)+sampleClickLength)), click)
sf.write('testAudioOnset.wav',insertedOnsetSignal,sr)

# insertedBeatSignal = y
# for beat in beatTimes:
#  	np.put(insertedBeatSignal,list(range(int(sr*beat), int(sr*beat)+sampleClickLength)), click)
# sf.write('testAudioBeat.wav',insertedBeatSignal,sr)

# print("check")
onsets = librosa.onset.onset_detect(y=y, sr=sr, normalize = 'True', units = 'time')
# onset_env = librosa.onset.onset_strength(y=y, sr=sr)

# pulse = librosa.beat.plp(onset_envelope=onset_env, sr=sr) #local pulse
# beats_plp = np.flatnonzero(librosa.util.localmax(pulse)) #local beats
# times = librosa.times_like(onset_env, sr=sr) # gives us the time value at each frame


# print("check2")
# times = librosa.times_like(onset_env, sr=sr)
# onsets2 = librosa.onset.onset_detect(onset_envelope = onset_env, sr=sr, normalize = 'True', units = 'time')
# print("check3")

plt.plot(times,onset_env,label = 'Onset strength')
plt.vlines(onsets, 0, onset_env.max(), color='r', alpha=0.9,
           linestyle='--', label='Onsets')
# plt.vlines(times[beats_plp], 0, onset_env.max(), alpha=0.5, color='b', linestyle='dotted', label='PLP Beats')
plt.xlim(0,20)
plt.show()
# beatTimes = times[beats_plp]

# fileBeatTime = open('beatTime.txt','w+')
# fileBeatTime.write('BeatTimes: \n')
# for beattime in beatTimes:
# 	fileBeatTime.write(str(beattime) + '\n')

# fileBeatTime.write('OnsetTimes: \n')
# for onsetTime in onsets:
# 	fileBeatTime.write(str(onsetTime) + '\n')
# fileBeatTime.close()




# file = open(writeFileName, "w+")
# for participantID in participantArray:
# 	print('Processing ' + participantID)
# 	soundFileName = datapath+participantID+'/' + condition + '/normalized_audio.wav'
# 	y, sr = librosa.load(soundFileName)
# 	onsets = librosa.onset.onset_detect(y=y, sr=sr, normalize = 'True')
# 	onset_times = np.array(librosa.frames_to_time(onsets))
# 	diff_onset = np.diff(onset_times)
# 	medianDuration = np.median(diff_onset)
# 	# avgDuration = np.mean(diff_onset)
# 	file.write(str(medianDuration)+ '\n')
# file.close()

# fileread = open(writeFileName,'r')
# lines = fileread.read().split('\n')
# lines = lines[0:-1]
# print(lines)
# fileread.close()
# finalarray = []
# for num in lines:
# 	finalarray.append(float(num))
# print(np.mean(finalarray))

# print(diff_onset)
# plt.hist(diff_onset, bins = 100, color = 'skyblue')
# plt.xlabel('Time')
# plt.ylabel('Frequency')
# plt.show()

# avgDuration = np.mean(diff_onset)


