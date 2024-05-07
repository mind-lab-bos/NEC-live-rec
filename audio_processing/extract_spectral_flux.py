import librosa
import soundfile as sf
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

datapath = '/Users/arun/Library/Mobile Documents/com~apple~CloudDocs/Documents/MINDLab/performance_study_final/Analysis/Participants/'
# analysispath = '/Users/arun/Library/Mobile Documents/com~apple~CloudDocs/Documents/performance study MINDLab/Analysis/'
participantfile = '/Users/arun/Library/Mobile Documents/com~apple~CloudDocs/Documents/MINDLab/performance_study_final/Analysis/participants.xlsx'
conditions = ['LiveFast','LiveSlow','RecordedFast','RecordedSlow']
hop_length = 225 # the sampling rate of the extracted feature signal will be sr/hop_length ** IMPORTANT
df = pd.read_excel(participantfile)
participantArray = df['Participant'].values
for participantID in participantArray:
	print("Processing participant: " + participantID)
	for condition in conditions:
		print("processing condition: " + condition)
		soundFileName = datapath+participantID+'/' + condition + '/normalized_audio.wav'
		y, sr = librosa.load(soundFileName, sr = 44100) # we use 44100 as sr
		onset_env = librosa.onset.onset_strength(y=y ,sr=sr,hop_length = hop_length) # calc spectral flux
		times = librosa.times_like(onset_env, sr=sr, hop_length= hop_length) # gives us the time value at each index of onset_env
		timediffs = np.round(np.diff(times), decimals = 8)
		finaltimediff = 0.00510204 # the onset strength function should have entries with this time step if final sr is 196
		print(len(times), len(onset_env), len(y)//sr//finaltimediff) # these 3 should be the same
		print(timediffs == timediffs[0]) # all time diffs should match the first
		if finaltimediff==timediffs[0] and np.all(timediffs == timediffs[0]):
			print("correct")
		else: 
			raise Exception("time diffs don't match " + str(finaltimediff) + \
				    "\n Your time diff was " + str(timediffs[0]))
		np.savetxt(datapath+participantID+'/'+condition+"/mir_onset_env.txt",onset_env)
print("SUCCESS")



