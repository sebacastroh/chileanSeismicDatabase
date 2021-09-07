import os
import pickle
import subprocess
import numpy as np
import multiprocessing
from PyPDF2 import PdfFileMerger

def compareMethodologies(n):
	subprocess.call(['python', 'compareMethodologies.py', '%i' %n])

if __name__ == '__main__':
	files = os.listdir('.')

	for filename in files:
	    if filename.endswith('.pdf'):
	        os.remove(filename)

	with open('p_waves.pkl', 'rb') as fopen:
	    p_waves = pickle.load(fopen)
	    
	valids = []
	for i,p_wave in enumerate(p_waves):
	    waveform, key, station_name, pos, status, mechanism, k = p_wave
	    if status == 'Valid' and mechanism == 'Automatic' and k != -1:
	        valids.append(i)

	nchoices = 60
	to_analyze = np.random.choice(valids, nchoices, replace=False)

	pool = multiprocessing.Pool(processes=15)
	pool.map(compareMethodologies, to_analyze)
	pool.close()

	merger = PdfFileMerger()
	for i in to_analyze:
	    merger.append('case%i.pdf' %i)

	merger.write('compareMethodologies.pdf')
	merger.close()

	for i in to_analyze:
	    os.remove('case%i.pdf' %i)
