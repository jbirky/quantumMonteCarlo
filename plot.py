import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
plt.style.use('classic')
rc('font', family='serif')
rc('figure', facecolor='w')
import os
import math
from math import sqrt, pi, exp
import argparse

def fill_between_steps(x, y1, y2=0, h_align='right', ax=None, **kwargs):

    if ax is None:
        ax = plt.gca()
    # First, duplicate the x values
    xx = x.repeat(2)[1:]
    # Now: the average x binwidth
    xstep = (x[1:] - x[:-1]).mean()
    # Now: add one step at end of row.
    xx = np.append(xx, xx.max() + xstep)

    # Make it possible to change step alignment.
    if h_align == 'mid':
        xx -= xstep / 2.
    elif h_align == 'right':
        xx -= xstep

    # Also, duplicate each y coordinate in both arrays
    y1 = y1.repeat(2)
    if type(y2) == np.ndarray:
        y2 = y2.repeat(2)

    # now to the plotting part:
    return ax.fill_between(xx, y1, y2=y2, **kwargs)


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Plot a file.')
	parser.add_argument("plot", action="store", type=str)
	parser.add_argument("type", action="store", type=str)
	args = parser.parse_args()

	if 'dist' in args.plot:
		folder = args.type
		file = folder + '/energies.dat'

		vals = []
		with open(file) as f:
			for line in f:
				vals.append(float(line))

		plt.hist(vals, facecolor='w')
		# plt.legend(loc='upper right')
		plt.savefig(folder + '/energies.png')
		plt.show()
		plt.close()

	if 'walkers' in args.plot:
		folder = args.type
		file = folder + '/energies.dat'

		vals = []
		with open(file) as f:
			for line in f:
				vals.append(float(line))
		steps = np.arange(0,len(vals))
		avg = round(np.mean(vals),3)
		std = round(np.std(vals)/math.sqrt(len(vals)),3)

		plt.figure(figsize=[16,5])
		plt.step(steps, vals, color='k', alpha=.7, label=r'$T=100$')
		plt.plot(steps, [avg for i in vals], color='r', label=r'$<E>=%s \pm %s$'%(str(avg), str(std)))
		plt.legend(loc='upper right')
		plt.xlim(steps[0], steps[-1])
		plt.xlabel('step', fontsize=15)
		plt.ylabel('Energy', fontsize=15)
		plt.title(folder + ' distribution', fontsize=20)
		plt.savefig(folder + '/walkers.png')
		plt.show()
		plt.close()

	if 'expected' in args.plot:
		folder = args.type
		file1 = folder + "/expected_energy.dat"
		file2 = folder + "/expected_error.dat"

		vals, errs = [], []
		with open(file1) as f:
			for line in f:
				vals.append(float(line))
		with open(file2) as f:
			for line in f:
				errs.append(float(line))

		vals, errs = np.array(vals), np.array(errs)
		steps = np.arange(0,len(vals))
		avg = round(np.mean(vals),3)

		plt.figure(figsize=[12,6])
		plt.plot(steps, vals, color='k', alpha=.7)
		plt.fill_between(steps, vals-errs, vals+errs, color='k', alpha=.3)
		# plt.legend(loc='upper right')
		plt.xlim(steps[0], steps[-1])
		plt.ylim(min(vals), max(vals))
		plt.xlabel('Tempertaure', fontsize=15)
		plt.ylabel(r'$<E>$', fontsize=15)
		plt.savefig(folder + '/expected.png')
		plt.title(folder + ' distribution', fontsize=20)
		plt.show()
		plt.close()