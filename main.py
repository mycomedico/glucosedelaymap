#   A Poincare plot is a type of non-linear scatter-plot of regular time series data plotted as a function of phase
#   such that for any given point:
#   the x-coordinate is the measured value at index n and y-coordinate is the measured value at index n + 1,
#   [(n,n+1),(n+1,n+2), (n+2,n+3),...]
# Currently runs as a script that generates poincare plot from dexcom clarity csv file
# configured for mg/dL
# a color density gradient avoids overplotting in a scatterplot https://www.data-to-viz.com/graph/density2d.html
# import csv data from libreview
# plot of normal cgm data
# surface plot?  https://yuchen52.medium.com/beyond-data-scientist-3d-plots-in-python-with-examples-2a8bd7aa654b
# next i need kivymd to make a gui
# https://kivymd.readthedocs.io/en/latest/components/charts/
# user input for time delay
# user input axis of graph
# user input for x y axis start and end, next version ...selector in test.py
# overlay normal distribution or side by side ... next version

import polars as pl
from csv import reader
from datetime import datetime
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from scipy.stats import gaussian_kde


def round_to_multiple(number, multiple):
    return multiple * round(number / multiple)


def file_sel():
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(title="CGM data selector",
                                           filetypes=(("csv files", "*.csv"), ("all files", "*.*")))
    return file_path
# https://datagy.io/python-round-to-multiple/

# returns a list of tuples (timedelta since first reading,reading)

def parse_csv(file):

    with open(file, 'rt') as f:
        csv_reader = reader(f)
        for line in csv_reader:

            #dexcom
            if line[0] == 'Index':
                headers= line
                break

            #libre
            else:
                headers = next(csv_reader)
                break

        #dexcom
        if headers[0] == 'Index':
            # first 10 lines of dexcom csv are descriptive, so skip
            for line in csv_reader:
                if line[0] == '10':
                    break
            # get initial values in order to determine time delta col
            for line in csv_reader:
                datetimeinitial = datetime.strptime(line[1], '%Y-%m-%dT%H:%M:%S')
                sensorglucoseinitial = int(line[7])
                dex_data = [(0, sensorglucoseinitial)]
                break
            # time delta = datetime - datetimeinitial
            for line in csv_reader:
                readingdelta = datetime.strptime(line[1], '%Y-%m-%dT%H:%M:%S') - datetimeinitial
                readingdelta_mins = readingdelta.total_seconds() / 60
                # data needs to be normalized by rounding to nearest multiple of 5
                readingdelta_normalized = round_to_multiple(readingdelta_mins, 5)
                dex_data.append((readingdelta_normalized, int(line[7])))
            return dex_data

        else:
            for line in csv_reader:
                if line[2] and line[4]:
                    # '%m-%d-%Y %I:%M %p' us
                    # '%d-%m-%Y %H:%M'eu
                    datetimeinitial = datetime.strptime(line[2], '%m-%d-%Y %I:%M %p')
                    sensorglucoseinitial = int(line[4])
                    libre_data = [(0, sensorglucoseinitial)]
                    break
                else:
                    pass

            for line in csv_reader:
                if line[2] and line[4]:
                    readingdelta = datetime.strptime(line[2], '%m-%d-%Y %I:%M %p') - datetimeinitial
                    readingdelta_mins = readingdelta.total_seconds() / 60
                    # data needs to be normalized by rounding to nearest multiple of 5
                    readingdelta_normalized = round_to_multiple(readingdelta_mins, 5)
                    libre_data.append((readingdelta_normalized, int(line[4])))
                else:
                    pass
            return libre_data

#data is list of tuples (delta from time0,GlucoseValue)
def delay_list(data, delay=25):
    # prune list of any values that don't have value at delay time
    i = 0
    delay = round_to_multiple(delay, 5)
    offset = int(delay / 5) + 5
    obs_times = [x[0] for x in data]
    obs_times_set = set(obs_times)
    glulist = [x[1] for x in data]
    poincarelist = []

    for r, g in data:
        # if r + delay is in obs_times, then append to poincarelist
        pdelay = r + delay
        if pdelay in obs_times_set:
            delayidex = obs_times.index(pdelay, i + 1, i + offset)
            poincarelist.append((g, glulist[delayidex]))
        i += 1
    return poincarelist


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    df = pl.read_csv('pbio.2005143.s010',sep='\t', has_header=True, parse_dates=True)
    filter_df = df.filter(pl.col('subjectId') == '1636-69-035')
    sub1 = filter_df.select(['DisplayTime','GlucoseValue'])
    date0 = filter_df.select('DisplayTime').row(0)[0]
    #mins = sub1.select(
    #    pl.col('DisplayTime').apply(lambda x: round_to_multiple( (x - date0).total_seconds()/60 , 5)).alias('minute')
    #)
    sub2 = sub1.with_columns(
        (pl.col('DisplayTime') - date0).alias('dduration')
    )
    sub3 = sub2.with_columns(
        pl.col('dduration').dt.minutes().alias('dminutes')
    )
    sub4 = sub3.with_columns(
        pl.col('dminutes').apply(lambda x: round_to_multiple(x, 5)).alias('deltamins')
    )
    sub5=sub4.select([pl.col('deltamins'),pl.col('GlucoseValue')])
    sub6 = [(row[0], row[1]) for row in sub5.iter_rows()]
    sub7 = delay_list(sub6)

    x = np.asarray([x[0] for x in sub7])
    y = np.asarray([x[1] for x in sub7])
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)
    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    # color options https://matplotlib.org/2.0.2/examples/color/colormaps_reference.html
    fig, ax = plt.subplots(1,2,sharex=True, sharey=True)
    fig.suptitle('Color Density Poincaré Plot')
    bx_ = ax[0].axline((0, 0), slope=1, color='black',linestyle=':', linewidth=.5,alpha=.7)
    ax[0].set_title("delay = 25")
    ax[0].text(25,375, "non-diabetic reference")
    ax_ = ax[0].scatter(x, y, c=z, s=50, cmap=cm.jet, alpha=.3, marker = '.')
    ax[0].set_xlabel("n")
    ax[0].set_ylabel("n + delay")
    ax[0].set_xlim(0, 400)
    ax[0].set_ylim(0, 400)
    ax[0].grid(color='green', linestyle='--', linewidth=0.25)

    csvimport = file_sel()
    cgmdata = parse_csv(csvimport)

    x1 = np.asarray([x[0] for x in delay_list(cgmdata)])
    y1 = np.asarray([x[1] for x in delay_list(cgmdata)])

    xy1 = np.vstack([x1, y1])
    z1 = gaussian_kde(xy1)(xy1)

    # Sort the points by density, so that the densest points are plotted last
    idx1 = z1.argsort()
    x1, y1, z1 = x1[idx1], y1[idx1], z1[idx1]

    # color options https://matplotlib.org/2.0.2/examples/color/colormaps_reference.html
    #fig, (ax) = plt.subplots()
    bx1_ = ax[1].axline((0, 0), slope=1, color='black',linestyle=':', linewidth=.5,alpha=.7)
    ax[1].set_title("delay = 25")
    #maybe include name of file
    ax[1].text(25,375, "csv file: ")

    ax1_ = ax[1].scatter(x1, y1, c=z1, s=50, cmap=cm.jet, alpha=.3, marker = '.')
    plt.colorbar(ax1_)
    ax[1].set_xlabel("n")
    ax[1].set_ylabel("n + delay")
    ax[1].set_xlim(0, 400)
    ax[1].set_ylim(0, 400)
    ax[1].grid(color='green', linestyle='--', linewidth=0.25)
    plt.show()