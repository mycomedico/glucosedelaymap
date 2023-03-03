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
#option to show pdf color bar
import polars as pl
import re
from csv import reader
from datetime import datetime
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from scipy.stats import gaussian_kde
from ntpath import basename

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
'''
deprecated in favor of polars
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
#needs to be replaced by polar joins
def delay_map(data, delay=25):
    # prune list of any values that don't have value at delay time
    delay = round_to_multiple(delay, 5)
    offset = int(delay / 5) + 5
    obs_times = [x[0] for x in data]
    obs_times_set = set(obs_times)
    glulist = [x[1] for x in data]
    poincarelist = []
    for i, (r, g) in enumerate(data):
        # if r + delay is in obs_times, then append to poincarelist
        pdelay = r + delay
        if pdelay in obs_times_set:
            delayidex = obs_times.index(pdelay, i + 1, i + offset)
            poincarelist.append((g, glulist[delayidex]))
    return poincarelist
'''

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    delays = input('Enter delay values you would like to plot (values must be a multiple of 5 and comma seperated, e.g. 5,15,25):')
    delays_list = delays.split(',')
    a = len(delays_list)
    control = input('Should a reference plot of non-diabetic control data be shown? (y/n): ')
    if control == ('y' or 'Y' or 'yes' or 'Yes' or 'YES'):
        control = True
        b = 2
    else:
        b = 1
        control = False
    #total number of plots    
    numplots = a * b
    pattern = 'Index.*'
    numplotslist = list(range(a * b))
    x = numplotslist.copy()
    y = numplotslist.copy()
    z = numplotslist.copy()
    xy = numplotslist.copy()
    idx = numplotslist.copy()
    # a = number of delay points, b= 1 for data alone, 2 w control
    fig, ax = plt.subplots(a,b,sharex=True, sharey=True)
    fig.suptitle('Color Density Poincaré Plot')
    if control:
        #every other plot is control
        controlplots = numplotslist[0::2]
        patientplots = numplotslist[1::2]

        # control data public cgm data from hall et al.  https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.2005143
        df = pl.read_csv('pbio.2005143.s010', sep='\t', has_header=True, parse_dates=True)
        #control data
        date0 = df.select('DisplayTime').row(0)[0]
        df = df.with_columns(
            (pl.col('DisplayTime') - date0).alias('duration')
        )
        df = df.with_columns(
            pl.col('duration').dt.minutes().alias('dminutes')
        )
        df = df.with_columns(
            pl.col('dminutes').apply(lambda x: round_to_multiple(x, 5)).alias('deltamins')
        )
        df = df.filter(pl.col('subjectId') == '1636-69-035')
        df = df.select([pl.col('deltamins'), pl.col('GlucoseValue')])
        df = df.unique(subset='deltamins', keep='first')

        #control data plots
        for index, plot in enumerate(controlplots):
            if len(controlplots) == 1:
                subp = ax[0]
            else:
                subp = ax[index, 0]

            delay = int(delays_list[index])
            dd = df.select('deltamins') + delay
            djoiny = df.join(dd, on='deltamins', how='inner')
            djoiny = djoiny.with_columns(
                pl.col('GlucoseValue').alias('GlucoseValue2')
            )
            # now in order to get the x coordinate, we need to subtract 15 from the y coordinate
            de = djoiny.select('deltamins') - delay
            djoinx = df.join(de, on='deltamins', how='inner')
            djoinx = djoinx.select(['GlucoseValue'])
            djoiny = djoiny.select(['GlucoseValue2'])
            djoinxy = pl.concat(
                [
                    djoinx,
                    djoiny,
                ],
                how='horizontal'
            )
            # color options https://matplotlib.org/2.0.2/examples/color/colormaps_reference.html
            x[plot] = djoinxy.get_column('GlucoseValue').to_numpy()
            y[plot] = djoinxy.get_column('GlucoseValue2').to_numpy()
            xy[plot] = np.vstack([x[plot], y[plot]])
            z[plot] = gaussian_kde(xy[plot])(xy[plot])
            idx[plot] = z[plot].argsort()
            x[plot], y[plot], z[plot] = x[plot][idx[plot]], y[plot][idx[plot]], z[plot][idx[plot]]
            if index == 0:
                subp.set_title("normal control")
            subp.set_xlabel("n")
            subp.set_ylabel("n + " + str(delays_list[index]) + " mins")
            subp.grid(color='green', linestyle='--', linewidth=0.25)
            subp.scatter(x[plot], y[plot], c=z[plot], s=50, cmap=cm.jet, alpha=.3, marker='.')
            subp.axline((0, 0), slope=1, color='black', linestyle=':', linewidth=.5, alpha=.7)
            #ax[plot].set_xlim(0, 400)
            #ax[plot].set_ylim(0, 400)

#this works for libre but need to manually input usa or euro, maybe query user for country??
#need to update for dexcom
        #load pt data
        csvimport = file_sel()
        di = pl.read_csv(csvimport,has_header=False)
        cgmtype = di.select('column_1').row(0)[0]
        #pattern = 'Index.*'


        if re.match(pattern, cgmtype):

            if cgmtype == 'Index':
                #this is a clarity file so....
                df = pl.read_csv(csvimport,has_header=True,try_parse_dates=True,null_values=['Low','High','NotComputable','SensorNotActive','SensorWarmup','OutOfRange','NoSensor','InvalidReading','SensorFailed','SensorInitializing','SensorCalibration'])

            else:
                df = pl.read_csv(csvimport, has_header=True, try_parse_dates=True, sep='\t',
                                 null_values=['Low', 'High', 'NotComputable', 'SensorNotActive', 'SensorWarmup',
                                              'OutOfRange', 'NoSensor', 'InvalidReading', 'SensorFailed',
                                              'SensorInitializing', 'SensorCalibration'])


            df = df.lazy().drop_nulls('Transmitter ID').collect()
            df = df.with_columns(
                (pl.col('Glucose Value (mg/dL)').cast(pl.Int64)).alias('glu'))
            df = df.lazy().drop_nulls('glu').collect()
            df = df.sort('Timestamp (YYYY-MM-DDThh:mm:ss)')
            date0 = df.select('Timestamp (YYYY-MM-DDThh:mm:ss)').row(0)[0]
            df = df.with_columns(
                (pl.col('Timestamp (YYYY-MM-DDThh:mm:ss)') - date0).dt.minutes().alias('duration')
            )

        else:
            #this is a libre file so....
            # we need a regex to detect whether or not 'PM or AM' is in the date column
            df = pl.read_csv(csvimport,has_header=True,try_parse_dates=True,skip_rows=1)
            datetest= df.select('Device Timestamp').row(0)[0]
            patternl = '.*[AP]M.*'

            if re.match(patternl, datetest):
                df = df.with_columns(
                    (pl.col('Device Timestamp').str.strptime(pl.Datetime, fmt='%m-%d-%Y %I:%M %p')).alias('timestamp'))
                #print('usa version')
            else:
                #print('euro version')
                df = df.with_columns(
                    (pl.col('Device Timestamp').str.strptime(pl.Datetime, fmt='%d-%m-%Y %H:%M')).alias('timestamp'))

            df = df.with_columns(
                (pl.col('Scan Glucose mg/dL').cast(pl.Int64).fill_null(0) + pl.col('Historic Glucose mg/dL')).alias('glu'))
            df = df.lazy().drop_nulls('glu').collect()
            df = df.sort('timestamp')
            date0 = df.select('timestamp').row(0)[0]
            df = df.with_columns(
                (pl.col('timestamp') - date0).dt.minutes().alias('duration')
            )
        df = df.with_columns(
            pl.col('duration').apply(lambda x: round_to_multiple(x, 5)).alias('deltamins')
        )
        df = df.select(['deltamins', 'glu'])
        #remove duplicates
        df = df.unique(subset='deltamins', keep='first')

        #patient data plots
        for index, plot in enumerate(patientplots):
            if len(patientplots) == 1:
                subpt = ax[1]
            else:
                subpt = ax[index, 1]

            delay = int(delays_list[index])

            # djoiny gives y coordinate in glu col, returned table is - deltamins, glu
            dd = df.select('deltamins') + delay
            djoiny = df.join(dd, on='deltamins', how='inner')
            djoiny = djoiny.with_columns(
                pl.col('glu').alias('glu2')
            )
            # now in order to get the x coordinate, we need to subtract 15 from the y coordinate
            de = djoiny.select('deltamins') - delay
            djoinx = df.join(de, on='deltamins', how='inner')
            djoinx = djoinx.select(['glu'])
            djoiny = djoiny.select(['glu2'])
            djoinxy = pl.concat(
                [
                    djoinx,
                    djoiny,
                ],
                how='horizontal'
            )
            if(len(djoinxy) <5):
                print("not enough data to plot")
                exit()
            # color options https://matplotlib.org/2.0.2/examples/color/colormaps_reference.html
            x[plot] = djoinxy.get_column('glu').to_numpy()
            y[plot] = djoinxy.get_column('glu2').to_numpy()
            xy[plot] = np.vstack([x[plot], y[plot]])
            z[plot] = gaussian_kde(xy[plot])(xy[plot])
            # Sort the points by density, so that the densest points are plotted last
            idx[plot] = z[plot].argsort()
            x[plot], y[plot], z[plot] = x[plot][idx[plot]], y[plot][idx[plot]], z[plot][idx[plot]]
            if index == 0:
                subpt.set_title("file: " + basename(csvimport),size=8)
            subpt.set_xlabel("n")
            subpt.set_ylabel("n + " + str(delays_list[index]) + " mins")
            subpt.grid(color='green', linestyle='--', linewidth=0.25)
            subpt.scatter(x[plot], y[plot], c=z[plot], s=50, cmap=cm.jet, alpha=.3, marker='.')
            subpt.axline((0, 0), slope=1, color='black', linestyle=':', linewidth=.5, alpha=.7)

    else:
    #patient data plots
        #fig, ax = plt.subplots(numplots,sharex=True, sharey=True)
        patientplots = numplotslist
        # load pt data
        csvimport = file_sel()
        di = pl.read_csv(csvimport, has_header=False)
        cgmtype = di.select('column_1').row(0)[0]

        if re.match(pattern, cgmtype):
            # this is a clarity file so....
            if cgmtype == 'Index':
                df = pl.read_csv(csvimport, has_header=True, try_parse_dates=True,
                                 null_values=['Low', 'High', 'NotComputable', 'SensorNotActive', 'SensorWarmup',
                                              'OutOfRange', 'NoSensor', 'InvalidReading', 'SensorFailed',
                                              'SensorInitializing', 'SensorCalibration'])

            else:
                df = pl.read_csv(csvimport, has_header=True, try_parse_dates=True, sep='\t',
                                 null_values=['Low', 'High', 'NotComputable', 'SensorNotActive', 'SensorWarmup',
                                              'OutOfRange', 'NoSensor', 'InvalidReading', 'SensorFailed',
                                              'SensorInitializing', 'SensorCalibration'])


            df = df.lazy().drop_nulls('Transmitter ID').collect()
            df = df.with_columns(
                (pl.col('Glucose Value (mg/dL)').cast(pl.Int64)).alias('glu'))
            df = df.lazy().drop_nulls('glu').collect()
            df = df.sort('Timestamp (YYYY-MM-DDThh:mm:ss)')
            date0 = df.select('Timestamp (YYYY-MM-DDThh:mm:ss)').row(0)[0]
            df = df.with_columns(
                (pl.col('Timestamp (YYYY-MM-DDThh:mm:ss)') - date0).dt.minutes().alias('duration')
            )

        else:
            # this is a libre file so....
            df = pl.read_csv(csvimport, has_header=True, try_parse_dates=True, skip_rows=1)
            datetest = df.select('Device Timestamp').row(0)[0]
            patternl = '.*[AP]M.*'

            if re.match(patternl, datetest):
                df = df.with_columns(
                    (pl.col('Device Timestamp').str.strptime(pl.Datetime, fmt='%m-%d-%Y %I:%M %p')).alias('timestamp'))
                # print('usa version')
            else:
                # print('euro version')
                df = df.with_columns(
                    (pl.col('Device Timestamp').str.strptime(pl.Datetime, fmt='%d-%m-%Y %H:%M')).alias('timestamp'))

            #df = pl.read_csv(csvimport, has_header=True, try_parse_dates=True, skip_rows=1)
            df = df.with_columns(
                (pl.col('Scan Glucose mg/dL').cast(pl.Int64).fill_null(0) + pl.col('Historic Glucose mg/dL')).alias('glu'))
            df = df.lazy().drop_nulls('glu').collect()
            df = df.sort('timestamp')
            date0 = df.select('timestamp').row(0)[0]
            df = df.with_columns(
                (pl.col('timestamp') - date0).dt.minutes().alias('duration')
            )

        df = df.with_columns(
            pl.col('duration').apply(lambda x: round_to_multiple(x, 5)).alias('deltamins')
        )
        df = df.select(['deltamins', 'glu'])
        # remove duplicates
        df = df.unique(subset='deltamins', keep='first')

        # patient data plots
        for index, plot in enumerate(patientplots):
            if len(patientplots) == 1:
                subpt = ax
            else:
                subpt = ax[index]
            delay = int(delays_list[index])

            # djoiny gives y coordinate in glu col, returned table is - deltamins, glu
            dd = df.select('deltamins') + delay
            djoiny = df.join(dd, on='deltamins', how='inner')
            djoiny = djoiny.with_columns(
                pl.col('glu').alias('glu2')
            )
            # now in order to get the x coordinate, we need to subtract 15 from the y coordinate
            de = djoiny.select('deltamins') - delay
            djoinx = df.join(de, on='deltamins', how='inner')
            djoinx = djoinx.select(['glu'])
            djoiny = djoiny.select(['glu2'])
            djoinxy = pl.concat(
                [
                    djoinx,
                    djoiny,
                ],
                how='horizontal'
            )
            if(len(djoinxy) <5):
                print("not enough data to plot")
                exit()
            # color options https://matplotlib.org/2.0.2/examples/color/colormaps_reference.html
            x[plot] = djoinxy.get_column('glu').to_numpy()
            y[plot] = djoinxy.get_column('glu2').to_numpy()
            xy[plot] = np.vstack([x[plot], y[plot]])
            z[plot] = gaussian_kde(xy[plot])(xy[plot])
            # Sort the points by density, so that the densest points are plotted last
            idx[plot] = z[plot].argsort()
            x[plot], y[plot], z[plot] = x[plot][idx[plot]], y[plot][idx[plot]], z[plot][idx[plot]]
            if index == 0:
                subpt.set_title("file: " + basename(csvimport), size=8)
            subpt.set_xlabel("n")
            subpt.set_ylabel("n + " + str(delays_list[index]) + " mins")
            subpt.grid(color='green', linestyle='--', linewidth=0.25)
            subpt.scatter(x[plot], y[plot], c=z[plot], s=50, cmap=cm.jet, alpha=.3, marker='.')
            subpt.axline((0, 0), slope=1, color='black', linestyle=':', linewidth=.5, alpha=.7)
    plt.show()

