#   A Poincare plot is a type of non-linear scatter-plot of regular time series data plotted as a function of phase
#   such that for any given point:
#   the x-coordinate is the measured value at index n and y-coordinate is the measured value at index n + 1,
#   [(n,n+1),(n+1,n+2), (n+2,n+3),...]
# Currently runs as a script that generates poincare plot from dexcom clarity csv file
# configured for mg/dL
# a color density gradient avoids overplotting in a scatterplot https://www.data-to-viz.com/graph/density2d.html

# surface plot?  https://yuchen52.medium.com/beyond-data-scientist-3d-plots-in-python-with-examples-2a8bd7aa654b
# https://kivymd.readthedocs.io/en/latest/components/charts/
# user input for x y axis start and end, next version ...selector in test.py
# overlay normal distribution or side by side ... next version
#Specify date range of analysis period... Start Date and Time, End Date and Time
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

    color_legend = input('Display color density function legend for each plot? (y/n): ')
    if color_legend == ('y' or 'Y' or 'yes' or 'Yes' or 'YES'):
        color_legend = True
    else:
        color_legend = False

    daterange = input('Would you like to select a range of dates to analyze? (n = use all data): ')
    if daterange == ('y' or 'Y' or 'yes' or 'Yes' or 'YES'):
        daterange = True
    else:
        daterange = False

    axes_scale = input('Auto-detect axes scale? (y/n): ')
    if axes_scale == ('n' or 'no' or 'NO' or 'No'):
        x_origin = int(input('Enter x-axis origin (e.g. 0) : '))
        xdomain = int(input('Enter x-axis domain (e.g. 300) : '))
        y_origin = int(input('Enter y-axis origin (e.g. 0) : '))
        yrange = int(input('Enter y-axis range (e.g. 300) : '))
        axes_scale = True
    else:
        axes_scale = False

    #total number of plots
    numplots = a * b
    pattern = 'Index.*'
    patternl = '.*[AP]M.*'
    numplotslist = list(range(numplots))
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
            con_plot = subp.scatter(x[plot], y[plot], c=z[plot], s=50, cmap=cm.jet, alpha=.3, marker='.')
            subp.axline((0, 0), slope=1, color='black', linestyle=':', linewidth=.5, alpha=.7)
            plt.colorbar(con_plot)


#this works for libre but need to manually input usa or euro, maybe query user for country??
#need to update for dexcom
        #load pt data
        csvimport = file_sel()
        di = pl.read_csv(csvimport,has_header=False)
        cgmtype = di.select('column_1').row(0)[0]

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
            df = df.with_columns(
                (pl.col('Timestamp (YYYY-MM-DDThh:mm:ss)')).alias('timestamp')
            )
        else:
            #this is a libre file so....
            # we need a regex to detect whether or not 'PM or AM' is in the date column
            df = pl.read_csv(csvimport,has_header=True,try_parse_dates=True,skip_rows=1)
            datetest= df.select('Device Timestamp').row(0)[0]

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
        if daterange:
            dstart = df.head(1).select('timestamp')
            dend = df.tail(1).select('timestamp')
            print('Date Range of uploaded file is from ' + str(dstart[0,0]) + ' to ' + str(dend[0,0]))
            analysis_start = input('Enter analysis start date/time (following format YYYY-MM-DD HH:mm): ')
            analysis_end = input('Enter analysis end date/time (same format): ')
            start_datetime = datetime.strptime(analysis_start, '%Y-%m-%d %H:%M')
            end_datetime = datetime.strptime(analysis_end, '%Y-%m-%d %H:%M')
            df = df.filter(
                    pl.col("timestamp").is_between(start_datetime, end_datetime)
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
            if axes_scale:
                subpt.set_xlim(x_origin, xdomain)
                subpt.set_ylim(y_origin, yrange)
            #subpt.set_xlim(0,200)
            #subpt.set_ylim(0,200)
            subpt.grid(color='green', linestyle='--', linewidth=0.25)
            exp_plot = subpt.scatter(x[plot], y[plot], c=z[plot], s=50, cmap=cm.jet, alpha=.3, marker='.')
            subpt.axline((0, 0), slope=1, color='black', linestyle=':', linewidth=.5, alpha=.7)
            plt.colorbar(exp_plot, label='probability density function')

    else:
    #patient data plots
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
            df = df.with_columns(
                (pl.col('Timestamp (YYYY-MM-DDThh:mm:ss)')).alias('timestamp')
            )
        else:
            # this is a libre file so....
            df = pl.read_csv(csvimport, has_header=True, try_parse_dates=True, skip_rows=1)
            datetest = df.select('Device Timestamp').row(0)[0]

            if re.match(patternl, datetest):
                df = df.with_columns(
                    (pl.col('Device Timestamp').str.strptime(pl.Datetime, fmt='%m-%d-%Y %I:%M %p')).alias('timestamp'))
                # print('usa version')
            else:
                # print('euro version')
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

        if daterange:
            dstart = df.head(1).select('timestamp')
            dend = df.tail(1).select('timestamp')
            print('Date Range of uploaded file is from ' + str(dstart[0, 0]) + ' to ' + str(dend[0, 0]))
            analysis_start = input('Enter analysis start date/time (following format YYYY-MM-DD HH:mm): ')
            analysis_end = input('Enter analysis end date/time (same format): ')
            start_datetime = datetime.strptime(analysis_start, '%Y-%m-%d %H:%M')
            end_datetime = datetime.strptime(analysis_end, '%Y-%m-%d %H:%M')
            df = df.filter(
                pl.col("timestamp").is_between(start_datetime, end_datetime)
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
            if axes_scale:
                subpt.set_xlim(x_origin, xdomain)
                subpt.set_ylim(y_origin, yrange)
            subpt.grid(color='green', linestyle='--', linewidth=0.25)
            exp_plot = subpt.scatter(x[plot], y[plot], c=z[plot], s=50, cmap=cm.jet, alpha=.3, marker='.')
            subpt.axline((0, 0), slope=1, color='black', linestyle=':', linewidth=.5, alpha=.7)
            plt.colorbar(exp_plot, label='probability density function')
    plt.show()
