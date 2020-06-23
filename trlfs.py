"""
Module for all the TRLFS functions.

energy_correct: divides every value by the energy provided. For that a energydata file is needed in a tab separated
                design with x = excitation wavelength and y = energy

bg_corrector: draws a line between the lowest intensity in the range from 0 - 575nm and the lowest point between 575 - 581 nm
                and subtracts this line from all points.

norm_data: normalizes the data as
            1) divided by max
            2) F2 corrected (i.e. the max value in the range of 600 to the end of the data set
            3) divided by the integrated are of the spectrum after interpolation between the data points


"""
import sif_reader
from sif_reader import extract_calibration
import os


def load_data(inp_data):
    if inp_data.split('.')[-1] == 'sif':
        data, info = sif_reader.np_open(inp_data)
        xdata = extract_calibration(info)
        ydata = data[0][0]
    else:
        xdata = []
        ydata = []
        with open(inp_data, 'r') as f:
            data = f.readlines()
        for line in data:
            try:
                splitted_line = line.split()
                x = float(splitted_line[0])
                y = float(splitted_line[1])
                xdata.append(x)
                ydata.append(y)
            except:
                try:
                    splitted_line = line.split()
                    x_with_comma = splitted_line[0].split(',')
                    x_with_dot = x_with_comma[0] + '.' + x_with_comma[1]
                    x = float(x_with_dot)
                    y = float(splitted_line[1])
                    xdata.append(x)
                    ydata.append(y)
                except:
                    splitted_line = line.split(',')
                    x = float(splitted_line[0])
                    y = float(splitted_line[1])
                    xdata.append(x)
                    ydata.append(y)
    return xdata, ydata

def load_multi_column_data(inp_data):
    xdata = []
    ydata = []
    number_of_y_columns = 0
    with open(inp_data, 'r') as f:
        data = f.readlines()
    for line in data:
        try:
            splitted_line = line.split()
            number_of_y_columns = len(splitted_line) - 1
            x = float(splitted_line[0])
            xdata.append(x)
            for i in range(number_of_y_columns):
                y = float(splitted_line[i + 1])
                ydata.append([])
                ydata[i].append(y)
        except:
            try:
                splitted_line = line.split()
                number_of_y_columns = len(splitted_line) - 1
                x_with_comma = splitted_line[0].split(',')
                x_with_dot = x_with_comma[0] + '.' + x_with_comma[1]
                x = float(x_with_dot)
                xdata.append(x)
                for i in range(number_of_y_columns):
                    y = float(splitted_line[i + 1])
                    ydata.append([])
                    ydata[i].append(y)
            except:
                splitted_line = line.split(',')
                number_of_y_columns = len(splitted_line) - 1
                x = float(splitted_line[0])
                xdata.append(x)
                for i in range(number_of_y_columns):
                    y = float(splitted_line[i + 1])
                    ydata.append([])
                    ydata[i].append(y)
    return xdata, ydata, number_of_y_columns

def lifetime_handler(inp_data):
    if inp_data.split('.')[-1] == 'sif':
        try:
            data, info = sif_reader.np_open(inp_data)
            data_fail_test = data[1][0] # Will fail if it is a single .sif spectrum
            type = 2
            return data, info, type
        except:
            data, info = sif_reader.np_open(inp_data)
            type = 1
            return data, info, type
    else:
        try:
            xdata, ydata, number_of_y_columns = load_multi_column_data(inp_data)
            data = ydata
            info = xdata
            type = 3
        except:
            data = []
            info = []
            type = 1
        return data, info, type



def load_wavelength_data(inp_data):
    wavelength_data = []
    with open(inp_data, 'r+') as w:
        for line in w:
            if len(line.split()) == 1:  # Check whether the data is single column x-values
                wavelength_data.append(float(line.split('\n')[0]))
            else:                       # If not just take first column
                wavelength_data.append(float(line.split()[0]))
    return wavelength_data

