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
from scipy import interpolate
from scipy import integrate


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
                    splitted_line = line.split(',')
                    x = float(splitted_line[0])
                    y = float(splitted_line[1])
                    xdata.append(x)
                    ydata.append(y)
                except:
                    splitted_line = line.split()
                    x_with_comma = splitted_line[0].split(',')
                    x_with_dot = x_with_comma[0] + '.' + x_with_comma[1]
                    x = float(x_with_dot)
                    y = float(splitted_line[1])
                    xdata.append(x)
                    ydata.append(y)
    return xdata, ydata

def lifetime_handler(inp_data):
    try:
        data, info = sif_reader.np_open(inp_data)
        data_fail_test = data[1][0] # Will fail if it is a single .sif spectrum
        type = 2
        return data, info, type
    except:
        data, info = sif_reader.np_open(inp_data)
        type = 1
        return data, info, type

def energy_corrector(data, i, xdata, energydata):
    out_energy_corr_data = []
    with open(energydata, 'r') as f:
        content = f.readlines()
        e_xdata = []
        e_ydata = []
        for line in content:
            splitted_line = line.split()
            x = float(splitted_line[0])
            y = float(splitted_line[1])
            e_xdata.append(x)
            e_ydata.append(y)
    if xdata != e_xdata:
        print('Warning, Energydata and x axis dont match!')
    else:
        for j,value in enumerate(data):
            out_energy_corr_data.append(value/e_ydata[j])
    return out_energy_corr_data

def bg_correction(data, xdata):
    out_bg_corr_data = []
    for i,value in enumerate(xdata):
        if value > 575:
            low_range_max = i
            break
    for j,value in enumerate(xdata):
        if value > 581:
            high_range_min = j
            break
    #y = ax + b
    low_min_x = xdata[data.index(min(data[0:low_range_max]))]
    high_min_x = xdata[data.index(min(data[high_range_min:]))]
    a = (min(data[high_range_min:])-min(data[0:low_range_max]))/(high_min_x-low_min_x)
    b = min(data[0:low_range_max])-a*low_min_x

    for i,value in enumerate(data):
        out_bg_corr_data.append(value-(a*xdata[i]+b))
    return out_bg_corr_data

def norm_data(data, decider, xdata=[], xlowlim=None, xuplim=None):
    out_norm_data = []
    if decider == 1:
        for value in data:
            out_value = value/max(data)
            out_norm_data.append(out_value)
        return out_norm_data

    elif decider == 2:
        for i,value in enumerate(xdata):
            if value > 600:
                low_value = i
                break
        for value in data:
            out_value = value/max(data[low_value:])
            out_norm_data.append(out_value)
        return out_norm_data

    elif decider ==3:
        if xlowlim != None:
            if xlowlim < xdata[0]:
                print(f'lower limit {xlowlim} out of range\n'
                f'min value is {xdata[0]}')
            elif xuplim > xdata[-1]:
                print(f'upper limit {xuplim} out of range\n'
                f'max value is {xdata[-1]}')
                xuplim = xdata[-1]
            for i,value in enumerate(xdata):
                if value > xlowlim:
                    low_value = i
                    break
            for j,value in enumerate(xdata):
                if value >= xuplim:
                    upper_value = j-1
                    break
            interpolate_function = interpolate.interp1d(xdata[low_value:upper_value], data[low_value:upper_value])
            integral = integrate.quad(interpolate_function, xdata[low_value], xdata[upper_value-1], limit=500)
            for value in data:
                out_norm_data.append(value/integral[0])
        else:
            interpolate_function = interpolate.interp1d(xdata, data)
            integral = integrate.quad(interpolate_function, xdata[0], xdata[-1], limit=500)

            for value in data:
                out_norm_data.append(value/integral[0])
        return out_norm_data

def deconvolute_data(data1, data2, fraction=1):
    out_deconv_data = []
    for i,value in enumerate(data1):
        out_value = float(value) - float(data2[i])*fraction
        out_deconv_data.append(out_value)
    return out_deconv_data

def Ne_corrector(xdata, Ne_data, ref_Ne_data):
    if Ne_data != ref_Ne_data:
        with open(ref_Ne_data, 'r') as f:
            content = f.readlines()
            Ne_ref_xdata = []
            Ne_ref_ydata = []
            for line in content:
                splitted_line = line.split()
                x = float(splitted_line[0])
                y = float(splitted_line[1])
                Ne_ref_xdata.append(x)
                Ne_ref_ydata.append(y)
        with open(Ne_data, 'r') as f:
            content = f.readlines()
            Ne_xdata = []
            Ne_ydata = []
            for line in content:
                splitted_line = line.split()
                x = float(splitted_line[0])
                y = float(splitted_line[1])
                Ne_xdata.append(x)
                Ne_ydata.append(y)
        for i,value in enumerate(Ne_ref_xdata):
            if value > 612:
                x_ref_lower_limit = i
                break
        for j, value in enumerate(Ne_ref_xdata):
            if value > 615:
                x_ref_upper_limit = j
                break

        ref_max = Ne_ref_ydata.index(max(Ne_ref_ydata[x_ref_lower_limit: x_ref_upper_limit]))

        for k,value in enumerate(Ne_xdata):
            if value > 612:
                x_lower_limit = k
                break
        for n, value in enumerate(Ne_xdata):
            if value > 615:
                x_upper_limit = n
                break
        Ne_max = Ne_ydata.index(max(Ne_ydata[x_lower_limit: x_upper_limit]))
        Ne_x_corr = []
        for value in xdata:
            Ne_x_corr.append(value + Ne_xdata[Ne_max] - Ne_ref_xdata[ref_max])

    else:
        Ne_x_corr = xdata

    return Ne_x_corr



def x_corrector(xdata):
    input_list = []
    for item in xdata:
        input_list.append(int(100*float(item)))
    output = []
    for value in input_list:
        if value in range(57200,57310):
            corr_value = value - 11
            output.append(float(corr_value/100))
        elif value in range(57310,57350):
            corr_value = value - 12
            output.append(float(corr_value/100))
        elif value in range(57350,57550):
            corr_value = value - 13
            output.append(float(corr_value/100))
        elif value in range(57550,57610):
            corr_value = value - 14
            output.append(float(corr_value/100))
        elif value in range(57610,57735):
            corr_value = value - 15
            output.append(float(corr_value/100))
        elif value in range(57735,57830):
            corr_value = value - 16
            output.append(float(corr_value/100))
        elif value in range(57830,57940):
            corr_value = value - 17
            output.append(float(corr_value/100))
        elif value in range(57940,58070):
            corr_value = value - 18
            output.append(float(corr_value/100))
        elif value in range(58070,58225):
            corr_value = value - 19
            output.append(float(corr_value/100))
        elif value in range(58225,58400):
            corr_value = value - 20
            output.append(float(corr_value/100))
        elif value in range(58400,59130):
            corr_value = value - 21
            output.append(float(corr_value/100))
        elif value in range(59130,59275):
            corr_value = value - 22
            output.append(float(corr_value/100))
        elif value in range(59275,60355):
            corr_value = value - 23
            output.append(float(corr_value/100))
        elif value in range(60355,60500):
            corr_value = value - 22
            output.append(float(corr_value/100))
        elif value in range(60500,60775):
            corr_value = value - 21
            output.append(float(corr_value/100))
        elif value in range(60775,61140):
            corr_value = value - 20
            output.append(float(corr_value/100))
        elif value in range(61140,61300):
            corr_value = value - 19
            output.append(float(corr_value/100))
        elif value in range(61300,61435):
            corr_value = value - 18
            output.append(float(corr_value/100))
        elif value in range(61435,61500):
            corr_value = value - 17
            output.append(float(corr_value/100))
        elif value in range(61500,61900):
            corr_value = value - 15
            output.append(float(corr_value/100))
        elif value in range(61900,62000):
            corr_value = value - 14
            output.append(float(corr_value/100))
        elif value in range(62000,62100):
            corr_value = value - 13
            output.append(float(corr_value/100))
        elif value in range(62100,62200):
            corr_value = value - 12
            output.append(float(corr_value/100))
        elif value in range(62200,62305):
            corr_value = value - 11
            output.append(float(corr_value/100))
        elif value in range(62305,62405):
            corr_value = value - 1
            output.append(float(corr_value/100))
        elif value in range(62405,62450):
            corr_value = value - 9
            output.append(float(corr_value/100))
        elif value in range(62450,62550):
            corr_value = value - 8
            output.append(float(corr_value/100))
        elif value in range(62550,62650):
            corr_value = value - 7
            output.append(float(corr_value/100))
        elif value in range(62650,62750):
            corr_value = value - 6
            output.append(float(corr_value/100))
        elif value in range(62750,62840):
            corr_value = value - 5
            output.append(float(corr_value/100))
        elif value in range(62840,62900):
            corr_value = value - 4
            output.append(float(corr_value/100))
        elif value in range(62900,62950):
            corr_value = value - 3
            output.append(float(corr_value/100))
        elif value in range(62950,63050):
            corr_value = value - 2
            output.append(float(corr_value/100))
        elif value in range(63050,63155):
            corr_value = value - 1
            output.append(float(corr_value/100))
        elif value in range(63155,63225):
            corr_value = value
            output.append(float(corr_value/100))
        elif value in range(62650,63255):
            corr_value = value + 1
            output.append(float(corr_value/100))
        elif value in range(63255,63405):
            corr_value = value + 2
            output.append(float(corr_value/100))
        elif value in range(63405,63450):
            corr_value = value + 3
            output.append(float(corr_value/100))
        elif value in range(63450,63575):
            corr_value = value + 4
            output.append(float(corr_value/100))
        elif value in range(63575,63625):
            corr_value = value + 5
            output.append(float(corr_value/100))
        elif value in range(63625,63725):
            corr_value = value + 6
            output.append(float(corr_value/100))
        elif value in range(63725,63775):
            corr_value = value + 7
            output.append(float(corr_value/100))
        elif value in range(63775,63875):
            corr_value = value + 8
            output.append(float(corr_value/100))
        elif value in range(63875,63925):
            corr_value = value + 9
            output.append(float(corr_value/100))
        elif value in range(63925,63975):
            corr_value = value + 1
            output.append(float(corr_value/100))
        elif value in range(63975,64030):
            corr_value = value + 11
            output.append(float(corr_value/100))
        elif value in range(64030,64090):
            corr_value = value + 12
            output.append(float(corr_value/100))
        elif value in range(64090,64175):
            corr_value = value + 13
            output.append(float(corr_value/100))
        elif value in range(64175,64225):
            corr_value = value + 14
            output.append(float(corr_value/100))
        elif value in range(64225,64275):
            corr_value = value + 15
            output.append(float(corr_value/100))
        elif value in range(64275,64325):
            corr_value = value + 16
            output.append(float(corr_value/100))
        elif value in range(64325,64375):
            corr_value = value + 17
            output.append(float(corr_value/100))
        elif value in range(64375,64425):
            corr_value = value + 18
            output.append(float(corr_value/100))
        elif value in range(64425,64525):
            corr_value = value + 19
            output.append(float(corr_value/100))
        elif value in range(64525,64625):
            corr_value = value + 20
            output.append(float(corr_value/100))
        elif value in range(64625,64675):
            corr_value = value + 21
            output.append(float(corr_value/100))
        elif value in range(64675,64720):
            corr_value = value + 22
            output.append(float(corr_value/100))
        elif value in range(64720,64740):
            corr_value = value + 23
            output.append(float(corr_value/100))
        elif value in range(64740,64775):
            corr_value = value + 24
            output.append(float(corr_value/100))
        elif value in range(64775,64825):
            corr_value = value + 25
            output.append(float(corr_value/100))
        elif value in range(64825,64875):
            corr_value = value + 26
            output.append(float(corr_value/100))
        elif value in range(64875,64925):
            corr_value = value + 27
            output.append(float(corr_value/100))
        elif value in range(64925,64975):
            corr_value = value + 28
            output.append(float(corr_value/100))
        elif value in range(64975,65090):
            corr_value = value + 29
            output.append(float(corr_value/100))
        elif value in range(65090,65125):
            corr_value = value + 3
            output.append(float(corr_value/100))
        elif value in range(65125,65175):
            corr_value = value + 31
            output.append(float(corr_value/100))
        elif value in range(65175,65225):
            corr_value = value + 32
            output.append(float(corr_value/100))
        elif value in range(65225,65275):
            corr_value = value + 33
            output.append(float(corr_value/100))
        elif value in range(65275,65360):
            corr_value = value + 34
            output.append(float(corr_value/100))
        elif value in range(65360,65410):
            corr_value = value + 35
            output.append(float(corr_value/100))
        elif value in range(65410,65450):
            corr_value = value + 36
            output.append(float(corr_value/100))
        elif value in range(65450,65500):
            corr_value = value + 37
            output.append(float(corr_value/100))
        elif value in range(65500,65550):
            corr_value = value + 38
            output.append(float(corr_value/100))
        elif value in range(65550,65610):
            corr_value = value + 39
            output.append(float(corr_value/100))
        elif value in range(65610,65640):
            corr_value = value + 4
            output.append(float(corr_value/100))
        elif value in range(65640,65710):
            corr_value = value + 41
            output.append(float(corr_value/100))
        elif value in range(65710,65750):
            corr_value = value + 42
            output.append(float(corr_value/100))
        elif value in range(65750,65810):
            corr_value = value + 43
            output.append(float(corr_value/100))
        elif value in range(65810,65835):
            corr_value = value + 44
            output.append(float(corr_value/100))
        elif value in range(65835,65845):
            corr_value = value + 45
            output.append(float(corr_value/100))
        elif value in range(65850,65900):
            corr_value = value + 46
            output.append(float(corr_value/100))
        elif value in range(65900,65925):
            corr_value = value + 47
            output.append(float(corr_value/100))
        elif value in range(65925,65975):
            corr_value = value + 48
            output.append(float(corr_value/100))
        elif value in range(65975,66025):
            corr_value = value + 49
            output.append(float(corr_value/100))
        elif value in range(66025,66075):
            corr_value = value + 50
            output.append(float(corr_value/100))
        elif value in range(66075,66125):
            corr_value = value + 51
            output.append(float(corr_value/100))
        elif value in range(66125,66375):
            corr_value = value + 52
            output.append(float(corr_value/100))
        elif value in range(66375,66475):
            corr_value = value + 51
            output.append(float(corr_value/100))
        elif value in range(66475,66525):
            corr_value = value + 50
            output.append(float(corr_value/100))
        elif value in range(66525,66675):
            corr_value = value + 49
            output.append(float(corr_value/100))
        elif value in range(66675,66725):
            corr_value = value + 48
            output.append(float(corr_value/100))
        elif value in range(66725,66775):
            corr_value = value + 47
            output.append(float(corr_value/100))
        elif value in range(66775,66805):
            corr_value = value + 46
            output.append(float(corr_value/100))
        else:
            output.append('Not in Range')
    return output