"""
This program reads spectroscopic x,y data (.asc, .txt, or AndorFile .sif) of single or multiple spectra.
This data can be background_corrected, integrated, saved and plotted.

Started: 10.09.2019
Last edit: 10.06.2020
Version: 1.1.0
Creator: Manuel Eibl
"""

import os
import shutil
import ctypes
import tkinter as tk
from tkinter import filedialog
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
matplotlib.use('TkAgg')
import numpy as np
from scipy import interpolate, integrate
from sif_reader import extract_calibration
import trlfs


# Handles all the data input, access, integration and output
class DataFiles():
    def __init__(self, dir='/', temp_dirs=[] , file_list=[], file_dict={}, integral_xdata=[], integral_ydata=[]):
        self.data_dir = dir
        self.temp_dirs = temp_dirs
        self.file_list = file_list  # list of all files with exact paths
        self.file_dict = file_dict  # dict, key_word = file name, value = path
        self.integral_xdata = integral_xdata
        self.integral_ydata = integral_ydata
        self.info = []

        #load last opened directory for quicker file search
        saved_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'log.dat')
        if os.path.exists(saved_path):
            with open(saved_path, 'r') as r:
                content = r.readlines()
            try:
                self.data_dir = content[0]
            except:
                pass

    def add_files_from_list(self, new_files=[]):
        for file in new_files:
            self.file_list.append(file)
            self.data_dir, file_name = os.path.split(file)
            self.file_dict[file_name] = file
        return self.file_list

    def get_list_of_all_files(self):
        return self.file_list

    def get_list_of_all_file_names(self):
        return self.file_dict.keys()

    def get_file_from_file_name(self, key_word):
        return self.file_dict[key_word]

    def get_list_of_files_from_indexes(self, index):
        selected_files = []
        try:
            for value in index:
                selected_files.append(self.file_list[value])
        except:
            selected_files.append(self.file_list[index])
        return selected_files

    def get_data_dir(self):
        return self.data_dir

    # Write last used directory for initialization above
    def write_data_dir(self):
        script_path = os.path.dirname(os.path.abspath('__file__'))
        log_file_name = 'log.dat'
        log_file = os.path.join(script_path, log_file_name)
        if os.path.split(self.data_dir)[1] == '.temp':
            self.data_dir = os.path.split(self.data_dir)[0]
        with open (log_file, 'w+') as f:
            f.write(self.data_dir)

    def short_cut_remove_file(self, event):
        selected_indexes = lister.curselection()
        for index in selected_indexes[::-1]:    # starting from back to avoid change of indexes during removing
            file_name = lister.get(index)
            file = self.get_file_from_file_name(file_name)
            del self.file_dict[file_name]
            for i, value in enumerate(self.file_list):
                if file == value:
                    self.file_list.pop(i)
        plot_all(self.file_list)

    def clear_all(self):
        self.file_list = []
        self.file_dict = {}
        self.integral_xdata = []
        self.integral_ydata = []
        self.info = []
        lister.delete(0, tk.END)

    def integrate_data(self):
        button_excitation_wavelengths['state'] = 'normal'
        self.integral_xdata = []
        self.integral_ydata = []
        line1_x_value = None
        line2_x_value = None
        if integration_lines.get_range() is not None:
            line1_x_value, line2_x_value = integration_lines.get_range()
        selected_item_indexes = lister.curselection()
        for i, item_index in enumerate(selected_item_indexes):
            file_name = lister.get(item_index)
            file = self.get_file_from_file_name(file_name)
            stripped_file_name = os.path.split(file)[1].split('.')[0]
            try:
                if stripped_file_name.split('_')[-1] == 'corr':
                    first_value = int(stripped_file_name.split('_')[-3])
                    if first_value > 30000:         # which means it is excitation data
                        xvalue = float(stripped_file_name.split('_')[-3]) / 100
                    else:                           # which means it is lifetime data
                        xvalue = int(stripped_file_name.split('_')[-3])
                else:
                    first_value = int(stripped_file_name.split('_')[-1])
                    if first_value > 30000:         # which means it is excitation data
                        xvalue = float(stripped_file_name.split('_')[-1]) / 100
                    else:                           # which means it is lifetime data
                        xvalue = int(stripped_file_name.split('_')[-1])
            except:
                xvalue = i + 1
            self.integral_xdata.append(xvalue)
            if background.background_corr_data == []:
                xdata, ydata = trlfs.load_data(file)
            else:
                xdata, ydata = background.background_corr_data[i]
            if line1_x_value != None:
                for j, value in enumerate(xdata):   # find xdata value closest to line1 x value
                    if value >= line1_x_value:
                        lower_value = j
                        break
                for k, value in enumerate(xdata):   # find xdata value closest to line2 x value
                    if value >= line2_x_value:
                        upper_value = k
                        break
            else:     # choose whole range
                lower_value = 0
                upper_value = -1
            interpolate_function = interpolate.interp1d(xdata[lower_value:upper_value],
                                                        ydata[lower_value:upper_value])
            integral_of_file = integrate.quad(interpolate_function, xdata[lower_value], xdata[upper_value - 1], limit=50)[0]
            self.integral_ydata.append(integral_of_file)
        integral_plotter(self.integral_xdata, self.integral_ydata)

    def integrate_all_data(self):
        button_excitation_wavelengths['state'] = 'normal'
        self.integral_xdata = []
        self.integral_ydata = []
        line1_x_value = None
        line2_x_value = None
        if integration_lines.get_range() is not None:
            line1_x_value, line2_x_value = integration_lines.get_range()
        all_files = lister.get(0, tk.END)
        for i, file_name in enumerate(all_files):
            file = self.get_file_from_file_name(file_name)
            stripped_file_name = os.path.split(file)[1].split('.')[0]
            try:
                if stripped_file_name.split('_')[-1] == 'corr':
                    first_value = int(stripped_file_name.split('_')[-3])
                    if first_value > 30000:         # which means it is excitation data
                        xvalue = float(stripped_file_name.split('_')[-3]) / 100
                    else:                           # which means it is lifetime data
                        xvalue = int(stripped_file_name.split('_')[-3])
                else:
                    first_value = int(stripped_file_name.split('_')[-1])
                    if first_value > 30000:         # which means it is excitation data
                        xvalue = float(stripped_file_name.split('_')[-1]) / 100
                    else:                           # which means it is lifetime data
                        xvalue = int(stripped_file_name.split('_')[-1])
            except:
                xvalue = i + 1
            self.integral_xdata.append(xvalue)
            if background.background_corr_data == []:
                xdata, ydata = trlfs.load_data(file)
            else:
                xdata, ydata = background.background_corr_data[i]
            if line1_x_value != None:
                for j, value in enumerate(xdata):   # find xdata value closest to line1 x value
                    if value >= line1_x_value:
                        lower_value = j
                        break
                for k, value in enumerate(xdata):   # find xdata value closest to line2 x value
                    if value >= line2_x_value:
                        upper_value = k
                        break
            else:     # choose whole range
                lower_value = 0
                upper_value = -1
            interpolate_function = interpolate.interp1d(xdata[lower_value:upper_value],
                                                        ydata[lower_value:upper_value])
            integral_of_file = integrate.quad(interpolate_function, xdata[lower_value], xdata[upper_value - 1], limit=50)[0]
            self.integral_ydata.append(integral_of_file)
        integral_plotter(self.integral_xdata, self.integral_ydata)

    def integral_output_maker(self):
        output = ''
        for i, x_value in enumerate(self.integral_xdata):
            output += str(x_value) + '\t' + str(self.integral_ydata[i]) + '\n'
        return output

    def save_integration_data(self):
        first_file_name = lister.get(0)
        if integration_lines.get_range() is not None:
            line1_x_value, line2_x_value = integration_lines.get_range()
        else:
            file = data_files.get_file_from_file_name(first_file_name)
            xdata, _ = trlfs.load_data(file)
            line1_x_value = xdata[0]
            line2_x_value = xdata[-1]
        out_file_name = os.path.split(first_file_name)[1].split('.')[0]
        if out_file_name.split('_')[0] != out_file_name:
            temp_out_file_name = ''
            for i in range(len(out_file_name.split('_')) - 1):
                temp_out_file_name += out_file_name.split('_')[i]
            out_file_name = temp_out_file_name
        output = self.integral_output_maker()
        init_file_name = f'{out_file_name} integration from {int(line1_x_value*100)/100} - {int(line2_x_value*100)/100}'
        save_file = filedialog.asksaveasfilename(initialdir=self.get_data_dir(), initialfile=init_file_name, defaultextension='.txt', filetypes=[('Text file', '*.txt'), ('Others', ('*.asc', '*dat'))])
        with open(str(save_file), 'w+') as o:
            o.write(str(output))

    def lifetime_converter(self, lifetime_file):
        output_files_list = []
        file_path, file_name = os.path.split(lifetime_file)
        temp_path = os.path.join(file_path, '.temp')
        if not os.path.exists(temp_path):
            os.mkdir(temp_path)
            ctypes.windll.kernel32.SetFileAttributesW(temp_path, 2) # hidden dir
            self.temp_dirs.append(temp_path)    # for deletion of temp dir on closing of GUI
        data, info, type = trlfs.lifetime_handler(lifetime_file)
        if type == 3:
            xdata, ydata, number_of_frames = trlfs.load_multi_column_data(lifetime_file)
            for i in range(number_of_frames):
                temp_file_name = file_name.split('.')[0] + '_' + str(int(i + 1)) + '.asc'
                temp_file = os.path.join(temp_path, temp_file_name)
                output = ''
                ydata = data[i]
                for j, x_value in enumerate(xdata):
                    output += str(x_value) + '\t' + str(ydata[j]) + '\n'
                with open(temp_file, 'w+') as f:
                    f.write(output)
                output_files_list.append(temp_file)
        else:
            xdata = extract_calibration(info)
            for i in range(info['NumberOfFrames']):
                temp_file_name = file_name.split('.')[0] + '_' + str(int(1 + (i*info['GateDelayStep']*10**-6))) + '.asc'
                temp_file = os.path.join(temp_path, temp_file_name)
                output = ''
                ydata = data[i][0]
                for j, x_value in enumerate(xdata):
                    output += str(x_value) + '\t' + str(ydata[j]) + '\n'
                with open(temp_file, 'w+') as f:
                    f.write(output)
                output_files_list.append(temp_file)
        return output_files_list

    def save_selected_files(self):
        output_files_list = []
        selected_item_indexes = lister.curselection()
        for item_index in selected_item_indexes:
            file_name = lister.get(item_index)
            temp_file_name = file_name.split('.')[0] + '.asc'
            file = data_files.get_file_from_file_name(file_name)
            file_path, _ = os.path.split(file)
            if os.path.split(file_path)[1] == '.temp':
                temp_path = file_path
                temp_file = os.path.join(temp_path, temp_file_name)
            else:
                temp_path = os.path.join(file_path, '.temp')
                temp_file = os.path.join(temp_path, temp_file_name)
                output_files_list.append(temp_file)
                if not os.path.exists(temp_path):
                    os.mkdir(temp_path)
                    ctypes.windll.kernel32.SetFileAttributesW(temp_path, 2) # hidden dir
                    self.temp_dirs.append(temp_path)    # for deletion of temp dir on closing of GUI
            output = ''
            xdata, ydata = trlfs.load_data(file)
            for i, x_value in enumerate(xdata):
                output += str(x_value) + '\t' + str(ydata[i]) + '\n'
            with open(temp_file, 'w+') as o:
                o.write(output)

        abs_path = filedialog.askdirectory(initialdir=self.get_data_dir())
        for temp_file in output_files_list:
            _, temp_filename = os.path.split(temp_file)
            abs_out_path = os.path.join(abs_path, temp_filename)
            shutil.copyfile(temp_file, abs_out_path)

    def save_all_files(self):
        output_files_list = []
        list_of_all_files = self.get_list_of_all_files()
        for file in list_of_all_files:
            file_path, file_name = os.path.split(file)
            if os.path.split(file_path)[1] == '.temp':
                temp_path = file_path
                temp_file = os.path.join(temp_path, file_name)
            else:
                temp_path = os.path.join(file_path, '.temp')
                temp_file = os.path.join(temp_path, file_name)
                output_files_list.append(temp_file)
                if not os.path.exists(temp_path):
                    os.mkdir(temp_path)
                    ctypes.windll.kernel32.SetFileAttributesW(temp_path, 2) # hidden dir
                    self.temp_dirs.append(temp_path)    # for deletion of temp dir on closing of GUI
            output = ''
            xdata, ydata = trlfs.load_data(file)
            for i, x_value in enumerate(xdata):
                output += str(x_value) + '\t' + str(ydata[i]) + '\n'
            with open(temp_file, 'w+') as o:
                o.write(output)

        abs_path = filedialog.askdirectory(initialdir=self.get_data_dir())
        for temp_file in output_files_list:
            _, temp_filename = os.path.split(temp_file)
            abs_out_path = os.path.join(abs_path, temp_filename)
            shutil.copyfile(temp_file, abs_out_path)

    def get_acquisition_info(self):
        selected_file_index = lister.curselection()
        if selected_file_index == () and data_files.info != []:
            selected_file = data_files.get_list_of_files_from_indexes([0])[0]
        else:
            if len(selected_file_index) > 1:
                print('More than 1 file selected')
                return
            selected_file = self.get_list_of_files_from_indexes(selected_file_index)[0]
            if selected_file.split('.')[-1] != 'sif' and data_files.info == []:
                print('Not a .sif file. Other file types currently not supported for Acquisition info')
                return
        self.display_acquisition_info(selected_file)

    def display_acquisition_info(self, selected_file):
        if self.info == []:
            data, info, type = trlfs.lifetime_handler(selected_file)
            _, file_name = os.path.split(selected_file)
            if type == 2:
                file_type = 'Lifetime measurement'
            else:
                file_type = 'Single spectrum'
        else:
            _, temp_file_name = os.path.split(selected_file)
            file_name = temp_file_name.rsplit('_', 1)[0] + '.sif'
            info = self.info
            file_type = 'Lifetime measurement'

        self.display_acquisition_window = tk.Toplevel(root)
        self.display_acquisition_window.title('Acquisition information')

        # File name
        label_file_name = tk.Label(self.display_acquisition_window, text='File name')
        label_file_name.grid(row=0, column=0)
        value_file_name = tk.Label(self.display_acquisition_window, text=file_name)
        value_file_name.grid(row=0, column=1)
        # Lifetime or Single spectrum
        label_file_type = tk.Label(self.display_acquisition_window, text='File type')
        label_file_type.grid(row=1, column=0)
        value_file_type = tk.Label(self.display_acquisition_window, text=file_type)
        value_file_type.grid(row=1, column=1)
        #Accumulations
        label_accumulations = tk.Label(self.display_acquisition_window, text='Accumulations')
        label_accumulations.grid(row=2, column=0)
        value_accumulations = tk.Label(self.display_acquisition_window, text=info['AccumulatedCycles'])
        value_accumulations.grid(row=2, column=1)
        # Grating
        label_grating = tk.Label(self.display_acquisition_window, text='Grating')
        label_grating.grid(row=3, column=0)
        value_grating = tk.Label(self.display_acquisition_window, text="{:.0f} lines/mm".format(info['Grating'])) # This might not be the grating! Check that
        value_grating.grid(row=3, column=1)
        # Slit
        label_slit = tk.Label(self.display_acquisition_window, text='Slit opening')
        label_slit.grid(row=4, column=0)
        value_slit = tk.Label(self.display_acquisition_window, text=str(info['SlitOpening']) + ' um')
        value_slit.grid(row=4, column=1)
        # Gain
        label_gain = tk.Label(self.display_acquisition_window, text='Gain')
        label_gain.grid(row=5, column=0)
        value_gain = tk.Label(self.display_acquisition_window, text=str(int(info['GainDAC'])))
        value_gain.grid(row=5, column=1)
        # Gate delay
        label_delay = tk.Label(self.display_acquisition_window, text='Gate delay')
        label_delay.grid(row=6, column=0)
        value_delay = tk.Label(self.display_acquisition_window, text=(str(float(info['GateDelay'])/1e6) + ' us'))
        value_delay.grid(row=6, column=1)
        # Gate width
        label_width = tk.Label(self.display_acquisition_window, text='Gate width')
        label_width.grid(row=7, column=0)
        value_width = tk.Label(self.display_acquisition_window, text=(str(float(info['GateWidth'])/1e9) + ' ms'))
        value_width.grid(row=7, column=1)

        i = 8
        # Number of frames and steps if appliccable
        if file_type == 'Lifetime measurement':
            # Number of Frames
            label_number_of_frames = tk.Label(self.display_acquisition_window, text='Number of frames')
            label_number_of_frames.grid(row=i, column=0)
            value_number_of_frames = tk.Label(self.display_acquisition_window, text=info['NumberOfFrames'])
            value_number_of_frames.grid(row=i, column=1)
            # Delay step size
            label_step_size = tk.Label(self.display_acquisition_window, text='Step size')
            label_step_size.grid(row=i+1, column=0)
            value_step_size = tk.Label(self.display_acquisition_window, text=(str(float(info['GateDelayStep'])/1e6) + ' us'))
            value_step_size.grid(row=i+1, column=1)

            i += 2
        # Exposure time
        label_exposure_time = tk.Label(self.display_acquisition_window, text='Exposure time')
        label_exposure_time.grid(row=i, column=0)
        value_exposure_time = tk.Label(self.display_acquisition_window, text=(str(info['ExposureTime']) + ' s'))
        value_exposure_time.grid(row=i, column=1)
        # Center wavelength
        label_center_wavelength = tk.Label(self.display_acquisition_window, text='Center wavelength')
        label_center_wavelength.grid(row=i+1, column=0)
        value_center_wavelength = tk.Label(self.display_acquisition_window, text=(str(info['CenterWavelength']) + ' nm'))
        value_center_wavelength.grid(row=i+1, column=1)
        # Detector type
        label_detector = tk.Label(self.display_acquisition_window, text='Detector')
        label_detector.grid(row=i+2, column=0)
        value_detector = tk.Label(self.display_acquisition_window, text=info['DetectorType'].split()[0])
        value_detector.grid(row=i+2, column=1)


# This class creates draggable lines to define the integration area
class DraggableLines():
    lock = None

    def __init__(self, ax, canvas, lines=[]):
        self.ax = ax
        self.canvas = canvas
        self.lines = lines

    def init_line(self):
        self.lines = []
        # Make default lines first
        file = data_files.get_list_of_all_files()[0]
        xdata, _ = trlfs.load_data(file)
        center_x_value_index = int(len(xdata) / 2)
        x_value_line1 = xdata[center_x_value_index] - 5
        x_value_line2 = xdata[center_x_value_index] + 5
        left_line = self.ax.axvline(x_value_line1)
        right_line = self.ax.axvline(x_value_line2)
        self.lines.append(left_line)
        self.lines.append(right_line)
        self.canvas.draw()

    def click(self, event): # For selection of line, could picker be used instead?
        width, _ = self.canvas.get_width_height()
        xlim_lower, xlim_upper = self.ax.get_xlim()
        pix_per_nm = width / (xlim_upper - xlim_lower)
        for line in self.lines:
            try:
                len(line.get_xdata())
                line_xdata = line.get_xdata()[0]
            except:
                line_xdata = line.get_xdata()
            if np.abs(((line_xdata - event.xdata)*pix_per_nm)) < 10.0:
                DraggableLines.lock = line

    def release(self, event):
        if DraggableLines.lock is None:
            return
        else:
            DraggableLines.lock = None

    def drag(self, event):
        if DraggableLines.lock is None:
            return
        chosen_line = DraggableLines.lock
        chosen_line.set_xdata(event.xdata)
        self.canvas.draw()

    def get_range(self):
        if self.lines == []:
            return None
        out_xdata = []
        for line in self.lines:
            try:
                len(line.get_xdata())
                out_xdata.append(line.get_xdata()[0])
            except:
                out_xdata.append(line.get_xdata())
        return out_xdata[0], out_xdata[1]


# This class creates draggable points as used in the BG-correction function
class DraggablePoints():
    lock = None
    def __init__(self, ax, canvas, item_indexes):
        self.ax = ax
        self.canvas = canvas
        self.item_indexes = item_indexes
        self.x_coord = 0
        self.dots = []
        self.lines = []
        self.background_lines = None
        self.init_data_range()

    def init_data_range(self):
        first_data_file = data_files.get_list_of_files_from_indexes(self.item_indexes[0])[0]
        self.xdata, self.ydata = trlfs.load_data(first_data_file)

    def draw_background_line(self):
        if len(self.dots) > 1:  # can only draw a line once 2 points exist
            self.dots_xdata = []
            self.dots_ydata = []
            for dot in self.dots:
                try:
                    self.dots_xdata.append(dot.get_xdata()[0])
                    self.dots_ydata.append(dot.get_ydata()[0])
                except:
                    self.dots_xdata.append(dot.get_xdata())
                    self.dots_ydata.append(dot.get_ydata())
            # sorting from lowest to highest x value
            sorted_dots_data = sorted(zip(self.dots_xdata, self.dots_ydata))
            sorted_xdata = []
            sorted_ydata = []
            for data in sorted_dots_data:
                sorted_xdata.append(data[0])
                sorted_ydata.append(data[1])
            if self.background_lines is None:
                self.background_lines = self.ax.plot(sorted_xdata, sorted_ydata)
            else:
                self.background_lines[0].set_data(sorted_xdata, sorted_ydata)

    def click(self, event):
        if event.dblclick:
            self.x_coord = event.xdata
            closest_point_data = self.get_closest_data_point(self.x_coord)
            dot = self.ax.plot([closest_point_data[0]], [closest_point_data[1]], marker='o', color='red', markersize=5, linestyle='none', zorder=5, picker=5)[0]
            self.dots.append(dot)
            self.draw_background_line()
            self.canvas.draw()
        else:
            width, _ = canvas.get_width_height()
            xlim_lower, xlim_upper = ax.get_xlim()
            pix_per_nm = width / (xlim_upper - xlim_lower)
            for dot in self.dots:
                try:
                    len(dot.get_xdata())
                    dot_xdata = dot.get_xdata()[0]
                except:
                    dot_xdata = dot.get_xdata()
                if np.abs(((dot_xdata - event.xdata) * pix_per_nm)) < 10.0:
                    DraggablePoints.lock = dot

    def release(self, event):
        if DraggablePoints.lock is None:
            return
        else:
            DraggablePoints.lock = None
            self.draw_background_line()

    def select(self, event):
        checker = str(event.mouseevent.button)
        if checker == '3':
            select_x_value = event.artist.get_xdata()[0]
            for i, dot in enumerate(self.dots):
                if dot.get_xdata()[0] == select_x_value:
                    dot.remove()
                    self.dots.pop(i)
            self.draw_background_line()
            self.canvas.draw()

    def drag(self, event):
        if DraggablePoints.lock is None:
            return
        chosen_dot = DraggablePoints.lock
        closest_point_data = self.get_closest_data_point(event.xdata)
        chosen_dot.set_xdata(closest_point_data[0])
        chosen_dot.set_ydata(closest_point_data[1])
        self.draw_background_line()
        self.canvas.draw()

    def get_closest_data_point(self, x_value_point):
        closest_data_point = {}
        for i, x_value_data in enumerate(self.xdata):
            diff = np.abs(x_value_data - x_value_point)
            if i == 0:
                closest_data_point[0] = (i, diff)
            else:
                if diff < closest_data_point[0][1]:
                    closest_data_point[0] = (i, diff)
        closest_index = closest_data_point[0][0]
        closest_point_data = [self.xdata[closest_index], self.ydata[closest_index]]
        return closest_point_data

    def get_bg_data(self):
        return self.dots


# Optional feature of correcting linear background with n-1 lines between n points
class BackgroundCorrection():
    def __init__(self):
        self.background_corr_data = []

    # BG-correction window initialization
    def base_line_correction_init(self):
        self.fi2, self.ax2 = plt.subplots()
        self.correction_window = tk.Toplevel(root)
        self.correction_window.title('Background correction')
        self.correction_canvas = FigureCanvasTkAgg(self.fi2, master=self.correction_window)
        self.background_corr_data = []

        for x in range(10):
            self.correction_window.columnconfigure(x, weight=1)
        for y in range(22):
            self.correction_window.rowconfigure(y, weight=1)
        self.correction_canvas.get_tk_widget().grid(row=2, column=0, rowspan=20, columnspan=10, sticky='nsew')
        ####################################################################
        # Make Toolbar above figure
        ###################################################################
        toolbarFrame = tk.Frame(self.correction_window)
        toolbarFrame.grid(row=1, column=0, columnspan=7, sticky='wnse')
        toolbar = NavigationToolbar2Tk(self.correction_canvas, toolbarFrame)
        button_bg_corr = tk.Button(self.correction_window, text='Apply correction to all', command=self.bg_correct_all)
        button_bg_corr.grid(row=0, column=0, sticky='nsew')
        self.data_plotter()
        self.draw_points()
        self.correction_canvas.draw()
        self.correction_window.mainloop()

    def data_plotter(self):
        self.selected_item_indexes = lister.curselection()
        for item_index in self.selected_item_indexes:
            file_name = lister.get(item_index)
            file = data_files.get_file_from_file_name(file_name)
            xdata, ydata = trlfs.load_data(file)
            self.ax2.plot(xdata, ydata, label=file_name)
            if len(self.selected_item_indexes) <= 10:
                ax.legend()

    def draw_points(self):
        self.bg_corr_line = DraggablePoints(self.ax2, self.correction_canvas, self.selected_item_indexes)
        self.cid_click = self.fi2.canvas.mpl_connect('button_press_event', self.bg_corr_line.click)
        self.cid_select = self.fi2.canvas.mpl_connect('pick_event', self.bg_corr_line.select)
        self.cid_drag = self.fi2.canvas.mpl_connect('motion_notify_event', self.bg_corr_line.drag)
        self.cid_release = self.fi2.canvas.mpl_connect('button_release_event', self.bg_corr_line.release)
        self.correction_canvas.draw()

    def bg_correct_all(self):
        self.disconnect()
        dots_xdata = self.bg_corr_line.dots_xdata
        dots_ydata = self.bg_corr_line.dots_ydata
        bg_functions = {}   # dict of functions, valid below the according x_value, i.e. the key
        for i, dot_xdata in enumerate(dots_xdata):
            if i > 0:
                slope = (dots_ydata[i] - dots_ydata[i-1])/(dots_xdata[i] - dots_xdata[i-1])
                y_intercept = dots_ydata[i - 1] - slope * dots_xdata[i - 1]
                bg_functions[dot_xdata] = [slope, y_intercept]
        function_delimiters = list(bg_functions.keys())  # list of the x_values of the dots, first bg function valid until first value etc.
        self.list_of_all_files = data_files.get_list_of_all_files()
        for file in self.list_of_all_files:
            _, file_name = os.path.split(file)
            xdata, ydata = trlfs.load_data(file)
            corrected_ydata = []
            range_of_points = list(range(0, len(function_delimiters)))
            for i, x_value in enumerate(xdata):
                for j in range_of_points:   # step through function_delimiters and find corresponding bg function
                    if x_value <= function_delimiters[j]:
                        slope, y_intercept = bg_functions[function_delimiters[j]]
                        corrected_ydata.append(ydata[i] - (slope * xdata[i] + y_intercept))
                        break
                    elif x_value > function_delimiters[-1]:
                        slope, y_intercept = bg_functions[function_delimiters[-1]]
                        corrected_ydata.append(ydata[i] - (slope * xdata[i] + y_intercept))
                        break
                    else:
                        j += 1
            self.background_corr_data.append([xdata, corrected_ydata])
        ax.cla()
        for i, xy_data in enumerate(self.background_corr_data):
            file_name = list(data_files.get_list_of_all_file_names())[i]
            ax.plot(xy_data[0], xy_data[1], label=file_name.split('.')[0])
        self.correction_window.destroy()
        self.save_bg_corr_data()
        if len(self.background_corr_data) <= 10:
            ax.legend()
        canvas.draw()

    # Needs to be stored in /.temp/ so it can be replotted at any time
    def save_bg_corr_data(self):
        temp_files_list = []
        temp_files_dict = {}
        data_path = data_files.get_data_dir()
        if os.path.split(data_path)[1] == '.temp':
            temp_path = data_path
        else:
            temp_path = os.path.join(data_path, '.temp')
            if not os.path.exists(temp_path):
                os.mkdir(temp_path)
                ctypes.windll.kernel32.SetFileAttributesW(temp_path, 2)
                data_files.temp_dirs.append(temp_path)
        for i, xy_data in enumerate(self.background_corr_data):
            _, file_name = os.path.split(self.list_of_all_files[i])
            if file_name.split('_')[-1].split('.')[0] == 'corr':
                temp_file_name = file_name
            else:
                temp_file_name = file_name.split('.')[0] + '_bg_corr.asc'
            temp_file = os.path.join(temp_path, temp_file_name)
            output = ''
            for j, x_value in enumerate(xy_data[0]):
                output += str(x_value) + '\t' + str(xy_data[1][j]) + '\n'
            with open(temp_file, 'w+') as o:
                o.write(output)
            temp_files_list.append(temp_file)
            temp_files_dict[file_name] = temp_file
        data_files.file_list = temp_files_list
        data_files.file_dict = temp_files_dict

    def disconnect(self):
        self.fi2.canvas.mpl_disconnect(self.cid_click)
        self.fi2.canvas.mpl_disconnect(self.cid_drag)
        self.fi2.canvas.mpl_disconnect(self.cid_release)
        self.fi2.canvas.mpl_disconnect(self.cid_select)


# Adds all data to the lister box and plots it
def plot_all(file_list):
    ax.cla()
    integration_lines.lines = []
    lister.delete(0, tk.END)
    for file in file_list:
        _, file_name = os.path.split(file)
        lister.insert(tk.END, file_name)
        xdata, ydata = trlfs.load_data(file)
        ax.plot(xdata, ydata, linewidth=1, label=file_name.split('.')[0])
    if len(file_list) <= 10:
        ax.legend()
    canvas.draw()


def re_plotter(ax):
    ax.cla()
    integration_lines.lines = []
    integration_lines.x1 = None
    integration_lines.x2 = None
    selected_item_indexes = lister.curselection()
    for item_index in selected_item_indexes:
        file_name = lister.get(item_index)
        file_name = lister.get(item_index)
        file = data_files.get_file_from_file_name(file_name)
        xdata, ydata = trlfs.load_data(file)
        ax.plot(xdata, ydata, label=file_name.split('.')[0])
    if len(selected_item_indexes) <= 10:
        ax.legend()
    canvas.draw()


# For double click re-plotting
def short_cut_re_plotter(event, ax):
    ax.cla()
    integration_lines.lines = []
    integration_lines.x1 = None
    integration_lines.x2 = None
    selected_item_indexes = lister.curselection()
    for selected_item_index in selected_item_indexes:
        file_name = lister.get(selected_item_index)
        file = data_files.get_file_from_file_name(file_name)
        xdata, ydata = trlfs.load_data(file)
        ax.plot(xdata, ydata, label=file_name)
    if len(selected_item_indexes) <= 10:
        ax.legend()
    canvas.draw()


def integral_plotter(integral_xdata, integral_ydata):
    ax.cla()
    ax.plot(integral_xdata, integral_ydata, marker='o')
    canvas.draw()


# Opens a selection window and loads files. It checks whether the chosen file is .sif lifetime data, which has to
# be separated into individual files first (done by lifetime_converter()) then loads and plots the data
def browse_files():
    root.filename = filedialog.askopenfilenames(initialdir=data_files.get_data_dir(), title='Select files', filetype=(('All files', '*.*'), ('Text files', ('*.asc', '*.txt')), ('Andor .sif files', '*.sif')))
    new_files = list(root.filename)
    files_to_add = []
    data_files.info = []    # If lifetime measurements were loaded earlier, its acquisition information is cleared
    if len(new_files) == 1:
        file = new_files[0]
        data, info, type = trlfs.lifetime_handler(file)
        if type == 2:
            new_files = []
            data_files.clear_all()
            files_to_add = data_files.lifetime_converter(file)
            data_files.info = info
        elif type == 3:     # Multi column ascii lifetime file
            new_files = []
            data_files.clear_all()
            files_to_add = data_files.lifetime_converter(file)
    old_files = data_files.get_list_of_all_files()

    if new_files != []:
        for file in new_files:
            if file not in old_files:
                files_to_add.append(file)
        file_list = data_files.add_files_from_list(files_to_add)
        data_files.write_data_dir()
        plot_all(file_list)
    else:
        file_list = data_files.add_files_from_list(files_to_add)
        data_files.write_data_dir()
        plot_all(file_list)


def load_x_values():
    root.filename = filedialog.askopenfilename(initialdir=data_files.get_data_dir(), title='Select files', filetype=(('All files', '*.*'), ('Text files', ('*.asc', '*.txt'))))
    xdata = trlfs.load_wavelength_data(root.filename)
    if len(xdata) != len(data_files.integral_xdata):
        print('Error, number of data points and x-values unequal')
        return
    data_files.integral_xdata = xdata
    integral_plotter(data_files.integral_xdata, data_files.integral_ydata)


# Terminates the program if user closes GUI window
def on_closing():
    for dir in data_files.temp_dirs:
        shutil.rmtree(dir, ignore_errors=True)
    root.quit()

data_files = DataFiles()
background = BackgroundCorrection()

# GUI definition
root = tk.Tk()
root.title('TRLFS Analysis Tool')

for x in range(15):
    root.columnconfigure(x, weight=1)
for y in range(22):
    root.rowconfigure(y, weight=1)

####################################################################
# Canvas definition
###################################################################
fi, ax = plt.subplots()

canvas = FigureCanvasTkAgg(fi)
canvas.get_tk_widget().grid(row=2, column=5, rowspan=20, columnspan=10, sticky='ewns')

####################################################################
# Make Toolbar above figure
###################################################################
toolbarFrame = tk.Frame(root)
toolbarFrame.grid(row=1, column=5, columnspan=7, sticky='w')
toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)

####################################################################
# List of data definition
###################################################################
button_plot = tk.Button(root, text='Plot', command= lambda: re_plotter(ax))
button_clear_list = tk.Button(root, text='Clear List', command=data_files.clear_all)
scrollbar = tk.Scrollbar(root, orient='vertical')
lister = tk.Listbox(root, selectmode=tk.EXTENDED, width=25, height=24, yscrollcommand=scrollbar.set)
scrollbar.config(command=lister.yview)
button_integrate = tk.Button(root, text='Integrate', command=data_files.integrate_data)
button_integrate_all = tk.Button(root, text='Integrate all', command=data_files.integrate_all_data)

button_plot.grid(row=2, column=0, columnspan=4, sticky='nsew')
button_clear_list.grid(row=3, column=0, columnspan=4, sticky='wnse')
scrollbar.grid(row=4, column=4, rowspan=16, sticky='wns')
lister.grid(row=4, column=0, columnspan=4, rowspan=16, pady=2, sticky='nsew')
button_integrate.grid(row=21, column=0, columnspan=2, sticky='nsew')
button_integrate_all.grid(row=21, column=2, columnspan=2, sticky='nsew')

integration_lines = DraggableLines(ax, canvas)
fi.canvas.mpl_connect('button_press_event', integration_lines.click)
fi.canvas.mpl_connect('motion_notify_event', integration_lines.drag)
fi.canvas.mpl_connect('button_release_event', integration_lines.release)

####################################################################
# Top menu bar definition
###################################################################
button_open = tk.Button(root, text='Open', command=browse_files)
button_open.grid(row=0, column=0, columnspan=1, sticky='ewns')

button_baseline = tk.Button(root, text='Baseline subtraction', command=background.base_line_correction_init)
button_baseline.grid(row=0, column=1, columnspan=2, sticky='ewns')

button_integrationArea = tk.Button(root, text='Choose integration range', command=integration_lines.init_line)
button_integrationArea.grid(row=0, column=3, columnspan=4, sticky='ewns')

button_save_int_data = tk.Button(root, text='Save integration data', command=data_files.save_integration_data)
button_save_int_data.grid(row=0, column=7, columnspan=2, sticky='ewns')

button_save_selected = tk.Button(root, text='Save selected', command=data_files.save_selected_files)
button_save_selected.grid(row=0, column=12, columnspan=2, sticky='wnse')

button_save_data = tk.Button(root, text='Save all', command=data_files.save_all_files)
button_save_data.grid(row=0, column=14, columnspan=1, sticky='wnse')

button_excitation_wavelengths = tk.Button(root, text="Load x-values", state=tk.DISABLED, command=load_x_values)
button_excitation_wavelengths.grid(row=1, column=12, columnspan=2, sticky='wnse')

button_meta_data = tk.Button(root, text='Acquisition info', command=data_files.get_acquisition_info)
button_meta_data.grid(row=1, column=14, sticky='wnse')

root.bind('<Double-Button-1>', lambda event, ax=ax: short_cut_re_plotter(event, ax))
root.bind('<Return>', lambda event, ax=ax: short_cut_re_plotter(event, ax))
root.bind('<Delete>', data_files.short_cut_remove_file)

# Activated when user closes GUI window
root.protocol("WM_DELETE_WINDOW", on_closing)
# Start the GUI
root.mainloop()