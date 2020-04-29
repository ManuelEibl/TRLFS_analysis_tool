"""
This program reads spectroscopic x,y data (.asc, .txt, or AndorFile .sif) of single or multiple spectra.


Start: 10.09.2019
Last edit: 28.04.2020
Creator: Manuel Eibl
"""

import os
import tkinter as tk
from tkinter import filedialog
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
matplotlib.use('TkAgg')
import numpy as np
from scipy import interpolate, integrate
import ctypes
from sif_reader import extract_calibration
import trlfs
import shutil

# Handles all the data input, access, integration and output
class DataFiles():
    def __init__(self, dir='/', temp_dirs=[] , file_list=[], file_dict={}, integral_xdata=[], integral_ydata=[]):
        self.dir = dir
        self.temp_dirs = temp_dirs
        self.file_list = file_list  # list of all files with exact paths
        self.file_dict = file_dict  # dict, key_word = file name, value = path
        self.integral_xdata = integral_xdata
        self.integral_ydata = integral_ydata

        #load last opened directory for quicker file search
        saved_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'log.dat')
        if os.path.exists(saved_path):
            with open(saved_path, 'r') as r:
                content = r.readlines()
            try:
                self.dir = content[0]
            except:
                pass

    def add_files(self, new_files=[]):
        for item in new_files:
            self.file_list.append(item)
            self.dir, key_word = os.path.split(item)
            self.file_dict[key_word] = item
        return self.file_list

    def get_all_files(self):
        return self.file_list

    def get_all_keys(self):
        return self.file_dict.keys()

    def get_files(self, key_word):
        return self.file_dict[key_word]

    def get_files_from_index(self, index):
        selected_files = []
        try:
            for value in index:
                selected_files.append(self.file_list[value])
        except:
            selected_files.append(self.file_list[index])
        return selected_files

    def get_dir(self):
        return self.dir

    # Write last used directory for initialization above
    def write_dir(self):
        path = os.path.dirname(os.path.abspath('__file__'))
        log_file_name = 'log.dat'
        log_path = os.path.join(path, log_file_name)
        if os.path.split(self.dir)[1] == '.temp':
            self.dir = os.path.split(self.dir)[0]
        with open (log_path, 'w+') as f:
            f.write(self.dir)

    def short_cut_remove_file(self, event):
        selected_index = lister.curselection()
        for item in selected_index[::-1]:
            key_word = lister.get(item)
            file = self.get_files(key_word)
            del self.file_dict[key_word]
            for i, value in enumerate(self.file_list):
                if file == value:
                    self.file_list.pop(i)
        plotter(self.file_list)

    def clear_list(self):
        self.file_list = []
        self.file_dict = {}
        self.integral_xdata = []
        self.integral_ydata = []
        lister.delete(0, tk.END)

    def data_integrate(self):
        self.integral_xdata = []
        self.integral_ydata = []
        x1 = None
        x2 = None
        if line1.get_range() is not None:
            x1, x2 = line1.get_range()
        items = lister.curselection()
        for i, item in enumerate(items):
            select = lister.get(item)
            file = self.get_files(select)
            try:
                xvalue = float(file.split('.')[0].split('_')[-1])
            except:
                xvalue = i + 1
            self.integral_xdata.append(xvalue)
            if background.bg_corr_data == []:
                xdata, ydata = trlfs.load_data(file)
            else:
                xdata, ydata = background.bg_corr_data[i]
            if x1 != None:
                for j, value in enumerate(xdata):
                    if value >= x1:
                        low_value = j
                        break
                for k, value in enumerate(xdata):
                    if value >= x2:
                        upper_value = k
                        break
            else:
                low_value = 0
                upper_value = -1
            interpolate_function = interpolate.interp1d(xdata[low_value:upper_value],
                                                        ydata[low_value:upper_value])
            integral = integrate.quad(interpolate_function, xdata[low_value], xdata[upper_value - 1], limit=50)[0]
            self.integral_ydata.append(integral)
        integral_plotter(self.integral_xdata, self.integral_ydata)

    def data_integrate_all(self):
        self.integral_xdata = []
        self.integral_ydata = []
        x1 = 0.0
        x2 = 0.0
        if line1.get_range() is not None:
            x1, x2 = line1.get_range()
        items = lister.get(0, tk.END)
        for i, item in enumerate(items):
            select = lister.get(i)
            file = data_files.get_files(select)
            file_name = os.path.split(file)[1].split('.')[0]

            try:
                if file_name.split('_')[-1] == 'corr':
                    first_value = int(file_name.split('_')[-3])
                    if first_value > 30000:
                        xvalue = float(file_name.split('_')[-3]) / 100
                    else:
                        xvalue = int(file_name.split('_')[-3])
                else:
                    first_value = int(file_name.split('_')[-1])
                    if first_value > 30000:
                        xvalue = float(file_name.split('_')[-1]) / 100
                    else:
                        xvalue = int(file_name.split('_')[-1])
            except:
                xvalue = i + 1
            self.integral_xdata.append(xvalue)
            if background.bg_corr_data == []:
                xdata, ydata = trlfs.load_data(file)
            else:
                xdata, ydata = background.bg_corr_data[i]
            if x1 + x2 != 0.0:
                for j, value in enumerate(xdata):
                    if value >= x1:
                        low_value = j
                        break
                for k, value in enumerate(xdata):
                    if value >= x2:
                        upper_value = k
                        break
            else:
                low_value = 0
                upper_value = -1
            interpolate_function = interpolate.interp1d(xdata[low_value:upper_value],
                                                        ydata[low_value:upper_value])
            integral = integrate.quad(interpolate_function, xdata[low_value], xdata[upper_value - 1], limit=50)[0]
            self.integral_ydata.append(integral)
        integral_plotter(self.integral_xdata, self.integral_ydata)

    def output_maker(self):
        output = ''
        for i, item in enumerate(self.integral_xdata):
            output += str(item) + '\t' + str(self.integral_ydata[i]) + '\n'
        return output

    def write_data(self):
        if line1.get_range() is not None:
            x1, x2 = line1.get_range()
        else:
            item = lister.get(0)
            file = data_files.get_files(item)
            xdata, _ = trlfs.load_data(file)
            x1 = xdata[0]
            x2 = xdata[-1]
        output = self.output_maker()
        init_file_name = f'Integration from {int(x1*100)/100} - {int(x2*100)/100}'
        abs_path = filedialog.asksaveasfilename(initialdir=self.get_dir(), initialfile=init_file_name, defaultextension='.txt', filetypes=[('Text file', '*.txt'), ('Others', ('*.asc', '*dat'))])
        with open(str(abs_path), 'w+') as o:
            o.write(str(output))

    def lifetime_converter(self, file):
        out_file_list = []
        file_path, file_name = os.path.split(file)
        temp_path = os.path.join(file_path, '.temp')
        if not os.path.exists(temp_path):
            os.mkdir(temp_path)
            ctypes.windll.kernel32.SetFileAttributesW(temp_path, 2)
            self.temp_dirs.append(temp_path)
        data, info = trlfs.lifetime_handler(file)
        xdata = extract_calibration(info)
        for i in range(info['NumberOfFrames']):
            temp_file_name = file_name.split('.')[0] + '_' + str(int(1 + (i*info['GateDelayStep']*10**-6))) + '.asc'
            temp_file_path = os.path.join(temp_path, temp_file_name)
            out = ''
            ydata = data[i][0]
            for j,value in enumerate(xdata):
                out += str(value) + '\t' + str(ydata[j]) + '\n'
            with open(temp_file_path, 'w+') as f:
                f.write(out)
            out_file_list.append(temp_file_path)

        return out_file_list

    def save_selected(self):
        out_file_list = []
        selected_items = lister.curselection()
        for item in selected_items:
            select = lister.get(item)
            file = data_files.get_files(select)
            file_path, file_name = os.path.split(file)
            if os.path.split(file_path)[1] == '.temp':
                temp_path = file_path
                temp_file_path = os.path.join(temp_path, file_name)
            else:
                temp_path = os.path.join(file_path, '.temp')
                temp_file_path = os.path.join(temp_path, file_name)
                out_file_list.append(temp_file_path)
                if not os.path.exists(temp_path):
                    os.mkdir(temp_path)
                    ctypes.windll.kernel32.SetFileAttributesW(temp_path, 2)
                    self.temp_dirs.append(temp_path)
            out = ''
            xdata, ydata = trlfs.load_data(file)
            for i, value in enumerate(xdata):
                out += str(value) + '\t' + str(ydata[i]) + '\n'
            with open(temp_file_path, 'w+') as o:
                o.write(out)

        abs_path = filedialog.askdirectory(initialdir=self.get_dir())
        for temp_file in out_file_list:
            _, temp_filename = os.path.split(temp_file)
            abs_out_path = os.path.join(abs_path, temp_filename)
            shutil.copyfile(temp_file, abs_out_path)

    def save_data(self):
        out_file_list = []
        all_files = self.get_all_files()
        for file in all_files:
            file_path, file_name = os.path.split(file)
            if os.path.split(file_path)[1] == '.temp':
                temp_path = file_path
                temp_file_path = os.path.join(temp_path, file_name)
            else:
                temp_path = os.path.join(file_path, '.temp')
                temp_file_path = os.path.join(temp_path, file_name)
                out_file_list.append(temp_file_path)
                if not os.path.exists(temp_path):
                    os.mkdir(temp_path)
                    ctypes.windll.kernel32.SetFileAttributesW(temp_path, 2)
                    self.temp_dirs.append(temp_path)
            out = ''
            xdata, ydata = trlfs.load_data(file)
            for i, value in enumerate(xdata):
                out += str(value) + '\t' + str(ydata[i]) + '\n'
            with open(temp_file_path, 'w+') as o:
                o.write(out)

        abs_path = filedialog.askdirectory(initialdir=self.get_dir())
        for temp_file in out_file_list:
            _, temp_filename = os.path.split(temp_file)
            abs_out_path = os.path.join(abs_path, temp_filename)
            shutil.copyfile(temp_file, abs_out_path)

# This class creates draggable lines to define the integration area
class DraggyLine():
    lock = None

    def __init__(self, ax, canvas, lines=[]):
        self.ax = ax
        self.canvas = canvas
        self.lines = lines

    def init_line(self):
        # Make default lines first
        file = data_files.get_all_files()[0]
        xdata, _ = trlfs.load_data(file)
        center_index = int(len(xdata) / 2)
        x1 = xdata[center_index] - 5
        x2 = xdata[center_index] + 5

        left_line = self.ax.axvline(x1)
        right_line = self.ax.axvline(x2)
        self.lines.append(left_line)
        self.lines.append(right_line)
        self.canvas.draw()

    def click(self, event):
        w, _ = self.canvas.get_width_height()
        xlim_lower, xlim_upper = self.ax.get_xlim()
        pix_per_nm = w / (xlim_upper - xlim_lower)
        for line in self.lines:
            try:
                len(line.get_xdata())
                line_xdata = line.get_xdata()[0]
            except:
                line_xdata = line.get_xdata()
            if np.abs(((line_xdata - event.xdata)*pix_per_nm)) < 10.0:
                DraggyLine.lock = line

    def release(self, event):
        if DraggyLine.lock is None:
            return
        else:
            DraggyLine.lock = None

    def drag(self, event):
        if DraggyLine.lock is None:
            return
        chosen_line = DraggyLine.lock
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
class DraggyPoint():
    lock = None
    def __init__(self, ax, canvas, items, x_coord=0, dots=None, lines=[], bg_lines=None):
        self.ax = ax
        self.canvas = canvas
        self.items = items
        self.x_coord = x_coord
        self.dots = dots
        if self.dots is None:
            self.dots = []
        self.lines = lines
        self.bg_lines = bg_lines
        self.init_items()

    def init_items(self):
        data_file = data_files.get_files_from_index(self.items[0])[0]
        self.xdata, self.ydata = trlfs.load_data(data_file)

    def draw_bg_line(self):
        if len(self.dots) > 1:
            self.dots_xdata = []
            self.dots_ydata = []
            for dot in self.dots:
                try:
                    self.dots_xdata.append(dot.get_xdata()[0])
                    self.dots_ydata.append(dot.get_ydata()[0])
                except:
                    self.dots_xdata.append(dot.get_xdata())
                    self.dots_ydata.append(dot.get_ydata())
            sorted_dots_data = sorted(zip(self.dots_xdata, self.dots_ydata))
            sorted_xdata = []
            sorted_ydata = []
            for data in sorted_dots_data:
                sorted_xdata.append(data[0])
                sorted_ydata.append(data[1])
            if self.bg_lines is None:
                self.bg_lines = self.ax.plot(sorted_xdata, sorted_ydata)
            else:
                self.bg_lines[0].set_data(sorted_xdata, sorted_ydata)

    def get_closest_point(self, x_coord):
        closest_point = {}
        for i, item in enumerate(self.xdata):
            diff = np.abs(item - x_coord)
            if i == 0:
                closest_point[0] = (i, diff)
            else:
                if diff < closest_point[0][1]:
                    closest_point[0] = (i, diff)
        closest_index = closest_point[0][0]
        closest_point_data = [self.xdata[closest_index], self.ydata[closest_index]]
        return closest_point_data

    def click(self, event):
        if event.dblclick:
            self.x_coord = event.xdata
            closest_point_data = self.get_closest_point(self.x_coord)
            dot = self.ax.plot([closest_point_data[0]], [closest_point_data[1]], marker='o', color='red', markersize=5, linestyle='none', zorder=5, picker=5)[0]
            self.dots.append(dot)
            self.draw_bg_line()
            self.canvas.draw()
        else:
            w, h = canvas.get_width_height()
            xlim_lower, xlim_upper = ax.get_xlim()
            pix_per_nm = w / (xlim_upper - xlim_lower)
            for dot in self.dots:
                try:
                    len(dot.get_xdata())
                    dot_xdata = dot.get_xdata()[0]
                except:
                    dot_xdata = dot.get_xdata()
                if np.abs(((dot_xdata - event.xdata) * pix_per_nm)) < 10.0:
                    DraggyPoint.lock = dot

    def release(self, event):
        if DraggyPoint.lock is None:
            return
        else:
            DraggyPoint.lock = None
            self.draw_bg_line()

    def select(self, event):
        checker = str(event.mouseevent.button)
        if checker == '3':
            select_x_coord = event.artist.get_xdata()[0]
            for i, dot in enumerate(self.dots):
                if dot.get_xdata()[0] == select_x_coord:
                    dot.remove()
                    self.dots.pop(i)
            self.draw_bg_line()
            self.canvas.draw()

    def drag(self, event):
        if DraggyPoint.lock is None:
            return
        chosen_dot = DraggyPoint.lock
        closest_point_data = self.get_closest_point(event.xdata)
        chosen_dot.set_xdata(closest_point_data[0])
        chosen_dot.set_ydata(closest_point_data[1])
        self.draw_bg_line()
        self.canvas.draw()

    def get_bg_data(self):
        return self.dots


# Optional feature of correcting linear background with n-1 lines between n points
class BackgroundCorrection():
    def __init__(self, bg_corr_data=[]):
        self.bg_corr_data = bg_corr_data

    # BG-correction window initialization
    def base_line_correction_init(self):
        self.fi2, self.ax2 = plt.subplots()
        self.correction_window = tk.Toplevel(root)
        self.correction_window.title('Background correction')
        self.correction_canvas = FigureCanvasTkAgg(self.fi2, master=self.correction_window)
        self.bg_corr_data = []

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

        button_bg_corr = tk.Button(self.correction_window, text='Apply correction to all', command=self.correct_all)
        button_bg_corr.grid(row=0, column=0, sticky='nsew')

        self.data_plotter()
        self.draw_points()

        self.correction_canvas.draw()
        self.correction_window.mainloop()

    def data_plotter(self):
        self.items = lister.curselection()
        for item in self.items:
            select = lister.get(item)
            file = data_files.get_files(select)
            xdata, ydata = trlfs.load_data(file)
            self.ax2.plot(xdata, ydata)

    def draw_points(self):
        self.u_line1 = DraggyPoint(self.ax2, self.correction_canvas, self.items)
        self.cid_click = self.fi2.canvas.mpl_connect('button_press_event', self.u_line1.click)
        self.cid_select = self.fi2.canvas.mpl_connect('pick_event', self.u_line1.select)
        self.cid_drag = self.fi2.canvas.mpl_connect('motion_notify_event', self.u_line1.drag)
        self.cid_release = self.fi2.canvas.mpl_connect('button_release_event', self.u_line1.release)

        self.correction_canvas.draw()

    def correct_all(self):
        self.disconnect()
        dots_xdata = self.u_line1.dots_xdata
        dots_ydata = self.u_line1.dots_ydata
        bg_functions = {}
        for i, dot_xdata in enumerate(dots_xdata):
            if i > 0:
                slope = (dots_ydata[i] - dots_ydata[i-1])/(dots_xdata[i] - dots_xdata[i-1])
                b = dots_ydata[i-1] - slope*dots_xdata[i-1]
                bg_functions[dot_xdata] = [slope, b]
        ranges = list(bg_functions.keys())
        self.all_files = data_files.get_all_files()
        for file in self.all_files:
            xdata, ydata = trlfs.load_data(file)
            corrected_ydata = []
            points = list(range(0, len(ranges)))
            for i, value in enumerate(xdata):
                j = 0
                for j in points:
                    if value <= ranges[j]:
                        slope, b = bg_functions[ranges[j]]
                        corrected_ydata.append(ydata[i] - (slope*xdata[i] + b))
                        break
                    elif value > ranges[-1]:
                        slope, b = bg_functions[ranges[-1]]
                        corrected_ydata.append(ydata[i] - (slope*xdata[i] + b))
                        break
                    else:
                        j += 1

            self.bg_corr_data.append([xdata, corrected_ydata])
        ax.cla()
        for i, data in enumerate(self.bg_corr_data):
            ax.plot(data[0], data[1])
        self.correction_window.destroy()
        self.save_bg_corr_data()
        canvas.draw()

    # Needs to be stored in /.temp/ so it can be replotted at any time
    def save_bg_corr_data(self):
        temp_files_list = []
        temp_files_dict = {}
        data_path = data_files.get_dir()
        if os.path.split(data_path)[1] == '.temp':
            temp_path = data_path
        else:
            temp_path = os.path.join(data_path, '.temp')
            if not os.path.exists(temp_path):
                os.mkdir(temp_path)
                ctypes.windll.kernel32.SetFileAttributesW(temp_path, 2)
                data_files.temp_dirs.append(temp_path)
        for i, data in enumerate(self.bg_corr_data):
            lister_name = list(data_files.get_all_keys())[i]
            temp_file = os.path.split(self.all_files[i])[1]
            if temp_file.split('_')[-1].split('.')[0] == 'corr':
                temp_file_name = temp_file
            else:
                temp_file_name = temp_file.split('.')[0] + '_bg_corr.asc'
            temp_file_path = os.path.join(temp_path, temp_file_name)
            out = ''
            for j, value in enumerate(data[0]):
                out += str(value) + '\t' + str(data[1][j]) + '\n'
            with open(temp_file_path, 'w+') as o:
                o.write(out)
            temp_files_list.append(temp_file_path)
            temp_files_dict[lister_name] = temp_file_path
        data_files.file_list = temp_files_list
        data_files.file_dict = temp_files_dict

    def disconnect(self):
        self.fi2.canvas.mpl_disconnect(self.cid_click)
        self.fi2.canvas.mpl_disconnect(self.cid_drag)
        self.fi2.canvas.mpl_disconnect(self.cid_release)
        self.fi2.canvas.mpl_disconnect(self.cid_select)


# Adds all data to the lister box and plots it
def plotter(file_list):
    ax.cla()
    lister.delete(0, tk.END)
    for file in file_list:
        _, file_name = os.path.split(file)
        lister.insert(tk.END, file_name)
        xdata, ydata = trlfs.load_data(file)
        ax.plot(xdata, ydata, linewidth=1)
    canvas.draw()


def re_plotter(ax):
    ax.cla()
    line1.x1 = None
    line1.x2 = None
    items = lister.curselection()
    for item in items:
        select = lister.get(item)
        file = data_files.get_files(select)
        xdata, ydata = trlfs.load_data(file)
        ax.plot(xdata, ydata)
    canvas.draw()


# For double click re-plotting
def short_cut_re_plotter(event, ax):
    ax.cla()
    line1.x1 = None
    line1.x2 = None
    items = lister.curselection()
    for item in items:
        select = lister.get(item)
        file = data_files.get_files(select)
        xdata, ydata = trlfs.load_data(file)
        ax.plot(xdata, ydata)
    canvas.draw()


def integral_plotter(integral_xdata, integral_ydata):
    ax.cla()
    ax.plot(integral_xdata, integral_ydata, marker='o')
    canvas.draw()


# Opens a selection window and loads files. It checks whether the chosen file is .sif lifetime data, which has to
# be separated into individual files first (done by lifetime_converter()) then loads and plots the data
def browseFiles():
    root.filename = filedialog.askopenfilenames(initialdir=data_files.get_dir(), title='Select files', filetype=(('All files', '*.*'), ('Text files', ('*.asc', '*.txt')), ('Andor .sif files','*.sif')))
    new_files = list(root.filename)
    files_to_add = []
    if len(new_files) == 1:
        file = new_files[0]
        lifetime = trlfs.lifetime_handler(file)
        if lifetime != None:
            new_files = []
            data_files.clear_list()
            files_to_add = data_files.lifetime_converter(file)
    old_files = data_files.get_all_files()
    if new_files != []:
        for file in new_files:
            if file not in old_files:
                files_to_add.append(file)
        file_list = data_files.add_files(files_to_add)
        data_files.write_dir()
        plotter(file_list)
    else:
        file_list = data_files.add_files(files_to_add)
        data_files.write_dir()
        plotter(file_list)


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
#root.geometry('1050x650')

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
button_clear_list = tk.Button(root, text='Clear List', command=data_files.clear_list)
scrollbar = tk.Scrollbar(root, orient='vertical')
lister = tk.Listbox(root, selectmode=tk.EXTENDED, width=25, height=24, yscrollcommand=scrollbar.set)
scrollbar.config(command=lister.yview)
button_integrate = tk.Button(root, text='Integrate', command=data_files.data_integrate)
button_integrate_all = tk.Button(root, text='Integrate all', command=data_files.data_integrate_all)

button_plot.grid(row=2, column=0, columnspan=4, sticky='nsew')
button_clear_list.grid(row=3, column=0, columnspan=4, sticky='wnse')
scrollbar.grid(row=4, column=4, rowspan=16, sticky='wns')
lister.grid(row=4, column=0, columnspan=4, rowspan=16, pady=2, sticky='nsew')
button_integrate.grid(row=21, column=0, columnspan=2, sticky='nsew')
button_integrate_all.grid(row=21, column=2, columnspan=2, sticky='nsew')

line1 = DraggyLine(ax, canvas)
fi.canvas.mpl_connect('button_press_event', line1.click)
fi.canvas.mpl_connect('motion_notify_event', line1.drag)
fi.canvas.mpl_connect('button_release_event', line1.release)

####################################################################
# Top menu bar definition
###################################################################
button_open = tk.Button(root, text='Open', command=browseFiles)
button_open.grid(row=0, column=0, columnspan=1, sticky='ewns')

button_baseline = tk.Button(root, text='Baseline subtraction', command=background.base_line_correction_init)
button_baseline.grid(row=0, column=1, columnspan=2, sticky='ewns')

button_integrationArea = tk.Button(root, text='Choose integration range', command=line1.init_line)
button_integrationArea.grid(row=0, column=3, columnspan=4, sticky='ewns')

button_save_int_data = tk.Button(root, text='Save integration data', command=data_files.write_data)
button_save_int_data.grid(row=0, column=7, columnspan=2, sticky='ewns')

button_save_selected = tk.Button(root, text='Save selected', command=data_files.save_selected)
button_save_selected.grid(row=0, column=12, columnspan=2, sticky='wnse')

button_save_data = tk.Button(root, text='Save all', command=data_files.save_data)
button_save_data.grid(row=0, column=14, columnspan=1, sticky='wnse')

root.bind('<Double-Button-1>', lambda event, ax=ax: short_cut_re_plotter(event, ax))
root.bind('<Return>', lambda event, ax=ax: short_cut_re_plotter(event, ax))
root.bind('<Delete>', data_files.short_cut_remove_file)

# Activated when user closes GUI window
root.protocol("WM_DELETE_WINDOW", on_closing)
# Start the GUI
root.mainloop()