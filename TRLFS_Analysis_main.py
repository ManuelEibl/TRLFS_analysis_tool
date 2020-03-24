"""
So here I would like to tackle a huge project.
Basically I want to translate Bins TRLFS script from MATLAB into python.

Start: 10.09.2019
Creator: Manuel Eibl
"""

# Here I should import all the additional python scripts which I will write as modules
# SIF reading
# Integration
import tkinter as tk
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

# GUI definition
root = tk.Tk()
root.geometry('1050x650')

button_open = tk.Button(text = '1. Open')
button_open.place(x=0, y=25)

button_Baseline = tk.Button(text = '2. Baseline subtraction')
button_Baseline.place(x=54, y=25)

button_IntegrationArea = tk.Button(text = '3. Choose integration range')
button_IntegrationArea.place(x=185, y=25)


f = Figure(figsize=(10,5), dpi=100)
a = f.add_subplot(111)
a.plot([1,2,3,4,5,6,7,8],[5,6,1,3,8,9,3,5])

canvas = FigureCanvasTkAgg(f)
canvas.get_tk_widget().pack(side=tk.TOP, expand=True)

toolbar = NavigationToolbar2Tk(canvas, root)

# Start the GUI
root.mainloop()
