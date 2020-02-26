# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 15:46:56 2019

@author: ash
"""

import tkinter

from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure

import numpy as np
import serial
import time
import math
# sphinx_gallery_thumbnail_number = 3
import matplotlib.pyplot as plt

ser = serial.Serial('com4', 115200)
file = open ("calib_bal", "w")
log = []
for i in range (1000):
        ser.readline().decode("utf-8")
temp_acquisition = 30
time_start = time.time()

for i in range (1000):
        ser.readline().decode("utf-8")
        
while time.time() < time_start+temp_acquisition :
    data = ser.readline().decode("utf-8")
    log.append(int(data[:-1].split(",")[2]))
    file.write(data[:-1])
ser.close()

root = tkinter.Tk()
root.wm_title("Embedding in Tk")


fig = Figure(figsize=(5, 4), dpi=100)
t = np.arange(0, 3, .01)
fig.add_subplot(111).plot(log)


canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.
canvas.draw()
canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)


toolbar = NavigationToolbar2Tk(canvas, root)
toolbar.update()
canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)


def on_key_press(event):
    print("you pressed {}".format(event.key))
    key_press_handler(event, canvas, toolbar)


canvas.mpl_connect("key_press_event", on_key_press)


def _quit():
    root.quit()     # stops mainloop
    root.destroy()  # this is necessary on Windows to prevent

button = tkinter.Button(master=root, text="Quit", command=_quit)
button.pack(side=tkinter.BOTTOM)

tkinter.mainloop()
# If you put root.destroy() here, it will cause an error if the window is
# closed with the window manager.