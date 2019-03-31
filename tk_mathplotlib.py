#import serial
from tkinter import *
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib import style
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import serial


root = Tk()
root.geometry('1200x700+200+100')
root.title('This is my root window')
root.state('zoomed')
root.config(background='#fafafa')

xar = []
yar = []
zar = []
war = []
cnt = 1

style.use('ggplot')
fig = plt.figure(figsize=(14, 4.5), dpi=100)
ax1 = fig.add_subplot(1, 1, 1)
ax1.set_ylim(0, 8140)
line, = ax1.plot(xar, yar ,'r', marker='o')
line2, =ax1.plot(xar,zar,'b',marker='o')
line3, =ax1.plot(xar,war,'y',marker='o')

ser = serial.Serial('com4', 115200)

def animate(i):
    # ser.reset_input_buffer()
    data = ser.readline().decode("utf-8")
    data_array = data.split(',')
    print (data_array)
    try :
        yvalue = int(data_array[1])
        xvalue = int(data_array[3])
        zvalue = int (data_array[2])
        
    except :
        print ("ya embrouille")
        pass
    else :
        yar.append(yvalue)
        zar.append(xvalue)        
        war.append(zvalue)        
        xar.append(i)
        line.set_data(xar, yar)
        line2.set_data(xar, zar)
        line3.set_data(xar, war)
        ax1.set_xlim(i-100, i+1)
        i=i+1
        if(i>100):                            #If you have 50 or more points, delete the first one from the array
            yar.pop(0)                       #This allows us to just see the last 50 data points
            xar.pop(0)
            zar.pop(0)
            war.pop(0)
    
plotcanvas = FigureCanvasTkAgg(fig, root)
plotcanvas.get_tk_widget().grid(column=1, row=1)
ani = animation.FuncAnimation(fig, animate, interval=1, blit=False)

root.mainloop()