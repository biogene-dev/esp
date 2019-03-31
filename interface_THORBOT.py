# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 09:20:39 2019

@author: ash
"""

# -*- coding: utf-8 -*-
 
from tkinter import * 

def awesome():
    def click_me():
        print(var.get())
        return

    root = Toplevel()
    root.title("a good try")
    var = IntVar()
    x = Checkbutton(root, text = "check me", variable = var)
    y = Button(root, text = "click me", command = click_me)
    x.pack()
    y.pack()
    root.mainloop()
    return



fenetre = Tk()

#NOM DE LA FENETRE
label = Label(fenetre, text="Hello World")
label.pack()

#BOUTON QUIT
bouton=Button(fenetre, text="Fermer", command=fenetre.quit)
bouton.pack()

# checkbutton
bouton = Checkbutton(fenetre, text="TEMPS REEL", command = awesome)
bouton.pack()

fenetre.mainloop()