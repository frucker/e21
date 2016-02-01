import tkMessageBox
import pickle
import tkMessageBox
from Tkinter import *
from ttk import *

def save_obj(obj, name):
    with open('.'+name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open('.'+name + '.pkl', 'rb') as f:
        return pickle.load(f)

def test():
    print 'Test'



    
def runwidget(params):
    Variables = {}
    test = Tk()
    test.title('Fuckit')
    
    t = ButtonWidget(test)
    cols = 0
    row = 1
    
    for i, key in enumerate(params.keys()):
        if isinstance(params[key], dict):
            row2 = 1
            cols = cols +1 
            for j, key2 in enumerate(sorted(params[key].keys())):
                Variables[key] = {}  
                if isinstance(key2, dict):
                    Variables = []
                else:
                    row2 = row2 +1
                    Label(t, text = key, font=("Helvetica", 16)).grid(row=1, column = cols)
                    Variable = IntVar()
                    C1 = Checkbutton(t, text = key2, variable = Variable,
                                                 onvalue = 1, offvalue = 0)
                    C1.grid(row = row2, column = cols , sticky = 'w')
                    Variables[key][key2] = Variable
        else:
            row = row +1
            Variable = IntVar()
            C1 = Checkbutton(t, text = key, variable = Variable, 
                                         onvalue = 1, offvalue = 0)
            C1.grid(row = row, column = 0 , sticky = 'w')
            Variables[key] = Variable
    Label(t, text = 'Toplevel', font=("Helvetica", 16)).grid(row=1, column = 0)
    test.mainloop()
    return Variables
    
class ButtonWidget(Frame):
    def __init__(self, parent):
        Frame.__init__(self, parent)
        B = Button(self, text ="Done!", command = self._close).grid()
        self.grid()
        #scrollbar = Scrollbar(self)
        #scrollbar.pack(side=RIGHT, fill=Y)
    def _close(self):
        self.master.destroy()
        
class SelectionWidget(Frame):  
    
    def __init__(self, parent, params, name):
        Frame.__init__(self, parent)
        self.params = params
        self.name = name
        
        
        self.x = 0
        flashCardText = Label(self, text = 'Please Choose', font=("Helvetica", 16)).grid(row=1, column = 1)
        B1 = Button(self, text ="New!", command = self.new).grid(row=2, column = 0)
        B2 = Button(self, text ="Old!", command = self.old).grid(row=2, column = 1)
        B3 = Button(self, text ="Edit Old!", command = self.edit_old).grid(row=2, column = 2)
        self.grid()
    def new(self):
        self._close()
        runwidget(self.params)
    def old(self):
        self._close()
        runwidget(self.params)
    def edit_old(self):
        self._close()
        runwidget(self.params)
        
        
        
    def _close(self):
        self.master.destroy()
        

def overview(params, name):
    
    try:
        Variables = load_obj(name)
        test = Tk()
        test.title('Fuckit')
        
        t = SelectionWidget(test, params, name)
        test.mainloop()
    except IOError:
        Variables = runwidget(params) 
    for keys in Variables.keys():
        Variables[keys] = Variables[keys].get()
    return Variables




    
