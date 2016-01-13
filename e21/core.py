# -*- coding: utf-8 -*-
#
# e21, (c) 2013, see AUTHORS. Licensed under the GNU GPL.
import collections
import matplotlib as mpl
import matplotlib.ticker
import matplotlib.pyplot as pl
#import seaborn as sns
import e21.core

# colors taken from bootstrap
COLORS = {
    'black': '#000000',
    'darker gray': '#222222',
    'dark gray': '#333333',  #text color
    'gray': '#555555',
    'light gray': '#999999',
    'lighter gray': '#eeeeee',
    'white': '#ffffff',
    'blue': '#049cdb',
    'dark blue': '#0064cd',
    'light blue': '#3a87ad',
    'lighter blue': '#d9edf7',
    'green': '#46a546',
    'light green': '#468847',
    'lighter green': '#dff0d8',
    'red': '#9d261d',
    'light red': '#b94a48',
    'lighter red': '#f2dede',
    'yellow': '#ffc40d',
    'light yellow': '#fcf8e3',
    'orange': '#f89406',
    'pink': '#c3325f',
    'purple': '#7a43b6',
}
RC = {
    'lines.linewidth': 2,
    'lines.linewidth': 2,
    'axes.edgecolor': COLORS['dark gray'],
    'axes.labelsize': 22,
    'axes.linewidth': 2,
    'grid.color': COLORS['dark gray'],
    'grid.linewidth': 1,
    'xtick.major.width': 1.5,
    'xtick.color': COLORS['dark gray'],
    'xtick.labelsize': 22,
    'ytick.major.width': 1.5,
    'ytick.labelsize': 22,
    'ytick.color': COLORS['dark gray'],
    'axes.color_cycle': (COLORS['light blue'], COLORS['green'], COLORS['light red'], COLORS['yellow'], COLORS['pink'], COLORS['black']),
    'axes.labelcolor': COLORS['dark gray'],
    'axes.grid': True,
    'legend.fontsize': 22,
    'figure.figsize': (16, 10),
    'legend.fancybox': True, 
}


def init_mpl(nbins=4, steps=[1, 2, 3, 5, 10]):
    mpl.rcParams.update(RC)
    def init(this):
        super(matplotlib.ticker.AutoLocator, this).__init__(
              nbins=nbins,
              steps=steps)
    matplotlib.ticker.AutoLocator.__init__ = init
    #sns.set_style("whitegrid",RC)
    #sns.set(font="Arial")

#-Measurement-classes----------------------------------------------------------
class Measurement(object):
    """Measurement base class."""
    def __init__(self, data, params):
        self.data = data
        self.params = params

    def __getattr__(self, item):
        return self.data[item]

    def __dir__(self):
        return self.keys()


    def _repr_html_(self):
        return "<h1>" + self.text + "</h1>"

    
    def __len__(self):
        key = self.data.keys()[0]
        return len(self.data[key])
    

class Experiment(object):
    """Minimal implementation of a experiment class."""
    def __init__(self, measurements={}):
        self._measurements = {}
        
    def __getitem__(self, item):
        return self._measurements[item]
    
    def __dir__(self):
        return self.keys()

    def __setitem__(self, item, value):
        self._measurements[item] = value

    def __delitem__(self, item):
        del(self._measurements[item])

    def __len__(self):
        return len(self._measurements)
    
#-Plottable-classes------------------------------------------------------------
class Plottable(object):
    """Mixin class providing a generic plot function."""
    def plot(self, y, x, axes=None, subplot_kw={}, fig_kw={}, **kw):
        """
        :returns: A matplotlib figure and axes instance.
        """
        x, y = lookup(self, x), lookup(self, y)
        if not axes:
            f, axes = pl.subplots(1, 1, subplot_kw=subplot_kw, **fig_kw)
        axes.plot(x, y, **kw)
        pl.legend()
        return axes.figure, axes


def lookup(obj, name):
    """Helper function to provide function and string lookup of attributes.

    E.g.::
        >>>class Test(object):
        ...    pass
        >>>t = Test()
        >>>t.value = 'myvalue'
        >>>lookup(t, 'value')
        myvalue
        >>>lookup(t, lambda x: x.value)
        myvalue
        >>lookup(t, [1, 2, 3])
        [1, 2, 3]
    """
    if isinstance(name, basestring):
        return getattr(obj, name)
    elif isinstance(name, collections.Callable):
        return name(obj)
    else:
        return name

init_mpl()
