# -*- coding: utf-8 -*-
#
# e21, (c) 2013, see AUTHORS. Licensed under the GNU GPL.
import matplotlib as mpl
import matplotlib.ticker
import matplotlib.pyplot as pl

from e21.utility import COLORS

#-Measurement-classes----------------------------------------------------------
class Measurement(object):
    """Measurement base class."""
    def __init__(self, data, params):
        self.data = data
        self.params = params

    def __getattr__(self, item):
        return self.data[item]

#-Plottable-classes------------------------------------------------------------
class Plottable(object):
    """Mixin class providing a generic plot function."""
    RC = {
        'ticker.nbins': 4, #number of ticker bins
        'ticker.steps': [1, 2, 3, 5, 10],
        'lines.linewidth': 2,
        'lines.linewidth': 2,
        'axes.edgecolor': COLORS['dark gray'],
        'axes.labelsize': 14,
        'axes.linewidth': 2,
        'grid.color': COLORS['dark gray'],
        'grid.linewidth': 1,
        'xtick.major.width': 1.5,
        'xtick.color': COLORS['dark gray'],
        'xtick.labelsize': 14,
        'ytick.major.width': 1.5,
        'ytick.labelsize': 14,
        'ytick.color': COLORS['dark gray'],
        'axes.color_cycle': (COLORS['light blue'], COLORS['green'], COLORS['light red'], COLORS['yellow'], COLORS['pink'], COLORS['black']),
        'axes.labelcolor': COLORS['dark gray'],
        'axes.grid': True,
        'legend.fontsize': 14,
        'figure.figsize': (8, 5)
    }

    def plot(self, y, x, axes=None, subplot_kw={}, fig_kw={}, **kw):
        """
        :returns: A matplotlib figure and axes instance.
        """
        x, y = lookup(x), lookup(y)
        if not axes:
            # monkey patch AutoLocator
            old_init = matplotlib.ticker.AutoLocator.__init__
            def init(this):
                super(AutoLocator, this).__init__(
                        this,
                        nbins=self.RC['ticker.nbins'],
                        steps=self.RC['ticker.steps'])
            matplotlib.ticker.AutoLocator.__init__ = init
            with mpl.rc_context(rc=self.RC):
                f, axes = pl.subplots(1, 1, subplot_kw=subplot_kw, **fig_kw)
                axes.plot(x, y, **kw)
            matplotlib.ticker.AutoLocator.__init__ = old_init # undo monkeypatch
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
        return obj
