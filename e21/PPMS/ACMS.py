# -*- coding: utf-8 -*-
#
# e21, (c) 2013, see AUTHORS. Licensed under the GNU GPL.
import e21.core
from e21.core import lookup
import numpy as np
from e21.sweet16 import Loader
import quantities as pq


def _correct_field(field, correction):
    return field - correction * pq.T

# Susceptibility-classes ------------------------------------------------------


class ACMS(e21.core.Measurement, e21.core.Plottable):

    @property
    def M(self):
        return self.data['M-DC (emu)']*pq.emu

    @property
    def real(self):
        return self.data["M' (emu)"]*pq.emu

    @property
    def imag(self):
        return self.data["M'' (emu)"]*pq.emu

    @property
    def time(self):
        return self.data['Time Stamp (sec)']*pq.s

    @property
    def field(self):
        return self.data['Magnetic Field (oe)']

    @property
    def mean_field(self):
        return np.mean(self.field)

    @property
    def mean_current(self):
        return np.round(np.mean(self.data['current']), 2)

    @property
    def temperature(self):
        return self.data['Temperature (K)']

    @property
    def mean_temperature(self):
        return np.mean(self.temperature)


class FieldScan(ACMS, e21.core.Plottable):

    def plot(self,  x='field', y1='M', y2='real', y3='imag', axes=None,
             subplot_kw={}, fig_kw={}, **kw):
        """A default implementation for a generic fieldscan plot."""
        xa, ya1 = lookup(self, x), lookup(self, y1)
        ya2, ya3 = lookup(self, y2), lookup(self, y3)
        xlabel = 'B ({0})'.format(xa.dimensionality) if x is 'field' else ''
        ylabel = 'U(V)' if y is 'real' else ''
        subplot_default = {
            'xlabel': xlabel,
            'ylabel': ylabel
        }
        subplot_default.update(subplot_kw)
        label = kw.pop('label', str(np.round(self.mean_temperature, 2)))
        subplot(3,1,1)
        plt1 = plot(x, y1, axes, subplot_default, fig_kw, label=label, **kw)
        subplot(3,1,2)
        plt2 = plot(x, y2, axes, subplot_default, fig_kw, label=label, **kw)
        subplot(3,1,3)
        plt3 = plot(x, y3, axes, subplot_default, fig_kw, label=label, **kw)
        return plt1, plt2, plt3


class TemperatureScan(Susceptibility, e21.core.Plottable):

    def plot(self, y='real', x='temperature',
             axes=None, subplot_kw={}, fig_kw={}, **kw):
        """A default implementation for a generic fieldscan plot."""
        xa, ya = lookup(self, x), lookup(self, y)
        xlabel = 'T ({0})'.format(
            xa.dimensionality) if x is 'temperature' else ''
        ylabel = r'${\mu}$V' if y is 'real' else ''
        subplot_default = {
            'xlabel': xlabel,
            'ylabel': ylabel
        }
        subplot_default.update(subplot_kw)
        label = kw.pop('label', '{0}'.format(self.field[0]))
        return super(TemperatureScan, self).plot(
            y, x, axes, subplot_default, fig_kw, label=label, **kw)



def create(data, params):
    """Takes data and param and creates a FieldScan or TemperatureScan."""
    mode = params['info']['command']['mode']
    if mode == 'BSWEEP' or mode == 'BSTEP':
        return FieldScan(data, params)
    elif mode == 'TSWEEP' or mode == 'TSTEP':
        return TemperatureScan(data, params)
    else:
        return Susceptibility(data, params)


class Experiment(e21.core.Experiment):

    @property
    def temperature_scans(self):
        return [
            x for x in sorted(
                self._measurements,
                key=lambda x: x.field[0]) if isinstance(
                x,
                TemperatureScan)]

    @property
    def field_scans(self):
        return [
            x for x in sorted(
                self._measurements,
                key=lambda x: x.mean_temperature) if isinstance(
                x,
                FieldScan)]          

    
