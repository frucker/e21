# -*- coding: utf-8 -*-
#
# e21, (c) 2013, see AUTHORS. Licensed under the GNU GPL.
import e21.core
from e21.core import lookup
from e21.PPMS.core import PPMS
import numpy as np
from e21.PPMS import Loader
import quantities as pq
from e21.utility import ProgressBar


def _correct_field(field, correction):
    return field - correction * pq.T

# Susceptibility-classes ------------------------------------------------------


class Susceptibility(PPMS):

    @property
    def experiment_type(self):
        return 'susceptibility'

    @property
    def real(self):
        return self.data["M' (emu)"]

    @property
    def imag(self):
        return self.data["M'' (emu)"]

    @property
    def M(self):
        return self.data['M-DC (emu)']

    @property
    def amplitude(self):
        return np.round(np.median(data['Amplitude (Oe)']),2)*10000*pq.T

    @property
    def frequency(self):
        return np.round(np.median(data['Frequency (Hz)']),2)*pq.Hz

class FieldScan(Susceptibility, e21.core.Plottable):

    def plot(self, y='real', x='field', axes=None,
             subplot_kw={}, fig_kw={}, **kw):
        """A default implementation for a generic fieldscan plot."""
        xa, ya = lookup(self, x), lookup(self, y)
        xlabel = 'B ({0})'.format(xa.dimensionality) if x is 'field' else ''
        ylabel = 'U(V)' if y is 'real' else ''
        subplot_default = {
            'xlabel': xlabel,
            'ylabel': ylabel
        }
        subplot_default.update(subplot_kw)
        label = kw.pop('label', str(np.round(np.mean(self.temperature), 2)))
        return super(FieldScan, self).plot(
            y, x, axes, subplot_default, fig_kw, label=label, **kw)


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


class AngleScan(Susceptibility, e21.core.Plottable):

    def plot(self, y='real', x='angle', axes=None,
             subplot_kw={}, fig_kw={}, **kw):
        """A default implementation for a generic fieldscan plot."""
        xa, ya = lookup(self, x), lookup(self, y)
        xlabel = 'Angle (Deg)' if x is 'angle' else ''
        ylabel = r'${\mu}$V' if y is 'real' else ''
        subplot_default = {
            'xlabel': xlabel,
            'ylabel': ylabel
        }
        subplot_default.update(subplot_kw)
        label = kw.pop('label', '{0}'.format(self.field[0]))
        return super(AngleScan, self).plot(
            y, x, axes, subplot_default, fig_kw, label=label, **kw)


def create(data, params):
    return Susceptibility(data, params)


def sort_experiment(Exp, key='mean_field'):
    sorted_list = sorted(
        Exp._measurements.items(),
        key=lambda x: float(
            x[1].__getattribute__(key)))
    sorted_exp = Experiment()
    for i in range(len(sorted_list)):
        sorted_exp[i] = sorted_list[i][1]
    return sorted_exp

def check_empty_files(filelist):
    """ Checks for empty files in list of files.
        Empty files are files containing only header lines.
        input: filelist
        output: filelist, list of non-empty files
    """ 

    s16 = Loader(mode='susceptibility')
    fl = []
    for n, files in enumerate(filelist):
        Test = s16(files)
        try:
            if Test.data['Temperature (K)'][0]:
                fl.append(files)
        except IndexError:
            print 'file empty: {}'.format(files)
    return fl

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

    def add_measurements(self, filelist):
        """ Adds measurements to Experiment

            filelist: list of measurement files
        """
        p = ProgressBar(len(filelist))
        s16 = Loader(mode='susceptibility')

        meas_len = len(self._measurements)
        filelist = check_empty_files(filelist)
        testfile = s16(filelist[0])
        
        for n,files in enumerate(filelist):
            p.animate(n+1)
            index = n + meas_len
            self._measurements[index] = s16(files)
            self._measurements[index].params['general']['filename'] = files
            n = n+ 1
     
