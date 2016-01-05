# -*- coding: utf-8 -*-
#
# e21, (c) 2013, see AUTHORS. Licensed under the GNU GPL.
import e21.core
from e21.core import lookup
import numpy as np
from e21.Nicos.loader import Loader
import quantities as pq


def _correct_field(field, correction):
    return field - correction * pq.T

# Susceptibility-classes ------------------------------------------------------

    

class Susceptibility(e21.core.Measurement, e21.core.Plottable):

    @property
    def real(self):
        return self.data['lia1_Y']

    @property
    def imag(self):
        return self.data['lia1_X']


    @property
    def chi(self):
        """Calculates the susceptibility of the X signal."""
        raise NotImplementedError()

    @property
    def field(self):
        return self.data['sample-field']

    @property
    def mean_field(self):
        try:
            return np.mean(self.field)
        except KeyError:
            return self.params['Device positions and sample environment state']['sample_magnet_value'].strip('T\n')

    @property
    def mean_temperature(self):
        try:
            return float(self.params['Device positions and sample environment state']['adr_control_thermometer_value'].strip('K\n'))
        except:
            pass
    @property
    def cernox_temperature(self):
        try:
            return self.data['temperature-lts_cernox_thermometer']
        except KeyError:
            return float(self.params['Device positions and sample environment state']['lts_cernox_thermometer_value'].strip('K\n'))
    @property
    def rox_temperature(self):
        try:
            return self.data['temperature-lts_rox_thermometer']
        except KeyError:
            return float(self.params['Device positions and sample environment state']['lts_rox_thermometer_value'].strip('K\n'))
    @property
    def sweep_rate(self):
        try:
            info = self.params['NICOS data file']['info'].split('(')
            rate= float(self.params['Device positions and sample environment state']['sample_magnet_rate'].strip('T/min\n'))
            if  info[0].strip() == 'sweep':
                if info[1].split(',')[0] == 'sample_magnet':
                    rate = info[1].split(',')[4].strip(']').strip()
            return rate
        except KeyError:
            return float(self.params['Device positions and sample environment state']['sample_magnet_rate'].strip('T/min\n'))



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

def check_empty_files(filelist, mode = 'susceptibility'):
    """ Checks for empty files in list of files.
        Empty files are files containing only header lines.
        input: filelist
        output: filelist, list of non-empty files
    """ 

    s16 = Loader(mode=mode)
    fl = []
    for n, files in enumerate(filelist):
        fl.append(files)
    return fl

def create(data, params):
    """Takes data and param and creates a FieldScan or TemperatureScan."""
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

    def add_measurements(self, filelist):
        """ Adds measurements to Experiment

            filelist: list of measurement files
        """
        s16 = Loader(mode='susceptibility')

        meas_len = len(self._measurements)
        filelist = check_empty_files(filelist)
        testfile = s16(filelist[0])

        
        for n, files in enumerate(filelist):
            index = n + meas_len

            self._measurements[index] = s16(files)
            self._measurements[index].params['filename'] = files
 

