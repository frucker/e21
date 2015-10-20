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

    

class Susceptibility(e21.core.Measurement, e21.core.Plottable):

    @property
    def real(self):
        return self.data['LI1_CH2']

    @property
    def imag(self):
        return self.data['LI1_CH1']

    @property
    def angle(self):
        try:
            return self.data['angle']
        except KeyError:
            return 0. * pq.deg

    @property
    def NV(self):
        return self.params['info']['command']['needle_valve_percentage']

    @property
    def time(self):
        return self.data['datetime']

    @property
    def heater(self):
        return self.data['heater']

    @property
    def needle_valve(self):
        return self.data['needle_valve']

    @property
    def chi(self):
        """Calculates the susceptibility of the X signal."""
        raise NotImplementedError()

    @property
    def field(self):
        return self.data['B_field']

    @property
    def mean_field(self):
        return np.mean(self.field)

    @property
    def mean_current(self):
        return np.round(np.mean(self.data['current']), 2)

    @property
    def field_corrected(self):
        """
        Correction for Remanent Field offset during Magnetic Field Sweeps.

        """
        if float(self.params['info']['command']['target_field_rate']) == 0.05:
            return _correct_field(self.field, 0.03)
        elif float(self.params['info']['command']['target_field_rate']) == 0.02:
            return _correct_field(self.field, 0.01)
        else:
            return self.field

    @property
    def target_temperature(self):
        # NOTE: The sweet 16 uses control_temp as the name for the target temp.
        return self.data['control_temp']

    @property
    def temperature(self):
        try:        
            return self.data['sample_temp_1']
        except KeyError:
            return self.data['MC_LS372_2']

    def temperature_stability(self):
        """Calculates the offset and standart deviation of the temperature
        from the target value.
        """
        dT = self.temperature - self.target_temperature
        return np.mean(dT), np.std(dT)

    @property
    def mean_temperature(self):
        return np.mean(self.temperature)

    @property
    def init_temp(self):
        return self.params['info']['command']['init_temperature']


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
        Test = s16(files) 
        try:
            if Test.data['datetime'][0]:
                fl.append(files)
        except IndexError:
            print 'file empty: {}'.format(files)
    return fl

def create(data, params):
    """Takes data and param and creates a FieldScan or TemperatureScan."""
    mode = params['info']['command']['mode']
    if mode == 'BSWEEP' or mode == 'BSTEP':
        return FieldScan(data, params)
    elif mode == 'TSWEEP' or mode == 'TSTEP':
        return TemperatureScan(data, params)
    elif mode == 'ASWEEP':
        return AngleScan(data, params)
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

    def add_measurements(self, filelist):
        """ Adds measurements to Experiment

            filelist: list of measurement files
        """
        s16 = Loader(mode='susceptibility')

        meas_len = len(self._measurements)
        filelist = check_empty_files(filelist)
        testfile = s16(filelist[0])

        if 'dropping_resistance' not in testfile.params['general'].keys():
            dropping_resistance = float(
                raw_input('Dropping Resistance [kOhms]: '))
        if 'amplification' not in testfile.params['general'].keys():
            amplification = float(raw_input('Amplification: '))

        for n, files in enumerate(filelist):
            index = n + meas_len

            self._measurements[index] = s16(files)
            self._measurements[index].params['general']['filename'] = files
            if 'amplification' not in self._measurements[
                    index].params['general'].keys():
                self._measurements[index].params['general'][
                    'amplification'] = amplification
            if 'current' not in self._measurements[index].data.keys():
                self._measurements[index].data['current'] = 0
            if 'dropping_resistance' in self._measurements[
                    index].params['general'].keys():
                self._measurements[index].params['general']['dropping_resistance'] = float(
                    self._measurements[index].params['general']['dropping_resistance'].strip('kOhm'))
            else:
                self._measurements[index].params['general'][
                    'dropping_resistance'] = dropping_resistance

