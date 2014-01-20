# -*- coding: utf-8 -*-
#
# e21, (c) 2013, see AUTHORS. Licensed under the GNU GPL.
import e21.core
from e21.core import lookup

#-Susceptibility-classes-------------------------------------------------------
class Susceptibility(e21.core.Measurement):
    @property
    def real(self):
        return self.data['LI1_CH2']

    def imag(self):
        return self.data['LI1_CH1']

    @property
    def chi(self):
        """Calculates the susceptibility of the X signal."""
        raise NotImplementedError()

    @property
    def field(self):
        return self.data['B_field']

    @property
    def target_temperature(self):
        # NOTE: The sweet 16 uses control_temp as the name for the target temp.
        return self.data['control_temp']

    @property
    def temperature(self):
        return self.data['sample_temp_1']

    def temperature_stability(self):
        """Calculates the offset and standart deviation of the temperature
        from the target value.
        """
        dT = self.temperature - self.target_temperature
        return np.mean(dT), np.std(dT)


class FieldScan(Susceptibility, e21.core.Plottable):
    def plot(self, y='real', x='field', axes=None, subplot_kw={}, fig_kw={}, **kw):
        """A default implementation for a generic fieldscan plot."""
        x, y = lookup(self, x), lookup(self, y)
        subplot_default = {
            'xlabel': 'B ({0})'.format(x.dimensionality),
            'ylabel': r'$\mathrm{{\chi}}$ ({0})'.format(y.dimensionality),
        }
        subplot_default.update(subplot_kw)
        label = kw.pop('label', str(np.mean(self.control_temperature)))
        return super(FieldScan, self).plot(y, x, axes, subplot_default, fig_kw, label=label, **kw)


class TemperatureScan(Susceptibility, e21.core.Plottable):
    def plot(self, y='real', x='temperature', axes=None, subplot_kw={}, fig_kw={}, **kw):
        """A default implementation for a generic fieldscan plot."""
        xa, ya = lookup(self, x), lookup(self, y)
        xlabel = 'T ({0})'.format(xa.dimensionality) if x is 'temperature' else ''
        ylabel = r'${\chi}$' if y is 'real' else ''
        subplot_default = {
            'xlabel': xlabel,
            'ylabel': ylabel
        }
        subplot_default.update(subplot_kw)
        label = kw.pop('label', '{0}'.format(self.field[0]))
        return super(TemperatureScan, self).plot(y, x, axes, subplot_default, fig_kw, label=label, **kw)


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
    def temperature_scan(self):
        return [x for x in self._measurements if isinstance(x, TemperatureScan)]

    @property
    def field_scan(self):
        return [x for x in self._measurements if isinstance(x, FieldScan)]

