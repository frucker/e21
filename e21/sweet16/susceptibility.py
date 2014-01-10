# -*- coding: utf-8 -*-
#
# e21, (c) 2013, see AUTHORS. Licensed under the GNU GPL.
import e21.core

#-Susceptibility-classes-------------------------------------------------------
class Susceptibility(e21.core.Measurement):
    @property
    def x(self):
        return self.data['LI1_CH1']

    def y(self):
        return self.data['LI1_CH2']

    @property
    def chi(self):
        """Calculates the susceptibility of the X signal."""
        raise NotImplementedError()

    @property
    def field(self):
        return self.data['H']

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


class FieldScan(Susceptibility, Plottable):
    def plot(self, y='chi', x='field', axes=None, subplot_kw={}, fig_kw={}, **kw):
        """A default implementation for a generic fieldscan plot."""
        subplot_default = {
            'xlabel': 'B ({0})'.format(x.dimensionality),
            'ylabel': r'$\mathrm{{\chi}}$ ({0})'.format(y.dimensionality),
        }
        subplot_default.update(subplot_kw)
        label = kw.pop('label', str(np.mean(self.control_temperature)))
        return super(FieldScan, self).plot(y, x, axes, subplot_default, fig_kw, label, **kw)


class TemperatureScan(Susceptibility, Plottable):
    def plot(self, y='chi', x='temperature', axes=None, subplot_kw={}, fig_kw={}, **kw):
        """A default implementation for a generic fieldscan plot."""
        subplot_default = {
            'xlabel': 'T ({0})'.format(x.dimensionality),
            'ylabel': r'$\mathrm{{\chi}}$ ({0})'.format(y.dimensionality),
        }
        subplot_default.update(subplot_kw)
        label = kw.pop('label', '{0}'.format(self.field[0]))
        return super(TemperatureScan, self).plot(y, x, axes, subplot_default, fig_kw, label, **kw)


def create(data, params):
    """Takes data and param and creates a FieldScan or TemperatureScan."""
    mode = params['info']['command'].mode
    if mode == 'BSWEEP' or mode == 'BSTEP':
        return FieldScan(data, params)
    elif mode == 'TSWEEP' or mode == 'TSTEP':
        return TemperatureScan(data, params)
    else:
        return Susceptibility(data, params)
