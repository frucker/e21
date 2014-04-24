# -*- coding: utf-8 -*-
#
# e21, (c) 2013, see AUTHORS. Licensed under the GNU GPL.
import e21.core
from e21.core import lookup
import numpy as np
from e21.sweet16 import Loader
import quantities as pq


def _correct_field(field, correction):
	return field - correction * np.abs(field)

#-Susceptibility-classes-------------------------------------------------------
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
    def field_corrected(self):
        """ Correction for Remanent Field offset during Magnetic Field Sweeps."""
        if float(self.params['info']['command']['target_field_rate']) == 0.05:
            return _correct_field(self.field, 0.02)
        else:
            return self.field

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

    @property
    def mean_temperature(self):
        return np.mean(self.temperature)

class FieldScan(Susceptibility, e21.core.Plottable):
    def plot(self, y='real', x='field', axes=None, subplot_kw={}, fig_kw={}, **kw):
        """A default implementation for a generic fieldscan plot."""
        xa, ya = lookup(self, x), lookup(self, y)
        xlabel = 'B ({0})'.format(xa.dimensionality) if x is 'field' else ''
        ylabel = 'U(V)' if y is 'real' else ''
        subplot_default = {
            'xlabel': xlabel,
            'ylabel': ylabel
        }
        subplot_default.update(subplot_kw)
        label = kw.pop('label', str(np.round(np.mean(self.temperature),2)))
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

class AngleScan(Susceptibility, e21.core.Plottable):
    def plot(self, y='real', x='angle', axes=None, subplot_kw={}, fig_kw={}, **kw):
        """A default implementation for a generic fieldscan plot."""
        xa, ya = lookup(self, x), lookup(self, y)
        xlabel = 'Angle (Deg)' if x is 'angle' else ''
        ylabel = r'${\chi}$' if y is 'real' else ''
        subplot_default = {
            'xlabel': xlabel,
            'ylabel': ylabel
        }
        subplot_default.update(subplot_kw)
        label = kw.pop('label', '{0}'.format(self.field[0]))
        return super(AngleScan, self).plot(y, x, axes, subplot_default, fig_kw, label=label, **kw)



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
        return [x for x in sorted(self._measurements, key=lambda x: x.field[0]) if isinstance(x, TemperatureScan)]

    @property
    def field_scans(self):
        return [x for x in sorted(self._measurements, key=lambda x: x.mean_temperature) if isinstance(x, FieldScan)]

    def add_measurements(self, filelist):
        """ Adds measurements to Experiment
            
            filelist: list of measurement files
        """
        s16 = Loader(mode='susceptibility')
        
        meas_len = len(self._measurements)
        testfile = s16(filelist[0])

        
        if not 'dropping_resistance' in testfile.params['general'].keys():
            dropping_resistance = float(raw_input('Dropping Resistance [kOhms]: '))   
        if not 'amplification' in testfile.params['general'].keys():
            amplification = float(raw_input('Amplification: '))  
        
        for n, files in enumerate(filelist):
            index = n + meas_len
            
            self._measurements[index] = s16(files)
            self._measurements[index].params['general']['filename'] = files
            if not 'amplification' in self._measurements[index].params['general'].keys():
               self._measurements[index].params['general']['amplification'] = amplification
            if 'dropping_resistance' in self._measurements[index].params['general'].keys():
               self._measurements[index].params['general']['dropping_resistance'] = float(self._measurements[index].params['general']['dropping_resistance'].strip('kOhm'))
            else:
                self._measurements[index].params['general']['dropping_resistance'] = dropping_resistance
 
    def MakeOverview(self, **kw):                     
        """
        Overview over all measurements including every type of measurement (so far B/T Sweeps)
        KW options:

                files: if true append column indicating paths of measurement files

        """
        html_table = '<h1> %s - Measurement Overview </h1>' % str(self[0].params['info']['sample']) #Create Table Header 
        html_table += '<strong> Sample: </strong> {} <br>'.format(self[0].params['info']['sample']) 
        
        html_table += '<strong> Measurement Option: </strong> Susceptibility <br> <br> '
        html_table += '<table class="table table-striped">'
        html_table += ('<tr> <th> Number </th> <th> Sweep Type </th> <th> Start </th> <th> Stop </th> <th> T / B </th> <th> B_init </th> <th> T_init </th><th> Ampl </th><th> dropping res </th><th> sweep [T,K/min] </th><th> filepath </th><th> reserve </th></tr>')
    
        for num in range(len(self)):
                       
           
            if (type(self[num]) == e21.sweet16.susceptibility.TemperatureScan):
                T_init = float(np.round(self[num].data['sample_temp_1'][0],6))
                T_final = float(np.round(self[num].data['sample_temp_1'][-1],6))
                B_field = float(np.round(np.median(self[num].data['B_field']),3))
                sigma = float(np.round(np.std(self[num].data['B_field']),3))
                html_table += ('<tr> <td>{}</td> <td> Tsweep </td>'
                                    '<td> {} K </td>' 
                                    '<td> {} K </td>' 
                                    '<td> {}+-{} T </td>'
                                    '<td> {} T </td>'
                                    '<td> {} K </td>'
                                    '<td> {} </td>'
                                    '<td> {} k&Omega; </td>'
                                    '<td> {} </td>'
                                    '<td> {} </td>'
                                    '<td> {} </td>'
                                    '</tr>').format(num,
                                    T_init,T_final,B_field,sigma, 
                                    self[num].params['info']['command']['init_field'],
                                    self[num].params['info']['command']['init_temperature'],
                                    self[num].params['general']['amplification'],
                                    self[num].params['general']['dropping_resistance'],
                                    self[num].params['info']['command']['target_temperature_rate'],
                                    self[num].params['general']['filename'],
                                    self[num].params['lock_in_1']['dyn_reserve'])

            elif (type(self[num]) == e21.sweet16.susceptibility.FieldScan):   
                B_init = float(np.round(self[num].data['B_field'][0],6))
                B_final = float(np.round(self[num].data['B_field'][-1],6))
                Temp = float(np.round(np.median(self[num].data['sample_temp_1']),3))               
                sigma = float(np.round(np.std(self[num].data['sample_temp_1']),3))
                html_table += ('<tr> <td>{}</td> <td> Field </td>'
                                    '<td> {} </td>' 
                                    '<td> {} </td>' 
                                    '<td> {}+-{} K </td>'
                                    '<td> {} T </td>'
                                    '<td> {} K </td>'
                                    '<td> {} </td>'
                                    '<td> {} k&Omega; </td>'
                                    '<td> {} </td>'
                                    '<td> {} </td>'
                                    '<td> {} </td>'
                                    '</tr>').format(num,
                                    B_init,B_final,Temp,sigma, 
                                    self[num].params['info']['command']['init_field'],
                                    self[num].params['info']['command']['init_temperature'],
                                    self[num].params['general']['amplification'],
                                    self[num].params['general']['dropping_resistance'],
                                    self[num].params['info']['command']['target_field_rate'],
                                    self[num].params['general']['filename'],
                                    self[num].params['lock_in_1']['dyn_reserve'])
                    
        html_table += '</table>'
        return html_table



