# -*- coding: utf-8 -*-
#
# e21, (c) 2013, see AUTHORS. Licensed under the GNU GPL.
import e21.core
from e21.core import lookup
from e21.sweet16 import Loader
import numpy as np
import operator
import quantities as pq

#-Transport-classes-------------------------------------------------------

class Transport(e21.core.Measurement, e21.core.Plottable):
    
    @property
    def real(self):
        return self.data['LI1_CH1']
    
    @property
    def imag(self):
        return self.data['LI1_CH2']
    
    @property
    def real_hall(self):
        return self.data['LI2_CH1']
    
    @property
    def imag_hall(self):
        return self.data['LI2_CH2']

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
    def geometry_factor(self):
        """
        Calculates geometry factor b*d/(I_ac*Ampl*l_xx) 
        in muOhm*cm/V.
        """
        return (1/self.excitation_current*self.sample_width*
                self.sample_thickness/self.contact_distance/
                self.amplification*10e8)

    @property
    def res(self):
        """Calculates the resistance of the X signal."""
        return self.real*self.geometry_factor

    @property
    def hall(self):
        """Calculates the resistance of the X signal."""
        return self.real_hall*self.geometry_factor
    
    @property
    def sample_thickness(self):
        return float(self.params['general']['sample_thickness'])*pq.mm

    @property
    def sample_width(self):
        return float(self.params['general']['sample_width'])*pq.mm

    @property
    def contact_distance(self):
        return float(self.params['general']['contact_distance'])*pq.mm

    @property
    def field(self):
        return self.data['B_field']

    @property
    def amplification(self):
        return float(self.params['general']['amplification'])

    @property
    def dropping_resistance(self):
        return float(self.params['general']['dropping_resistance'])

    @property
    def excitation_amplitude(self):
        return float(self.params['lock_in_1']['osc_amplitude'].strip('V'))*pq.V

    @property
    def excitation_current(self):
        return self.excitation_amplitude/self.dropping_resistance

    @property
    def target_temperature(self):
        # NOTE: The sweet 16 uses control_temp as the name for the target temp.
        return self.data['control_temp']

    @property
    def temperature(self):
        return self.data['sample_temp_1']

    def temperature_stability(self):
        """
        Calculates the offset and standart deviation of the temperature
        from the target value.
        """
        dT = self.temperature - self.target_temperature
        return np.mean(dT), np.std(dT)

    @property
    def mean_temperature(self):
        return np.mean(self.temperature)

class FieldScan(Transport, e21.core.Plottable):
    def plot(self, y='real', x='field',
            axes=None, subplot_kw={}, fig_kw={}, **kw):
        """A default implementation for a generic fieldscan plot."""
        xa, ya = lookup(self, x), lookup(self, y)
        xlabel = 'B ({0})'.format(xa.dimensionality) if x is 'B_field' else ''
        ylabel = 'V' if y is 'real' else ''
        subplot_default = {
            'xlabel': xlabel,
            'ylabel': ylabel
        }
        subplot_default.update(subplot_kw)
        label = kw.pop('label', str(np.mean(self.target_temperature)))
        return super(FieldScan, self).plot(y, x, axes,
                                          subplot_default,
                                          fig_kw, label=label, **kw)

    def plot_res(self, y='res', x='field', axes=None,
                subplot_kw={}, fig_kw={}, **kw):
        """A default implementation for a generic fieldscan plot."""
        xa, ya = lookup(self, x), lookup(self, y)
        xlabel = 'B ({0})'.format(xa.dimensionality) if x is 'field' else ''
        ylabel = '$ \\rho_{xx}$ ($\mu\Omega\cdot$cm)' if y is 'res' else ''
        subplot_default = {
            'xlabel': xlabel,
            'ylabel': ylabel
        }
        subplot_default.update(subplot_kw)
        label = kw.pop('label', str(np.round(np.mean(self.target_temperature),2)))
        return super(FieldScan, self).plot(y, x, axes,
                                          subplot_default,
                                          fig_kw, label=label, **kw)

    def plot_hall(self, y='hall', x='field', axes=None,
                 subplot_kw={}, fig_kw={}, **kw):
        """A default implementation for a generic fieldscan plot."""
        xa, ya = lookup(self, x), lookup(self, y)
        xlabel = 'B ({0})'.format(xa.dimensionality) if x is 'field' else ''
        ylabel = '$ \\rho_{xy}$ ($\mu\Omega\cdot$cm)' if y is 'hall' else ''
        subplot_default = {
            'xlabel': xlabel,
            'ylabel': ylabel
        }
        subplot_default.update(subplot_kw)
        label = kw.pop('label', str(np.round(np.mean(self.target_temperature),2)))
        return super(FieldScan, self).plot(y, x, axes,
                                          subplot_default,
                                          fig_kw, label=label, **kw)

class TemperatureScan(Transport, e21.core.Plottable):
    def plot(self, y='real', x='temperature', axes=None,
            subplot_kw={}, fig_kw={}, **kw):
        """A default implementation for a generic temperature scan plot."""
        xa, ya = lookup(self, x), lookup(self, y)
        xlabel = 'T ({0})'.format(xa.dimensionality) if x is 'temperature' else ''
        ylabel = 'V' if y is 'real' else ''
        subplot_default = {
            'xlabel': xlabel,
            'ylabel': ylabel
        }
        subplot_default.update(subplot_kw)
        label = kw.pop('label', '{0}'.format(self.field[0]))
        return super(TemperatureScan, self).plot(y, x, axes,
                                                subplot_default,
                                                fig_kw, label=label, **kw)

    def plot_res(self, y='res', x='temperature', axes=None,
                subplot_kw={}, fig_kw={}, **kw):
        """A default implementation for a generic fieldscan plot."""
        xa, ya = lookup(self, x), lookup(self, y)
        xlabel = 'T ({0})'.format(xa.dimensionality) if x is 'temperature' else ''
        ylabel = '$\\rho_{xx}$ ($\mu\Omega\cdot$cm)' if y is 'res' else ''
        subplot_default = {
            'xlabel': xlabel,
            'ylabel': ylabel
        }
        subplot_default.update(subplot_kw)
        label = kw.pop('label', str(np.round(np.mean(self.field),2)))
        return super(TemperatureScan, self).plot(y, x, axes,
                                                subplot_default,
                                                fig_kw, label=label, **kw)

    def plot_hall(self, y='hall', x='temperature', axes=None,
                 subplot_kw={}, fig_kw={}, **kw):
        """A default implementation for a generic fieldscan plot."""
        xa, ya = lookup(self, x), lookup(self, y)
        xlabel = 'T ({0})'.format(xa.dimensionality) if x is 'temperature' else ''
        ylabel = '$\\rho_{xy}$ ($\mu\Omega\cdot$cm)' if y is 'hall' else ''
        subplot_default = {
            'xlabel': xlabel,
            'ylabel': ylabel
        }
        subplot_default.update(subplot_kw)
        label = kw.pop('label', str(np.round(np.mean(self.field),2)))
        return super(TemperatureScan, self).plot(y, x, axes,
                                                subplot_default,
                                                fig_kw, label=label, **kw)
    
def create(data, params):
    """Takes data and param and creates a FieldScan or TemperatureScan."""
    mode = params['info']['command']['mode']
    if mode == 'BSWEEP' or mode == 'BSTEP':
        return FieldScan(data, params)
    elif mode == 'TSWEEP' or mode == 'TSTEP':
        return TemperatureScan(data, params)
    else:   
        return Transport(data, params)

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
        s16 = Loader(mode='transport')
        
        meas_len = len(self._measurements)
        
        testfile = s16(filelist[0])

        if not 'contact_distance' in testfile.params['general'].keys():
            contact_distance = float(raw_input('Contact distance [mm]: '))
        if not 'sample_width' in testfile.params['general'].keys():
            sample_width = float(raw_input('Sample width [mm]: '))
        if not 'sample_thickness' in testfile.params['general'].keys():
            sample_thickness = float(raw_input('Sample thickness [mm]: '))
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
            if not 'contact_distance' in self._measurements[index].params['general'].keys():
                self._measurements[index].params['general']['contact_distance'] = contact_distance
            if not 'sample_width' in self._measurements[index].params['general'].keys():
                self._measurements[index].params['general']['sample_width'] = sample_width
            if not 'sample_thickness' in self._measurements[index].params['general'].keys():
                self._measurements[index].params['general']['sample_thickness'] = sample_thickness

    def MakeOverview(self, **kw):    
        """
        Overview over all measurements including every type of measurement (so far B/T Sweeps)
        KW options:

                files: if true append column indicating paths of measurement files

        """
        
        html_table = '<table border = "1">'
        #html_table += '<tr><h1> %s - Measurement Overview </h1></tr>'%str(self[0].params['info']['sample']) #Create Table Header 
        #html_table += '<strong> Sample: </strong> {} <br>'.format(self[0].params['info']['sample']) 
        html_table += '<strong> Contact Distance: </strong> {}<br>'.format(self[0].contact_distance)
        html_table += '<strong> Measurement Option: </strong> Transport <br> <br> '
        html_table += '<p style="color:blue;margin-left:20px;">This is a paragraph.</p>  '
        html_table += '<table>'
        html_table += ('<tr> <th> Number </th> <th> Sweep Type </th> <th> Start </th> <th> Stop </th> <th> T / B </th> <th> B_init </th> <th> T_init </th>'
                      '<th> Amplification </th><th> dropping_resistance </th><th> filepath </th><th> reserve </th></tr>')
    
        for num in range(len(self)):
                       
           
            if (type(self[num]) == e21.sweet16.transport.TemperatureScan):
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
                                    '</tr>').format(num,
                                    T_init,T_final,B_field,sigma, 
                                    self[num].params['info']['command']['init_field'],
                                    self[num].params['info']['command']['init_temperature'],
                                    self[num].params['general']['amplification'],
                                    self[num].params['general']['dropping_resistance'],
                                    self[num].params['general']['filename'],
                                    self[num].params['lock_in_1']['dyn_reserve'])

            elif (type(self[num]) == e21.sweet16.transport.FieldScan):   
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
                                    '</tr>').format(num,
                                    B_init,B_final,Temp,sigma, 
                                    self[num].params['info']['command']['init_field'],
                                    self[num].params['info']['command']['init_temperature'],
                                    self[num].params['general']['amplification'],
                                    self[num].params['general']['dropping_resistance'],
                                    self[num].params['lock_in_1']['dyn_reserve'])
                    
        html_table += '</table>'
        return html_table






