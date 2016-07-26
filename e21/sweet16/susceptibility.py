# -*- coding: utf-8 -*-
#
# e21, (c) 2013, see AUTHORS. Licensed under the GNU GPL.
import e21.core
from e21.sweet16.core import Sweet16
from e21.core import lookup
import numpy as np
from e21.sweet16 import Loader
import quantities as pq
from e21.utility import ProgressBar


def _correct_field(field, correction):
    return field - correction * pq.T

# Susceptibility-classes ------------------------------------------------------

class Susceptibility(Sweet16):

    @property
    def experiment_type(self):
        return 'susceptibility'

    @property
    def amplification(self):
        return self.params['general']['amplification']

    @property
    def dropping_resistance(self):
        return self.params['general']['dropping_resistance']

    @property
    def real(self):
        return self.data['LI1_CH2']

    @property
    def imag(self):
        return self.data['LI1_CH1']

    @property
    def chi(self):
        """Calculates the susceptibility of the X signal."""
        raise NotImplementedError()

        
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
    """ 
	
		ACHTUNG! Dauert viel zu lange...
		Checks for empty files in list of files.
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
	
def check_empty_files(filelist, mode = 'susceptibility'):
    """ 
	
		ACHTUNG! Dauert viel zu lange...
		Checks for empty files in list of files.
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
        #filelist = check_empty_files(filelist)
        testfile = s16(filelist[0])
        p = ProgressBar(len(filelist))

        if 'dropping_resistance' not in testfile.params['general'].keys():
            dropping_resistance = float(
                raw_input('Dropping Resistance [kOhms]: '))
        if 'amplification' not in testfile.params['general'].keys():
            amplification = float(raw_input('Amplification: '))
        emptyfiles = '\nEmpty files:\n'
        i = 0
        for n, file in enumerate(filelist):
			
			# i counting only non-empty files
			
			# Count import of all files (including empty ones)
			p.animate(n+1)
			# Check if imported file is empty, cheap trick..
			try:
				Test = s16(file)
				if Test.data['datetime'][0]:
					index = i + meas_len
					self._measurements[index] = Test
					self._measurements[index].params['general']['filename'] = file
					if 'amplification' not in self._measurements[
							index].params['general'].keys():
						self._measurements[index].params['general'][
							'amplification'] = amplification
					if 'current' not in self._measurements[index].data.keys():
						if 'k2440_current' in self._measurements[index].data.keys():
							self._measurements[index].data['current'] = self._measurements[index].data['k2440_current']
						else:
							self._measurements[index].data['current'] = 0 
					if 'dropping_resistance' in self._measurements[
							index].params['general'].keys():
						self._measurements[index].params['general']['dropping_resistance'] = float(
							self._measurements[index].params['general']['dropping_resistance'].strip('kOhm'))
					else:
						self._measurements[index].params['general'][
							'dropping_resistance'] = dropping_resistance
					# only count if file was not empty	
					i = i+1
			except IndexError:
				emptyfiles = emptyfiles + file + '\n'	
        print emptyfiles

def measurement_details(Exp):
    '''
    returns a function g(i) which can be used with ipython.widges.interact

    input parameters:

        Exp: e21 Experiment class dicitonary

    output:
    
        g(i): function that displays a HTML table containing comprehensive
              measurement information on Measurement Exp[i]

    ''' 


    def g(i = 0):
        s = '<h2> Measurement Details for Measurement Number {}</h2>'.format(i)
        s += '<h3> General information: </h3>'
        s += '<table class="table table-hover">'
        s += ('<tr> <td><strong> Date </strong></td>'
             ' <td> {} </td></tr>'.format(Exp[i].params['info']['date']))
        s += ('<tr> <td><strong> Sample </strong></td>'
             '<td> {} </td></tr>'.format(Exp[i].params['info']['sample']))
        s += ('<tr> <td><strong> Time </strong></td>'
             '<td> {} </td></tr>'.format(Exp[i].params['info']['time']))
        s += ('<tr> <td><strong> User </strong></td>'
             '<td> {} </td></tr>'.format(Exp[i].params['info']['user']))
        s += ('<tr> <td><strong> Mode </strong></td>'
             '<td> {} </td></tr>'.format(Exp[i].params['mode']))
        s += ('<tr> <td><strong> Amplification </strong></td> <td> {} </td>'
             '</tr>'.format(Exp[i].params['general']['amplification']))
        s += ('<tr> <td><strong> Dropping Resistance </strong></td> <td>'
             ' {} k$\Omega$ </td></tr>'.format(
                 Exp[i].params['general']['dropping_resistance']))
        s += ('<tr> <td><strong> Filename </strong></td> <td> {} </td>'
             '</tr>'.format(Exp[i].params['general']['filename']))
        try: 
            if Exp[i].params['general']['contact_distance']:
                s += ('<tr> <td><strong> Contact Distance </strong></td>'
                     '<td> {} </td></tr>'.format(Exp[i].params['general']['contact_distance']))
                s += ('<tr> <td><strong> Sample Thickness </strong></td>'
                     '<td> {} </td></tr>'.format(Exp[i].params['general']['sample_thickness']))
                s += ('<tr> <td><strong> Sample Width </strong></td> <td> {} </td>'
                     '</tr>'.format(Exp[i].params['general']['sample_width']))
        except:
            pass
        s += '</table>'     
        s += '<h3> Lockin Information: </h3>'
        s += '<table class="table table-hover">' 
        try:
            keys = Exp[i].params['lock_in_1'].items()
            keys_sorted = sorted(keys, key=itemgetter(0))       
            for j in range(len(keys_sorted)):
                s += '<tr> <td><strong> {} </strong></td> <td> {} </td>'.format(keys_sorted[j][0],keys_sorted[j][1])
                if 'lock_in_2' in Exp[i].params.keys():
                    s +='<td> {} </td>'.format(Exp[i].params['lock_in_2'][keys_sorted[j][0]])
                if 'lock_in_3' in Exp[i].params.keys():
                    s +='<td> {} </td></tr> '.format(Exp[i].params['lock_in_3'][keys_sorted[j][0]])
                else:
                    s += '</tr>' 
        except KeyError:
            pass
            
        s += '</table>'
        s += '<h3> Temperature and Field Info: </h3>'
        s += '<table class="table table-hover">' 
        s += ('<tr> <td><strong> Mean Temperature </strong></td>'
             ' <td> {} </td></tr>'.format(Exp[i].mean_temperature))
        s += ('<tr> <td><strong> Mean Field </strong></td>'
             ' <td> {} </td></tr>'.format(Exp[i].mean_field))
        s += '</table>'
        s += '<h3> Command File: </h3>'
        s += '<table class="table table-hover">'
        keys = Exp[i].params['info']['command'].items()
        keys_sorted = sorted(keys, key=itemgetter(0))
        for j in range(len(keys_sorted)):
            s += '<tr> <td><strong> {} </strong></td> <td> {} </td></tr>'.format(keys_sorted[j][0],keys_sorted[j][1])      
        s = display(HTML(s))
    return g

