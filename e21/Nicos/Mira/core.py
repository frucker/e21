import e21.core
from e21.core import lookup
import numpy as np
from e21.Nicos.Mira import Loader
import quantities as pq


class Mira(e21.core.Measurement, e21.core.Plottable):

# Mandatory Properties for new Instrument cores

    @property
    def field(self):
        try:        
            return self.data['dct2'] 
        except KeyError: 
            try:
                return [float(self.params['devices']['dct2_value'].strip('V\n'))*pq.V]*len(self)
            except KeyError:
                return [0*pq.T]*len(self)
    
    @property
    def mean_field(self):
        return np.mean(self.field)


    @property
    def mean_temperature(self):
        return np.median(self.temperature)

    @property
    def sweep_rate(self):
        return 0.0 # Not in file
        
    @property
    def target_temperature_rate(self):
        return self.sweep_rate

    @property
    def target_field_rate(self):
        return self.sweep_rate

    @property
    def temperature(self):
        try:
            return self.data['Ts']
        except:
            return [float(self.params['devices']
                   ['Ts_value'].strip('K\n'))]*len(self)


    @property
    def sample(self):
        return self.params['Sample and alignment']['Sample_samplename']

    @property
    def init_field(self):
        return self.field[0]

    @property
    def init_temperature(self):
        return self.temperature[0]
 
    @property
    def sweep_type(self):
        try:
            info = self.params['general']['info'].split('(')
            typ = info[0].strip()
            if typ == 'sweep':
                sweep = info[1].split(',')[0]
                if sweep in['dct2']:
                    return 'Bsweep'
            elif typ == 'qscan':
                return 'Qscan'         
            else:   
                return info[0].strip()
        except KeyError:
            return 'None'

    def sweep_rate(self):
        try:
            info = self.params['general']['info'].split('(')
            if self.sweep_type == 'Qscan':
                return 'to be done'
            elif self.sweep_type == 'Bsweep':
                return 'to be done'
            else:
                return 0.0
        except:
            return 0
 
    @property
    def mean_angle(self):
        try:
            return round(np.median(self.angle),2)
        except KeyError:
            return 0. * pq.deg

    @property
    def angle(self):
        try:
            return self.data['angle']
        except KeyError:
            # If not present in Data, return array of length of measurement
            # with Angle 0.
            return [0*pq.deg]*len(self)

    @property
    def filename(self):
        return self.params['general']['filename']
