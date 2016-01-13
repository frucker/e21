import e21.core
from e21.core import lookup
import numpy as np
from e21.Nicos.Dryo import Loader
import quantities as pq


class Dryo(e21.core.Measurement, e21.core.Plottable):

# Mandatory Properties for new Instrument cores

    @property
    def field(self):
        try:        
            return self.data['sample-field']
        except KeyError: 
            return float(self.params['devices']['sample_magnet_value'].strip('T\n'))
    @property
    def mean_field(self):
        try:
            return np.mean(self.field)
        except KeyError:
            return self.params['devices']['sample_magnet_value'].strip('T\n')

    @property
    def mean_temperature(self):
        return np.median(self.temperature)

    @property
    def cernox_temperature(self):
        try:
            return self.data['temperature-lts_cernox_thermometer']
        except KeyError:
            return float(self.params['devices']['lts_cernox_thermometer_value'].strip('K\n'))
    @property
    def rox_temperature(self):
        try:
            return self.data['temperature-lts_rox_thermometer']
        except KeyError:
            return float(self.params['devices']['lts_rox_thermometer_value'].strip('K\n'))

    @property
    def sweep_rate(self):
        try:
            info = self.params['general']['info'].split('(')
            if self.sweep_type == 'Tscan':
                return info[1].split(',')[4].strip(']').strip()
            elif self.sweep_type == 'Bscan':
                return info[1].split(',')[4].strip(']').strip()
            else:
                return 0.0

        except:
            print '{}'.format(self.filename)
        
    @property
    def target_temperature_rate(self):
        return self.sweep_rate

    @property
    def target_field_rate(self):
        return self.sweep_rate

    @property
    def temperature(self):
        if self.sweep_type == 'Tscan':
            info = self.params['general']['info'].split('(')
            sweep = info[1].split(',')[0]
            if sweep == 'adr_servo':
                return self.data['temperature-adr_servo']
            else:
                return self.cernox_temperature           
        else:
            if np.median(self.cernox_temperature) < 1.5:
                return self.rox_temperature
            else:
                return self.cernox_temperature

    @property
    def sample(self):
        # No sample information in header
        return 'no Info'

    @property
    def init_field(self):
        return float(self.params['devices']['sample_magnet_value'].split('(')[0].strip('T\n'))

    @property
    def init_temperature(self):
        if float(self.params['devices']['lts_rox_thermometer_value'].strip('K\n')) > 4:
            return float(self.params['devices']['lts_rox_thermometer_value'].strip('K\n'))
        else:
            return float(self.params['devices']['lts_cernox_thermometer_value'].strip('K\n'))

 
    @property
    def sweep_type(self):
        try:
            info = self.params['general']['info'].split('(')
            typ = info[0].strip()
            if typ == 'sweep':
                sweep = info[1].split(',')[0]
                if sweep in['htc_servo','adr_servo']:
                    return 'Tscan'
                elif sweep == 'adr_magnet':
                    return 'ADR Magnet'
                elif sweep == 'sample_magnet':
                    return 'Bscan'
            else:   
                return info[0].strip()
        except KeyError:
            return 'None'
 
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
