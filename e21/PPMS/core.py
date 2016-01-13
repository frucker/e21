import e21.core
from e21.core import lookup
import numpy as np
from e21.PPMS import Loader
import quantities as pq
import warnings

class PPMS(e21.core.Measurement, e21.core.Plottable):

# Mandatory Properties for new Instrument cores
    @property
    def sample(self):
        try:
            warnings.warn('File Folder is used as Sample Name')
            return self.params['info']['filepath'].split('/')[-2]
            
            
        except:
            return ''
            warnings.warn('No sample folder found')

    @property
    def init_field(self):
        return float(self.data['Magnetic Field (Oe)'][0])/10000*pq.T

    @property
    def field(self):
        return [i/10000 for i in self.data['Magnetic Field (Oe)']]*pq.T

    @property
    def mean_field(self):
        return np.mean(self.field)

    @property
    def init_temperature(self):
        return float(self.data['Temperature (K)'][0])*pq.K

    @property
    def target_temperature_rate(self):
        return 'not included'

    @property
    def target_field_rate(self):
        return 'not included'

    @property
    def sweep_type(self):
        if self.filename.startswith('T'):
            return 'Tsweep'
        elif self.filename.startswith('B'):
            return 'Bsweep'
        elif self.filename.startswith('F'):
            return 'Bsweep'
        else:
            return 'Scan'

    @property
    def temperature(self):
        return self.data['Temperature (K)']

    @property      
    def temperature_stability(self):
        """Calculates the offset and standart deviation of the temperature
        from the target value.
        """
        dT = self.temperature - self.target_temperature
        return np.mean(dT), np.std(dT)

    @property
    def mean_angle(self):
        return 0. * pq.deg

    @property
    def angle(self):
        return [0*pq.deg]*len(self)

    @property
    def filename(self):
        return self.params['general']['filename'].split('/')[-1]

    @property
    def time(self):
        return self.data['Time Stamp (sec)']
