import e21.core
from e21.core import lookup
import numpy as np
from e21.sweet16 import Loader
import quantities as pq

class Sweet16(e21.core.Measurement, e21.core.Plottable):

# Mandatory Properties for new Instrument cores

    @property
    def sample(self):
        return self.params['info']['sample']
    @property
    def init_field(self):
        return float(self.params['info']['command']['init_field'].strip(' T'))*pq.T

    @property
    def field(self):
        return self.data['B_field']

    @property
    def mean_field(self):
        return np.mean(self.field)

    @property
    def init_temperature(self):
        return float(self.params['info']['command']['init_temperature'].strip(' K'))*pq.K

    @property
    def target_temperature_rate(self):
        return self.params['info']['command']['target_temperature_rate']

    @property
    def target_field_rate(self):
        return self.params['info']['command']['target_field_rate']

    @property
    def sweep_type(self):
        return self.params['info']['command']['mode']

    @property
    def temperature(self):
        try:        
            return self.data['LS340_temp']
        except KeyError:
            try: 
                return self.data['sample_temp_1']
            except KeyError:     
                try:
                    return self.data['MC_LS372_1']
                except KeyError:
                    try:
                        return self.data['MC_temp']
                    except KeyError:
                        return 0
      

    @property
    def mean_angle(self):
        try:
            return round(median(self.data['angle']),2)
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
        return self.params['general']['filename'].split('/')[-1]

# Optional Properties

    @property
    def amplitude(self):
        return self.params['info']['command']['K6221 current']

    @property
    def frequency(self):
        return self.params['info']['command']['k6221 frequency']

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
    def mean_current(self):
        return np.round(np.mean(self.data['current']), 3)
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

    @property
    def NV(self):
        return self.params['info']['command']['needle_valve_percentage']

