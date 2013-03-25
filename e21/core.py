# -*- coding: utf-8 -*-
#
# e21, (c) 2013, see AUTHORS. Licensed under the GNU GPL.
import jinja2
import markdown2


env = jinja2.Environment(loader=jinja2.PackageLoader('e21', 'templates'))


def markdown(value):
    return markdown2.markdown(value)

def temperature_stability(value):
    offset, deviation = value
    return u'{0:.3f} {1} \u00B1 {2:.3f} {3}'.format(offset.item(), offset.dimensionality.unicode, deviation.item(), deviation.dimensionality.unicode)

# register markdown filter
env.filters['markdown'] = markdown
env.filters['temp_stability'] = temperature_stability


class Html(object):
    def _repr_html_(self):
        try:
            path = '{0}.{1}.html'.format(type(self).__module__, type(self).__name__)
            template = env.get_template(path)
            return template.render(obj=self)
        except jinja2.TemplateNotFound:
            return None


class Sample(Html):
    def __init__(self, name, description=None, properties=None):
        self.name = name
        self.description = description
        self.properties = properties or {}


class Experiment(object):
    def __init__(self, date, description=None):
        self.date = date
        self.description = description

    def _repr_html_(self):
        template = env.get_template('core.Experiment.html')
        return template.render(obj=self)


class Measurement(object):
    def __init__(self, data, params):
        self.params = params
        self._data = data

    def __getitem__(self, item):
        return self._data[item]

    def _repr_html_(self):
        template = env.get_template('core.Measurement.html')
