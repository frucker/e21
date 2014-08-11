# -*- coding: utf-8 -*-
#
# e21, (c) 2013, see AUTHORS. Licensed under the GNU GPL.
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='e21',
    version=__import__('e21').__version__,
    packages=['e21'],
    include_package_data=True
)
