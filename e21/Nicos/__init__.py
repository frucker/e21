# -*- coding: utf-8 -*-
#
# e21, (c) 2013, see AUTHORS. Licensed under the GNU GPL.
from e21.Nicos.loader import Loader
from e21.Nicos.susceptibility import create as Nicos_create

Loader.CREATOR['susceptibility'] = Nicos_create

del Nicos_create
