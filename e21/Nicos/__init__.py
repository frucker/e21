# -*- coding: utf-8 -*-
#
# e21, (c) 2013, see AUTHORS. Licensed under the GNU GPL.
from e21.Nicos.loader import Loader
from e21.Nicos.Dryo.susceptibility import create as dryo_create

Loader.CREATOR['susceptibility'] = dryo_create

del dryo_create
