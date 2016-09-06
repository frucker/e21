# -*- coding: utf-8 -*-
#
# e21, (c) 2013, see AUTHORS. Licensed under the GNU GPL.
from e21.mischer.loader import Loader
from e21.mischer.vcm import create as vcm_create

Loader.CREATOR['vcm'] = vcm_create

del vcm_create