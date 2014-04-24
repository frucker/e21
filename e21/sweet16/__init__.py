# -*- coding: utf-8 -*-
#
# e21, (c) 2013, see AUTHORS. Licensed under the GNU GPL.
from e21.sweet16.loader import Loader
from e21.sweet16.susceptibility import create as susceptibility_create
from e21.sweet16.transport import create as transport_create

Loader.CREATOR['susceptibility'] = susceptibility_create
Loader.CREATOR['transport'] = transport_create

del transport_create
del susceptibility_create
