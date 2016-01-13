# -*- coding: utf-8 -*-
#
# e21, (c) 2013, see AUTHORS. Licensed under the GNU GPL.
from e21.Nicos.loader import Loader
from e21.Nicos.Mira.miraHe3 import create as mira_create

Loader.CREATOR['miraHe3'] = mira_create

del mira_create
