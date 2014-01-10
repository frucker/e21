# -*- coding: utf-8 -*-
#
# e21, (c) 2013, see AUTHORS. Licensed under the GNU GPL.
from e21.sweet16.loader import Loader
import e21.sweet16.susceptibility

Loader.CREATOR['susceptibility'] = e21.sweet16.susceptibility.create

del create
