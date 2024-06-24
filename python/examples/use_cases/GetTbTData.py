#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 10:15:54 2024

@author: Teresia Olsson, teresia.olsson@helmholtz-berlin.de
"""

import epics
from functools import lru_cache
import time
import pandas as pd
#import datetime as dt

class GetTbTData():
    # Private.

    def __init__(self, bpm_list):

        self._online    = False
        self._bpm_list  = bpm_list
        self._data_list = ['X', 'Y', 'MT']

        # Initialise the BPMs
        for bpm in self._bpm_list:
            epics.caput(bpm+':signals:ddc_synth.SCAN', 'I/O Intr')
            epics.caput(bpm+':signals:ddc_synth.ACQM', 'Event')
            epics.caput(bpm+':signals:ddc_synth.OFFS', -100)

        # What does this do?
        def callbackGetTbt(pvname, value, **kwargs):
            global wait_for_data
            wait_for_data = False

        timestamp = epics.PV(bpm_list[0]+':signals:ddc_synth.MT')
        timestamp.add_callback(callbackGetTbt)

    # Public.

    def getTbTdata(self):
        lru_cache(maxsize=None)

        global wait_for_data
        # wait for trigger
        wait_for_data = True
        while wait_for_data:
            time.sleep(0.1)

        def epics_PV(pv_name: str) -> epics.PV:
            return epics.PV(pv_name)

        result = pd.DataFrame([])
        for bpm in self._bpm_list:
            xpv = epics_PV(bpm+':signals:ddc_synth.X')
            ypv = epics_PV(bpm+':signals:ddc_synth.Y')
            tpv = epics_PV(bpm+':signals:ddc_synth.MT')

            x = xpv.get()#_with_metadata()
            y = ypv.get()#_with_metadata()
            t = tpv.get()#_with_metadata()

        actresult = \
            pd.DataFrame([[bpm, x, y, t]], columns=["BPM_name", 'X', 'Y', "MT"])
        result = pd.concat([result, actresult], axis=1)

        return result 
