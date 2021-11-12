#!/usr/bin/env python

def customizeEcalSev(process):
    process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")

    return process
