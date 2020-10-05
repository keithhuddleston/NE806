#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 13:33:45 2020

@author: keith
"""
import openmc as mc

fn = 'U238endf.txt'

u238_endf = mc.data.endf.Evaluation(fn)
print(u238_endf)
f = open(fn)

# rm_resonance = u238_endf.resonances.ranges[0]
# covariance = u238_endf.resonance_covariance.ranges[0].covariance
# f = open('U238endf.txt', 'rt')
# mc.data.endf.get_cont_record(f)