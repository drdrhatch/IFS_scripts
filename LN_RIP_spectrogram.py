#!/usr/bin/env python
# -*- coding: utf-8 -*-
#Created by Max T. Curie 09/22/2020

from LN_tools import get_suffix
from LN_tools import B1_ky_f_spectrum_Z_sum

#*****************************************************************
#*******************Beginning of the User block*******************

max_Z0=0.035    #in the unit of meter
min_Z0=-0.035   #in the unit of meter

Outboard_mid_plane=False  #change to True if one wants to only want to look at outboard mid-plane
plot=True
show=True
csv_output=True

pic_path='pic'        #path one want to store the picture and video in
csv_path='csv'        #path one want to store the picture and video in
iterdb_file_name='/global/u1/m/maxcurie/max/Cases/DIIID175823_250k/DIIID175823.iterdb'

#*********************End of the User block***********************
#*****************************************************************

suffix=get_suffix()
uni_freq,amplitude_frequency_uni_ky_sum,amplitude_frequency_uni=\
    B1_ky_f_spectrum_Z_sum(suffix,iterdb_file_name,\
        min_Z0,max_Z0,Outboard_mid_plane,\
        plot,show,csv_output,pic_path,csv_path)

