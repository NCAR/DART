#!/usr/bin/env python
# -*- coding: utf-8 -*-


#little tool to generate chemistry set to perturb the boundary condition with a normal law
#you should call the script this way: ./random.py spread ensemble_size
#also specify the good filew and filer pathways
#author: Jerome Barre barre@ucar.edu  (2013)
# modified by Arthur P. Mizzi
#
# to exit use exit()
#
import numpy as np
from sys import argv
#
spread=float(argv[1])
nens=int(argv[2])
dart_pert_dir=str(argv[3])
run_dir=str(argv[4])
#
z_stdn=2.58
z_stnd=1.96
#
filew = open(run_dir+"/set0", "r")
for i in range(nens):
   filer = open(run_dir+"/set"+str(i+1),"w")
   filer.write('spc_map =\n')
   for lig in file(run_dir+'/set0'):
      ligne=lig.split()
      if len(ligne)==3:
         moz=ligne[0]
         wrf=ligne[2]
         icnt=0
         coeff=np.random.normal(0,1,1)            
         coef=1+coeff*spread
         while -z_stnd>coeff or z_stnd<coeff or coef<=0:
            icnt+=1
            if icnt>10:
               print 'Tail cutoff error '
               exit()
            coeff=np.random.normal(0,1,1)            
            coef=1+coeff*spread
         if str(moz)=="'BC1" or str(moz)=="'BC2" or str(moz)=="'OC1" or str(moz)=="'OC2" or str(moz)=="'SEAS_1" or str(moz)=="'SEAS_2" or str(moz)=="'SEAS_3" or str(moz)=="'SEAS_3" or str(moz)=="'SEAS_4" or str(moz)=="'DUST_1" or str(moz)=="'DUST_2" or str(moz)=="'DUST_5":
            wrf_dec=wrf.split('*')
            coef=coef*float(wrf_dec[0])
            print '     '+moz+' -> '+'%.2f'%float(coef)+'*'+wrf_dec[1]
            filer.write('     '+moz+' -> '+'%.2f'%float(coef)+'*'+wrf_dec[1]+'\n')
         elif str(moz)=="'DUST_3" or str(moz)=="'DUST_4":
            wrf_el=wrf.split('+')
            wrf_dec1=wrf_el[0].split('*')
            coef1=coef*float(wrf_dec1[0])
            wrf_dec2=wrf_el[1].split('*')
            coef2=coef*float(wrf_dec2[0])
            print '     '+moz+' -> '+'%.2f'%float(coef1)+'*'+wrf_dec1[1]+'+'+'%.2f'%float(coef2)+'*'+wrf_dec2[1]
            filer.write('     '+moz+' -> '+'%.2f'%float(coef1)+'*'+wrf_dec1[1]+'+'+'%.2f'%float(coef2)+'*'+wrf_dec2[1]+'\n')
         else:
            print '     '+moz+' -> '+'%.2f'%float(coef)+'*'+wrf
            filer.write('     '+moz+' -> '+'%.2f'%float(coef)+'*'+wrf+'\n')
   filer.write('/')
   filer.close()
filew.close()
