#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 17:02:49 2018

@author: jose
"""

import os
import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt
import linecache
import glob

ruta=os.path.dirname(os.path.abspath(__file__)) #Obtengo ruta actual de donde esta el archivo
rutaRegistros=str(ruta)+'/registross/*.txt'
print rutaRegistros
def interpol(ruta,graficar=False):
	nombre=linecache.getline(ruta, 1)[19:39]
	nombre+= "\n"
	linea=linecache.getline(ruta, 2)
	linea=linea[20:]
	dt=1./float(linea[:linea.find("m")])

	linea=linecache.getline(ruta, 3)

	n=int(linea[28:])
	d=dt*n
	a = sp.loadtxt(ruta)
	ti = sp.arange(0,d,dt)
	from scipy.integrate import simps
	Il = simps(a**2,ti)
	int_arias = (sp.pi/2*9.81)*Il
    	
    

	a=sp.loadtxt(ruta,comments='#')
	linea=linecache.getline(ruta, 4)
	estacion=linea[12:16]
	print estacion

	t=sp.arange(0,dt*n,dt)
	i=interpolate.interp1d(t,a,kind='cubic')
	inter=i(t)



	if graficar:
		plt.plot(t,a,'o',t,inter,'-')
		plt.show()

	amax=max(abs(a/9.81))
	print amax    
    

	return {"funcion":i,"amax":amax,"tmax":dt*n, "nombre":nombre,"estacion":estacion}

print 'Eligiendo registros que cumplan 0.55<amax>0.35, guardandolos en .../registros/seleccion.txt '
listaRutas=(glob.glob(rutaRegistros))
arch=open('registross/seleccion.txt','w')
count=0
for route in listaRutas:
	try:
		i=interpol(route)
		if i['amax']>0.2 and i['amax']<0.8:
			arch.write(i['nombre']+i['estacion'])
			count+=1
	except:
		"Error"
print "Seleccionados {} registros".format(count)
arch.close()
metadatos = {}
metadatos['Fecha']='2012-05-19'
metadatos['Hora']='08:35:08.998393'
metadatos['Epi_Lat']=-25.74
metadatos['Epi_Lon']=-70.86
metadatos['Epi_Profundidad']=83
metadatos['M']=6.1
metadatos['Estacion_Lat']=-25.163
metadatos['Estacion_Lon']=-69.590
metadatos['Estacion_Nombre']='G002'
metadatos['Componente']='HNE'
metadatos['PGA']=0
metadatos['PGV']=0
metadatos['PGD']=0
metadatos['Duracion']=0
for i in metadatos:
    print str(i)+' : '+str(metadatos[i])
    
