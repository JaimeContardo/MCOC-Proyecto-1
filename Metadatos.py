#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 17:36:51 2018

@author: jaimecontardo
"""
import linecache
import numpy as np
import scipy as sp

def Base_Datos_Sismicos(ruta):
    metadatos = {}
    linea = linecache.getline(ruta,1)
    fecha = linea[20:30]
    hora = linea[31:39]
    linea = linecache.getline(ruta,5)
    e_Lat_estacion = float(linea [11:18])
    e_Lon_estacion = float(linea[29:])
    linea = linecache.getline(ruta,7)
    e_Lat_epi = float(linea [11:18])#Cambiar el rango de lectura
    e_Lon_epi = float(linea[29:])#Cambiar el rango de lectura
    linea = linecache.getline(ruta,8)
    Profundidad = float(linea [11:18])#Cambiar el rango de lectura
    Magnitud = float(linea[29:])#Cambiar el rango de lectura
    linea = linecache.getline(ruta,4)
    nombre = linea[12:16]
    componente = linea[29:32]
    linea = linecache.getline(ruta,2)
    tasa = float(linea[20:25])
    linea = linecache.getline(ruta,3)
    cantidad = float(linea[28:])
    metadatos['Fecha']= fecha
    metadatos['Hora']= hora
    metadatos['Epi_Lat']=e_Lat_epi
    metadatos['Epi_Lon']=e_Lon_epi
    metadatos['Epi_Profundidad']=Profundidad
    metadatos['M']=Magnitud
    metadatos['Estacion_Lat']= e_Lat
    metadatos['Estacion_Lon']= e_Lon
    metadatos['Estacion_Nombre']= nombre
    metadatos['Componente']= componente
    metadatos['PGA']='****'
    metadatos['PGV']='****'
    metadatos['PGD']='****'
    metadatos['Duracion']=cantidad/tasa #Esto lo hice segun la cantidad de datos muestreados y la frecuencia
    for i in metadatos:
        print str(i)+' : '+str(metadatos[i])
    t = sp.arange(0,float(cantidad/tasa),float(1/tasa))
    # agregar la lectura de vector de aceleraciones
    # falta imprimir vectores tiempo y aceleracion 
    return metadatos


    
