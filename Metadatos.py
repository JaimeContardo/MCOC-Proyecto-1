#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 17:36:51 2018

@author: jaimecontardo
"""
import linecache

def diccionario_metadatos(ruta):
    metadatos = {}
    linea = linecache.getline(ruta,1)
    fecha = linea[20:30]
    hora = linea[31:39]
    linea = linecache.getline(ruta,5)
    e_Lat = float(linea [11:18])
    e_Lon = float(linea[29:])
    linea = linecache.getline(ruta,4)
    nombre = linea[12:16]
    componente = linea[29:32]
    linea = linecache.getline(ruta,2)
    tasa = float(linea[20:25])
    linea = linecache.getline(ruta,3)
    cantidad = float(linea[28:])
    metadatos['Fecha']= fecha
    metadatos['Hora']= hora
    metadatos['Epi_Lat']='****'
    metadatos['Epi_Lon']='****'
    metadatos['Epi_Profundidad']='****'
    metadatos['M']='****'
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


    
