import linecache
import scipy as sp
from IA import PGV_PGD, int_arias

#Parametros : ruta (es la ruta del archivo), c (un contador para guardar el archvio creado)
def Base_Datos_Sismicos(ruta,c):
    metadatos = {}
    
    linea = linecache.getline(ruta,1)
    fecha = linea[20:30]
    hora = linea[31:39]
    
    linea = linecache.getline(ruta,5)
    est_Lat = float(linea [11:18])
    est_Lon = float(linea[29:])
    
    linea = linecache.getline(ruta,7)
    epi_Lat = float(linea [15:21])
    epi_Lon = float(linea[36:])
    
    linea = linecache.getline(ruta,8)
    Profundidad = float(linea [15:17])
    Magnitud = float(linea[28:])
    
    linea = linecache.getline(ruta,4)
    nombre = linea[12:16]
    componente = linea[29:32]
    
    linea = linecache.getline(ruta,2)
    tasa = float(linea[20:25])
    
    linea = linecache.getline(ruta,3)
    cantidad = float(linea[28:])
    
    t = sp.arange(0,float(cantidad/tasa),float(1/tasa))
    a = sp.loadtxt(ruta,comments='#')
                   
    print 'Vector de Tiempo : ', t[0:],'\n'
    print 'Vector de Acelearaciones : ', a[0:],'\n'
    
    amax=max(abs(a)/9.81)*100
    
    metadatos['Fecha']= fecha
    metadatos['Hora']= hora
    metadatos['Epi_Lat']=epi_Lat
    metadatos['Epi_Lon']=epi_Lon
    metadatos['Epi_Profundidad']=Profundidad
    metadatos['M']=Magnitud
    metadatos['Estacion_Lat']= est_Lat
    metadatos['Estacion_Lon']= est_Lon
    metadatos['Estacion_Nombre']= nombre
    metadatos['Componente']= componente
    metadatos['PGA']= amax
    metadatos['PGV']=PGV_PGD(ruta)[0]
    metadatos['PGD']=PGV_PGD(ruta)[1]
    metadatos['Duracion']=int_arias(ruta) #Esto lo hice segun la cantidad de datos muestreados y la frecuencia
    
    for i in metadatos:
        print str(i)+' : '+str(metadatos[i])+'\n'
    
    sp.savez('registro_{}.npz'.format(c),metadatos=metadatos,t=t,a=a)
    
    return metadatos, t , a

Base_Datos_Sismicos('/home/jose/Documents/registross/20130130-201540-GO03-HNE.txt',01)
# Para obtener los metadatos, el vector t o el vector a, se debe seguir esta forma:
# metadatos , t , a = Base_Datos_Sismicos(ruta)
    
