"""
Created on Sun Aug 19 22:26:07 2018

@author: jose
"""
import linecache
import scipy as sp
import matplotlib.pyplot as plt

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
    #Esto lo hice segun la cantidad de datos muestreados y la frecuencia
    
    for i in metadatos:
        print str(i)+' : '+str(metadatos[i])+'\n'
    
    sp.savez('registro_{}.npz'.format(c),metadatos=metadatos,t=t,a=a)
    a = sp.loadtxt(ruta)
 
    dt = 1./200.
    Nt = a.size
     
     
    t = sp.arange(0, dt*Nt, dt)
     
    Ia = sp.zeros(Nt)
    v = sp.zeros(Nt)
    d = sp.zeros(Nt)
     
    v[1:] = sp.cumsum(a[1:] + a[0:-1])*dt/2
    d[1:] = sp.cumsum(v[1:] + v[0:-1])*dt/2
     
     
     
    g = 9.806
     
    a2 = a**2
    da2 = (a2[0:-1] + a2[1:])*dt/2
     
    Ia[1:] = sp.cumsum(da2)*sp.pi/(2*g)
     
    Ia_inf = Ia.max()
     
    i_PGA = sp.argmax(abs(a))
    t_PGA = t[i_PGA]
    PGA = (a[i_PGA])
     
    i_PGV = sp.argmax(abs(v))
    t_PGV = t[i_PGV]
    PGV = (v[i_PGV])
     
    i_PGD = sp.argmax(abs(d))
    t_PGD = t[i_PGD]
    PGD = (d[i_PGD])
     
     
    i_05 = sp.argmin( abs(Ia - 0.05*Ia_inf) )
    i_95 = sp.argmin( abs(Ia - 0.95*Ia_inf) )
     
    t_05 = t[i_05]

     
    t_95 = t[i_95]
   
     
    D_5_95 = t_95 - t_05
     
    
    metadatos['PGV']=PGV
    metadatos['PGD']=PGD
    metadatos['Duracion']= D_5_95
     
    plt.figure().set_size_inches([9,6])
     
    plt.subplot(3,1,1)
    plt.plot(t,a/g)
    plt.axvline(t_05, color="k", linestyle="--")
    plt.axvline(t_95, color="k", linestyle="--")
    plt.text(t_PGA,PGA/g, "PGA={0:0.3f}g".format(abs(PGA)/g))
    plt.plot(t_PGA,PGA/g, "ob")
    plt.ylim([-0.8,0.8])
    plt.grid(True)
    plt.ylabel("Acc, $a$ (g)")
     
    plt.subplot(3,1,2)
    plt.plot(t,v*100)
    plt.axvline(t_05, color="k", linestyle="--")
    plt.axvline(t_95, color="k", linestyle="--")
    plt.text(t_PGV,PGV*100, "PGV={0:0.3f}cm/s".format(abs(PGV)*100))
    plt.plot(t_PGV,PGV*100, "ob")
    plt.ylim([-15,15])
    plt.grid(True)
    plt.ylabel("Vel, $v$ (cm/s)")
     
     
    plt.subplot(3,1,3)
    plt.plot(t,d*100)
    plt.axvline(t_05, color="k", linestyle="--")
    plt.axvline(t_95, color="k", linestyle="--")
    plt.text(t_PGD,PGD*100, "PGD={0:0.3f}cm".format(abs(PGD)*100))
    plt.plot(t_PGD,PGD*100, "ob")
    plt.ylim([-15,15])
    plt.grid(True)
    plt.xlabel("Tiempo, $t$ (s)")
    plt.ylabel("Dis, $d$ (cm)")
     
    plt.subplot(3,1,1)
    plt.title(ruta[:-4] + "    $D_{{5-95}} = {0:5.2f}$s".format(D_5_95))
     
    plt.tight_layout()
     
    plt.show()
    
    return metadatos, t , a
    
Base_Datos_Sismicos('/home/jose/Documents/registross/20160210-003304-CO06-HNN.txt',21)
# Para obtener los metadatos, el vector t o el vector a, se debe seguir esta forma:
# metadatos , t , a = Base_Datos_Sismicos(ruta)
# Para obtener los metadatos, el vector t o el vector a, se debe seguir esta forma:
# metadatos , t , a = Base_Datos_Sismicos(ruta)
    
