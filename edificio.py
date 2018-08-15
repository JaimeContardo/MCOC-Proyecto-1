# Metodos Computacionales en Ingenieria
# Facultad de Ingenieria y Ciencias Aplicadas - UANDES
# Segundo Semestre 2018
# Prof. J. Abell (jaabell@miuandes.cl)
#
# Proyecto 1 - Diseno de sistema de disipacion de energia en edificio 
# 
# Ejemplo 1  
# Lenguaje Python 2.7
#
# En este ejemplo se demuestra como escribir un integrador de ecuaciones diferenciales
# de segundo orden, usando el concepto de espacio de estado y el metodo de Euler. Se compara
# la solucion del metodo de Euler con el integrador RK45 de Scipy, via la interfaz solve_ivp.
#
# Se aplica el integrador a resolver la respuesta de un oscilador de 1 GDL a una condicion 
# inicial (sin forzamiento) considerando amortiguamiento lineal y un modelo de 
# amortiguamiento friccional. 
# =======================================================================================

#Importar librerias relevantes
import scipy as sp
import numpy as np

import matplotlib.pyplot as plt
# import scipy.integrate
from scipy.integrate import solve_ivp

from regis import *



datos=sp.load('mck.npz')
#Parametros

M= np.matrix(datos['Mmatriz'])      # [kg]  Masa del oscilador
K=np.matrix(datos['Kmatriz'])        # [N/m] Rigidez del oscilador       # [N]   Capacidad maxima friccional del disipador
CAP=datos['CAPvector']
C=sp.matrix(datos['Cmatriz'])
# C=np.zeros((20,20))


vr = 0.0001     # [m/s] Velocidad de referencia para la aproximacion de la friccion via tanh. 
d0 = 10         # [m]   Condicion inicial de desplazamiento
v0 = 0          # [m/s] Condicion inicial de velocidad
dt = 0.001       # [s]   Paso de integracion a usar
tmax = 600  # [s]   Tiempo maximo de integracion 


#Definimos la matriz de estado A
Mi=M.I
A = np.block(
    [
    [np.zeros((20,20)),np.eye(20)],
    [-sp.matmul(Mi,K),-sp.matmul(Mi,C)]
    ]
    )

inte=interpol('/home/felipe/Desktop/MCOC-Proyecto-1/20180702-092950-PB15-HLE.txt')[0]

#Definimos la funcion del lado derecho de la EDO de primer orden
#a resolver zp = fun(t,z).  zp es la derivada temporal de z.
def fun(t,z):
    #---- Esto sirve solo para reportar en que tiempo vamos
    if t > fun.tnextreport:
        print "  {} at t = {}".format(fun.solver, fun.tnextreport)
        fun.tnextreport += 1

    #------- Lo que sigue es lo que importa
    if fun.solver=="RK45":
        z=np.squeeze(z)

    Famort = sp.zeros(40)   #vector de fuerza friccional de amortiguamiento
    
    # tangh=sp.zeros(20)
    Famort[0]=- (
        CAP[0] * (1./M[0,0]) * sp.tanh( (z[20]/vr) )
        )

    # tangh[0]=-Mi[0,0]*CAP[0]*sp.tanh(z[20]/vr)
    for i in range(1,20):
        Famort[i]=-(1./M[i,i])*CAP[i]*sp.tanh((z[i+20]-z[i-1+20])/vr)

    Ft=sp.zeros(40)
    if t<226.81:
        Ft[20:]=inte(t)
    # Famort[20:]=tangh

    # Pueba
    # print t
    # Famort[20:]=-sp.matmul(Mi,tangh)
    # Prueba

    # CT=CAP*tangh
    # CT=np.matrix(CT).T

    # new=(Mi*CT).T
    # Famort[20:]= -new
    # Famort=np.matrix(Famort)
    
    # print sp.matmul(A,z)+Famort

    # q=raw_input()
    # return np.squeeze(np.asarray(sp.matmul(A,z) + Famort)) #+ funcion de fuerzas de t
    return sp.matmul(A,z)+Famort+Ft
#Vector de tiempo y su largo (Nt)
t = sp.arange(0, tmax, dt)
Nt = len(t)



#Inicializar una matriz z para guardar las solucion discretizada
z_euler = sp.zeros((40,Nt+1))
z_RK45 = sp.zeros((40,Nt+1))

#Condicion inciial en t = 0, i = 0. 
z0 = sp.array([
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                v0,v0, v0,v0, v0,v0, v0,v0, v0,v0, v0,v0, v0,v0, v0,v0, v0,v0, v0,v0
            ])


z_euler[:,0] = z0
z_RK45[:,0] = z0


print "Integrando con Euler"
fun.tnextreport = 0
fun.solver = "Euler"
   
i = 1
ti = dt 
while (ti < tmax):
    z_euler[:,i] = dt * fun(ti, z_euler[:,i-1]) + z_euler[:,i-1]
    ti += dt
    i += 1



print "Integrando con RK45"
fun.tnextreport = 0
fun.solver = "RK45"
solucion_rk45 = solve_ivp(fun, [0., tmax], z0, method='RK45', t_eval=t, vectorized=False )
z_RK45[:,1:] = solucion_rk45.y
# z_RK45=z_euler  



#Graficar solucion en desplazamiento y velocidad para ambos metodos
plt.figure()

for z, lab in zip([z_euler, z_RK45], ["Euler", "RK45"]):

    u = z[19,:]
    v = z[39,:]

    #Extraer desplazamientos y velocidades
    u = z[19,:-1]
    v = z[39,:-1]
    dmax=max(abs(u))
    plt.subplot(2,1,1)
    plt.plot(t, u, label=lab)
    plt.ylim([-1.5*dmax, 1.5*dmax])
    plt.xlim([0, tmax])
    plt.ylabel("Desplazamiento, $u = z_1$ (m)")
    plt.grid(True)

    vmax = max(abs(v))
    plt.subplot(2,1,2)
    plt.plot(t, v)
    plt.ylabel("Velocidad, $\dot{u} = z_2$ (m/s)")
    plt.xlabel("Tiempo, $t$ (s)")
    plt.ylim([-1.5*vmax, 1.5*vmax])
    plt.xlim([0, tmax])
    plt.grid(True)

plt.subplot(2,1,1)
plt.legend()
plt.suptitle("Solucion por metodo de Euler")

plt.show()
