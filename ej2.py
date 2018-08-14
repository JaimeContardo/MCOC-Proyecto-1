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


datos=sp.load('mck.npz')
#Parametros

M= np.matrix(datos['Mmatriz'])      # [kg]  Masa del oscilador
K=np.matrix(datos['Kmatriz'])        # [N/m] Rigidez del oscilador       # [N]   Capacidad maxima friccional del disipador
CAP=datos['CAPvector']
C=sp.matrix(datos['Cmatriz'])

vr = 0.001     # [m/s] Velocidad de referencia para la aproximacion de la friccion via tanh. 
d0 = 10         # [m]   Condicion inicial de desplazamiento
v0 = 0          # [m/s] Condicion inicial de velocidad
dt = 0.01       # [s]   Paso de integracion a usar
tmax = 200    # [s]   Tiempo maximo de integracion 


#Definimos la matriz de estado A
Mi=M.I
A = np.block(
    [
    [np.zeros((20,20)),np.eye(20)],
    [-Mi*K,-Mi*C]
    ]
    )

#Definimos la funcion del lado derecho de la EDO de primer orden
#a resolver zp = fun(t,z).  zp es la derivada temporal de z.
def fun(t,z):
    #---- Esto sirve solo para reportar en que tiempo vamos
    if t > fun.tnextreport:
        print "  {} at t = {}".format(fun.solver, fun.tnextreport)
        fun.tnextreport += 1

    #------- Lo que sigue es lo que importa
    Famort = sp.zeros(40)   #vector de fuerza friccional de amortiguamiento
    
    tangh=sp.zeros(20)
    tangh[0]=sp.tanh(z[20]/vr)
    
    for i in range(1,20):
        tangh[i]=sp.tanh((z[i+20]-z[i-1+20])/vr)

    CT=CAP*tangh
    CT=np.matrix(CT).T

    new=(Mi*CT).T
    Famort[20:]= -new
    Famort=np.matrix(Famort)
    
    # print sp.matmul(A,z)+Famort

    # q=raw_input()
    return np.squeeze(np.asarray(sp.matmul(A,z) + Famort)) #+ funcion de fuerzas de t
 
#Vector de tiempo y su largo (Nt)
t = sp.arange(0, tmax, dt)
Nt = len(t)



#Inicializar una matriz z para guardar las solucion discretizada
z_euler = sp.zeros((40,Nt+1))
z_RK45 = sp.zeros((40,Nt+1))

#Condicion inciial en t = 0, i = 0. 
z0 = sp.array([
                0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,7.5,8.,8.5,9.,9.5,
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
# solucion_rk45 = solve_ivp(fun, [0., tmax], z0, method='RK45', t_eval=t, vectorized=False )
# z_RK45[:,1:] = solucion_rk45.y
z_RK45=z_euler  



#Graficar solucion en desplazamiento y velocidad para ambos metodos
plt.figure()

for z, lab in zip([z_euler, z_RK45], ["Euler", "RK45"]):

    u = z[0,:]
    v = z[20,:]

    #Extraer desplazamientos y velocidades
    u = z[0,:-1]
    v = z[20,:-1]

    plt.subplot(2,1,1)
    plt.plot(t, u, label=lab)
    plt.ylim([-1.5*d0, 1.5*d0])
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
