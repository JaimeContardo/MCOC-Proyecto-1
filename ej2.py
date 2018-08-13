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

#Parametros
xi = 0          # Razon de amortiguamiento critico
m = 1.  
M=np.matrix([[m,0],[0,m]])        # [kg]  Masa del oscilador
k = 10. 
K=np.matrix([[2*k,-k],[-k,k]])        # [N/m] Rigidez del oscilador
Cap = 20.       # [N]   Capacidad maxima friccional del disipador
CAP=np.matrix([[2*Cap,-Cap],[-Cap,Cap]])
vr = 0.0001     # [m/s] Velocidad de referencia para la aproximacion de la friccion via tanh. 
d0 = 10         # [m]   Condicion inicial de desplazamiento
v0 = 0          # [m/s] Condicion inicial de velocidad
dt = 0.01       # [s]   Paso de integracion a usar
tmax = 10       # [s]   Tiempo maximo de integracion 

#Calcular algunos parametros dinamicos relevantes
wn = sp.sqrt(k/m)
fn = wn/2/sp.pi
Tn = 1/fn
c = 2*xi*wn*m    #constante de amortiguamiento lineal

#Reportar
print "wn = ",wn
print "fn = ",fn
print "Tn = ",Tn


#Definimos la matriz de estado A
A = np.block([
    [np.zeros((2,2)),np.eye(2)],
    [-M.I*K,np.zeros((2,2))]
    ])

#Definimos la funcion del lado derecho de la EDO de primer orden
#a resolver zp = fun(t,z).  zp es la derivada temporal de z.
def fun(t,z):
    #---- Esto sirve solo para reportar en que tiempo vamos
    if t > fun.tnextreport:
        print "  {} at t = {}".format(fun.solver, fun.tnextreport)
        fun.tnextreport += 1
    #------- Lo que sigue es lo que importa
    Famort = sp.zeros((4))   #vector de fuerza friccional de amortiguamiento
    Famort[2:] = -sp.tanh(z[2:]/vr) #N/2+1:

    Famort[2:]=sp.matmul(-M.I*CAP,Famort[2:])
    return sp.matmul(A,z) + Famort
 
#Vector de tiempo y su largo (Nt)
t = sp.arange(0, tmax, dt)
Nt = len(t)



#Inicializar una matriz z para guardar las solucion discretizada
z_euler = sp.zeros((4,Nt+1))
z_RK45 = sp.zeros((4,Nt+1))

#Condicion inciial en t = 0, i = 0. 
z0 = sp.array([d0, d0, v0,v0])


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




#Graficar solucion en desplazamiento y velocidad para ambos metodos
plt.figure()

for z, lab in zip([z_euler, z_RK45], ["Euler", "RK45"]):

    u = z[0,:]
    v = z[1,:]

    #Extraer desplazamientos y velocidades
    u = z[0,:-1]
    v = z[1,:-1]

    plt.subplot(2,1,1)
    plt.plot(t, u, label=lab)
    plt.ylim([-1.5*d0, 1.5*d0])
    plt.xlim([0, tmax])
    plt.ylabel("Despazamiento, $u = z_1$ (m)")
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
