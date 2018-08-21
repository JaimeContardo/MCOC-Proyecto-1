import scipy as sp
import matplotlib.pyplot as pyplot
from scipy.integrate import solve_ivp
def int_arias(ruta):	
	a = sp.loadtxt(ruta)
	Nt= a.size
	dt = 1./200

	t = sp.arange(0,dt*Nt,dt)
	Ia = sp.zeros(Nt)
	a2 = a**2
	da2 = (a2[0:-1]+a2[1:])*dt/2
	Ia[1:]=sp.cumsum(da2)*sp.pi/(2*9.81)

	Ia_inf = Ia.max()

	i_05 = sp.argmin (abs(Ia-(Ia_inf*0.05)))
	i_95 = sp.argmin (abs(Ia-(Ia_inf*0.95)))

	t_05 = t[i_05]
	t_95=t[i_05]
	Ia_05 = Ia[i_05]

	t_95 = t[i_95]
	Ia_95 = Ia[i_95]
	D = t_95-t_05 
	#plt.subplot(2,1,1)
	#plt.plot(t,a)
	#plt.subplot(2,1,2)
	#plt.plot(t,Ia)
	#plt.show()
	return D

def PGV_PGD(ruta):
	a = sp.loadtxt(ruta)
	Nt = a.size
	dt =1./200
	t = sp.arange(0,dt*Nt,dt)
	int_aceleracion = sp.zeros(Nt)
	da = (a[0:-1]+a[1:])*dt/2
	pgv = max(da)
	V = da
	Nt = V.size
	dt =1./200
	t = sp.arange(0,dt*Nt,dt)
	int_velocidad = sp.zeros(Nt)
	dv = (V[0:-1]+V[1:])*dt/2
	pgd = max(dv)
	pgv_pgd = [pgv,pgd]
	return pgv_pgd

