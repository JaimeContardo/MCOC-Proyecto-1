import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt
import linecache
import glob

ruta='/home/felipe/Desktop/MCOC-Proyecto-1/registros/20180702_092951_PB15/20180702-092950-PB15-HLE.txt'

def interpol(ruta,graficar=False):

	linea=linecache.getline(ruta, 2)
	linea=linea[20:]
	dt=1./float(linea[:linea.find("m")])

	linea=linecache.getline(ruta, 3)

	n=int(linea[28:])

	a=sp.loadtxt(ruta,comments='#')


	t=sp.arange(0,dt*n,dt)
	i=interpolate.interp1d(t,a,kind='cubic')
	inter=i(t)

	if graficar:
		plt.plot(t,a,'o',t,inter,'-')
		plt.show()

	amax=max(abs(a))

	return i,amax



listaRutas=(glob.glob("/home/felipe/Desktop/MCOC-Proyecto-1/registros/20180702_092951_PB15/*.txt"))

for ruta in listaRutas:
	i=interpol(ruta)
	if i[1]>0.2:
		print i[1]