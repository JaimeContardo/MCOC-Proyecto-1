from numpy import *
import numpy as np
import scipy as sp

# 
# 
# 
masa = []
largos = [12.*3.5,12.*3.5,12.*3.5,12.*3.5,8.*3.5,8.*3.5,8.*3.5,8.*3.5,4.*3.5,4.*3.5,4.*3.5,4.*3.5,4.*3.5,4.*3.5,4.*3.5,4.*3.5,4.*3.5,4.*3.5,4.*3.5,4.*3.5]
i = 0 
while i<20:
    masa.append(6.2*largos[i]*1000)
    i+=1
M =  np.transpose(masa)*np.identity(20)
#print masa #Vector de masas

print M #matriz de masa

# 
# 
# 

columnas = []
medidas= array([600.,700.,800.,900.,1000.])
largos = [4000.,2800.]
youngmod = 23.5 #KN/mm2
c = 0
while c<5:
    columnas.append(12.*youngmod*(medidas[c])**4/12.) #KN/mm
    c+=1
columnas = array(columnas)
pisos = array([[6.,2.,0.,1.,4.],[4.,0.,1.,0.,4.],[0.,0.,1.,4.,0.],[0.,1.,4.,0.,0.],[5.,0.,0.,0.,0.]])
k_pisos = np.zeros(20)

k_pisos[0:4] = np.sum((columnas*pisos[0]))
k_pisos[4:8] = np.sum((columnas*pisos[1]))
k_pisos[8:12] = np.sum((columnas*pisos[2]))
k_pisos[12:16] = np.sum((columnas*pisos[3]))
k_pisos[16:20] = np.sum((columnas*pisos[4]))

k_pisos[0]/= (largos[0]**3)
k_pisos[1:]/=(largos[1]**3)


kmatrix = np.zeros((20,20)) 
#print k_pisos #Vector de rigidez por piso

i = 0
while i < 19: 
    kmatrix[i][i] = k_pisos[i] + k_pisos[i+1] 
    kmatrix[i][i+1] = -k_pisos[i+1] 
    kmatrix[i+1][i] = -k_pisos[i+1] 
    i += 1
kmatrix[19][19] = k_pisos[19]
kmatrix*=10**6 #N/M
print kmatrix

# 
# 
# 

f1 = 0.2 
f2 = 2.
xi = 0.025
a0 = (4*np.pi*f1*f2*(f1-f2)*xi)/(f1**2-f2**2)
a1 = xi*(f1-f2)/(np.pi*(f1**2-f2**2))
c= a0*M+a1*kmatrix


print c

# 
# 
# 

caps= array([150.,250.,500.,800.]) #kN, capacidades disponibles
c_por_piso=np.zeros(20)
# c_por_piso=np.array([0,caps[3],caps[1],caps[2],2*caps[0],caps[3],0,2*caps[0],caps[1],caps[0],
#                      caps[2],0,2*caps[0],caps[0],2*caps[0],caps[1],caps[1],caps[0],caps[0],caps[0]])#cantidad de disipadores y que tipo por piso
c_por_piso*=1000 #a N

print "Capacidad total instalada: " + str(np.sum(c_por_piso))+" N"
print c_por_piso #

# 
# 
# 

sp.savez('mck.npz',Mmatriz=M,Cmatriz=c,Kmatriz=kmatrix, CAPvector=c_por_piso)
