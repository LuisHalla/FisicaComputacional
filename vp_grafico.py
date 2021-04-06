from matplotlib import pyplot as plt
import numpy as np
import os

os.system("./dominante.exe")

plt.style.use('dark_background')

data = np.loadtxt('vp_datos.txt')

x = data[:,0]
y = data[:,1]
z = data[:,2]

fig = plt.figure()

a = fig.add_subplot(121)
a.plot(x,y, color='yellow')

a.set_title("Power Method", fontsize=30)
a.set_xlabel("Iteraciones", fontsize=25)
a.set_ylabel('$\lambda_{PM}$', fontsize=25)

b = fig.add_subplot(122)
b.plot(x,z, color='yellow')

b.set_title("Comparacion con calculo exacto", fontsize=30)
b.set_xlabel("Iteraciones", fontsize=25)
b.set_ylabel('$|\lambda_{exacto} - \lambda_{PM}$|', fontsize=25)



plt.show()
