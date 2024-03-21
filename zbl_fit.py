import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

dict = {'H':1, 'He':2, 'Li':3, 'Be':4, 'B':5, 'C':6,'N':7, 'O':8,'F':9,'Ne':10,'Na':11, 
        'Mg':12, 'Al':13, 'Si':14, 'P':15,'S':16,'Cl':17,'Ar':18,'K':19,'Ca':20,'Sc':21,
        'Ti':22, 'V':23, 'Cr':24, 'Mn':25, 'Fe':26, 'Co':27, 'Ni':28, 'Cu':29, 'Zn':30, 
        'Ga':31, 'Ge':32, 'As':33, 'Se':34, 'Br':35, 'Kr':36,'Rb':37, 'Sr':38, 'Y':39, 
        'Zr':40, 'Nb':41, 'Mo':42, 'Tc':43, 'Ru':44, 'Rh':45,'Pd':46, 'Ag':47, 'Cd':48, 
        'Hf':72, 'Ta':73, 'W':74, 'Pt':78, 'Au':79}


arr_ele = input('请输入所有元素，以空格隔开\n')
file = str(input('请输入数据文件名\n'))
r_e = np.loadtxt('%s'%file)
inner = 1
outer = 2

element = [str(n) for n in arr_ele.split()]
n = len(element)

p0 =      [0.3,   50.0,  0.1,   0.3,   0.6,   50.0,  0.5,  50.0]
bounds = ([0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0], \
          [90.0,  900.0, 90.0,  900.0, 90.0,  900.0, 90.0, 900.0])

if n == 1:
    print('error')
elif n == 2:
    e1, e2 = element[0], element[1]
    Z1, Z2 = dict['%s'%e1], dict['%s'%e2]
    def quadratic(r, p1, p2, p3, p4, p5, p6, p7, p8):
        a = 0.46848 / (Z1 ** 0.23 + Z2 ** 0.23) 
        A = 14.399645 * Z1 * Z2
        x = r / a
        return A / r * (p1 * np.exp(-p2 * x) + p3 * np.exp(-p4 * x) + p5 * np.exp(-p6 * x) + p7 * np.exp(-p8 * x)) 
    r = r_e[:,0]
    y = r_e[:,1]
    popt, pcov = curve_fit(quadratic, r, y, p0 = p0, bounds = bounds, maxfev=1000000)
    zbl = np.array([inner, outer] + list(popt))
    np.savetxt('zbl.in.txt', zbl, fmt='%f')
    np.savetxt('%s-%s.txt'%(e1,e2), np.column_stack((np.linspace(0.4, 2, 100),quadratic(np.linspace(0.4, 2, 100), *popt))), delimiter=' ', fmt='%.4f %.4f')
    plt.scatter(r, y, label = 'DFT')
    plt.plot(np.linspace(np.min(r), np.max(r), 100), quadratic(np.linspace(np.min(r), np.max(r), 100), *popt), color='r', label=' fittedZBL')
    plt.xlabel('r (ANG.)')
    plt.ylabel('E (eV)')
    plt.legend()
    plt.savefig('%s-%s.png'%(e1, e2))
    plt.close()
else:
    num = 0
    zbl = []
    for i in range(n):
        for j in range(n):
            if j < i:
                pass
            else:
                e1, e2 = element[i], element[j]
                Z1, Z2 = dict['%s'%e1], dict['%s'%e2]
                def quadratic(r, p1, p2, p3, p4, p5, p6, p7, p8):
                    a = 0.46848 / (Z1 ** 0.23 + Z2 ** 0.23) 
                    A = 14.399645 * Z1 * Z2
                    x = r / a
                    return A / r * (p1 * np.exp(-p2 * x) + p3 * np.exp(-p4 * x) + p5 * np.exp(-p6 * x) + p7 * np.exp(-p8 * x)) 
                r = r_e[:,2*num]
                y = r_e[:,2*num + 1]
                num += 1
                popt, pcov = curve_fit(quadratic, r, y, p0 = p0, bounds = bounds, maxfev=1000000)
                zbl.append([inner, outer] + list(popt))
                np.savetxt('%s-%s.txt'%(e1,e2), np.column_stack((r,y,quadratic(r, *popt))), delimiter=' ', fmt='%.4f %.4f %.4f')
                plt.scatter(r, y, label = 'DFT')
                plt.plot(np.linspace(np.min(r), np.max(r), 100), quadratic(np.linspace(np.min(r), np.max(r), 100), *popt), color='r', label=' fittedZBL')
                plt.xlabel('r (ANG.)')
                plt.ylabel('E (eV)')
                plt.legend()
                plt.savefig('%s-%s.png'%(e1, e2))
                plt.close()
    np.savetxt('zbl.in.txt', np.array(zbl), fmt='%4d%4d%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f')     

