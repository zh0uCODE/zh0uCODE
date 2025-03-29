import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import csv
#temperature values:
vc = []
CAP1 = []
CAP2 = []
with open('CSV File For Proj - Sheet3.csv','r') as csv_file:
    csv_reader = csv.reader(csv_file)
    for line in csv_reader:
        vc.append(line[0])
        CAP1.append(line[1]) #this is good
        CAP2.append(line[3])
    for singleQuote in vc[:]: #remove single quotes (voltage)
        if singleQuote == '':
            vc.remove(singleQuote)
    for singleQuote in CAP1[:]: #remove single quotes (CAP1)
        if singleQuote == '':
            CAP1.remove(singleQuote)
    for singleQuote in CAP2[:]: #remove single quotes (CAP2)
        if singleQuote == '':
            CAP2.remove(singleQuote)
vc.remove('vc/V')
CAP1.remove('BCD-CAP7')
CAP1.remove('CAP7#1/pf')
CAP2.remove('CAP7#2/pf')
for i in range(len(vc)):
    vc[i] = float(vc[i])
for i in range(len(CAP1)):
    CAP1[i] = float(CAP1[i])
for i in range(len(CAP2)):
    CAP2[i] = float(CAP2[i])



print(vc)
print(CAP1)
print(CAP2)
plt.xlabel("Voltage")
plt.ylabel("CAP")
plt.title("CAP vs Voltage")
plt.plot(vc, CAP1)
plt.plot(vc, CAP2)
c = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] #21 order good!
def poly_calc(temp, order, coeff):
    sum = coeff[0]
    x = temp

    for i in range(order):
        sum += coeff[i + 1] * x
        x *= temp

    return sum #good!!!
def polyfit(n, x, y, poly_n, p):
    tempx = [1.0] * n
    sumxx = [0.0] * (2 * poly_n + 1)
    tempy = y[:]
    sumxy = [0.0] * (poly_n + 1)
    ata = [0.0] * ((poly_n + 1) * (poly_n + 1))

    for i in range(n):
        tempx[i] = 1
        tempy[i] = y[i]

    for i in range(2 * poly_n + 1):
        sumxx[i] = 0.0
        for j in range(n):
            sumxx[i] += tempx[j]
            tempx[j] *= x[j]

    for i in range(poly_n + 1):
        sumxy[i] = 0.0
        for j in range(n):
            sumxy[i] += tempy[j]
            tempy[j] *= x[j]

    for i in range(poly_n + 1):
        for j in range(poly_n + 1):
            ata[i * (poly_n + 1) + j] = sumxx[i + j]

    gauss_solve(poly_n + 1, ata, p, sumxy)
def gauss_solve(n, a, x, b):
    index = list(range(n))
    c = [0.0] * n

    for i in range(n):
        index[i] = i

    for i in range(n):
        c[i] = 0
        for j in range(n):
            if abs(a[i * n + j]) > c[i]:
                c[i] = abs(a[i * n + j])

    for j in range(n):
        max_val = 0
        for i in range(j, n):
            temp = abs(a[index[i] * n + j]) / c[index[i]]
            if temp > max_val:
                max_val = temp
                k = i

        temp = index[j]
        index[j] = index[k]
        index[k] = temp

        for i in range(j + 1, n):
            temp = a[index[i] * n + j] / a[index[j] * n + j]
            a[index[i] * n + j] = temp
            for k in range(j + 1, n):
                a[index[i] * n + k] -= temp * a[index[j] * n + k]

    for i in range(n):
        x[i] = b[index[i]]

    for i in range(1, n):
        for j in range(i):
            x[i] -= a[index[i] * n + j] * x[j]

    x[n - 1] = x[n - 1] / a[index[n - 1] * n + (n - 1)]
    for i in range(n - 2, -1, -1):
        for j in range(i + 1, n):
            x[i] -= a[index[i] * n + j] * x[j]
        x[i] = x[i] / a[index[i] * n + i]
print(poly_calc(len(c), len(c)-1, c))
a = len(vc)
polyfit(a, vc, CAP1, len(c)-1, c)
for i in range(len(c)):
    print(c[i])
print(c)
order = len(c) - 1
print(order)
x=0
for i in range(len(c)):
    print("+ (v(p,n)^"+str(x)+")*("+str(c[i])+") +")
    x=x+1
voltageValues = np.linspace(-3.3,3.3,500)
capacitanceValues = [poly_calc(v, order, c) for v in voltageValues]
plt.plot(voltageValues, capacitanceValues) #plots a line!
plt.xlim(-3.3,3.3)
plt.ylim(40,200)
plt.show()




