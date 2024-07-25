import matplotlib.pyplot as plt
import numpy as np


def poly_calc(temp, order, coeff):
    sum = coeff[0]
    x = temp

    for i in range(order):
        sum += coeff[i + 1] * x
        x *= temp

    return sum #good!!!
arr = [10, 2, 3, 4] #graph good # Example coefficients
order = len(arr) - 1

# Generate a range of temperature values
temp_values = np.linspace(-10, 10, 400)

# Evaluate the polynomial at each temperature value
poly_values = [poly_calc(temp, order, arr) for temp in temp_values]

# Plot the results
plt.plot(temp_values, poly_values, label='Polynomial')
plt.xlabel('Temperature')
plt.ylabel('Polynomial Value')
plt.title('Polynomial Evaluation')
plt.legend()
plt.grid(True)
plt.show()
print(poly_calc(5, 2, arr))
def LU(A, b, solution):
    NVAR = 4
    mat = [[0.0] * NVAR for _ in range(NVAR)]
    x = solution[:]  # Make a copy of the solution list
    for i in range(NVAR):
        x[i] = 0.0

    # Initialize mat matrix with A values
    i = 0
    for j in range(NVAR):
        mat[i][j] = A[i * NVAR + j]

    if A[0] == 0:
        print("Error! A[0] must not be 0!")
        return None

    j = 0
    for i in range(1, NVAR):
        mat[i][j] = A[i * NVAR + j] / A[0]

    for i in range(1, NVAR):
        for j in range(1, NVAR):
            if i <= j:
                temp = 0.0
                for n in range(i):
                    temp += mat[n][j] * mat[i][n]
                mat[i][j] = A[i * NVAR + j] - temp
            else:
                temp = 0.0
                for n in range(j):
                    temp += mat[n][j] * mat[i][n]
                if mat[j][j] == 0:
                    print(f"Error! Matrix[{j}][{j}]==0, NO solution.")
                    return None
                mat[i][j] = (A[i * NVAR + j] - temp) / mat[j][j]

    for i in range(NVAR):
        temp = 0.0
        for j in range(i):
            temp += mat[i][j] * x[j]
        x[i] = b[i] - temp

    for i in range(NVAR - 1, -1, -1):
        temp = 0.0
        for j in range(NVAR - 1, i, -1):
            temp += mat[i][j] * x[j]
        if mat[i][i] == 0:
            print(f"Error! Matrix[{i}][{i}]==0, NO solution.")
            return None
        x[i] = (x[i] - temp) / mat[i][i]

    return x

def TestLU():
    x = 25.0
    y = 35.0
    z = 45.0
    u = 55.0
    A = [
        1, 2.0, 3.0, 4.0,
        1.0, 1.0, 3.0, 4.0,
        22 * 1.0, 2 * 2.0, 2 * 3.0, 2 * 4.0,
        1.0, 2.0, 3.0, -4.0
    ]
    b = [0.0] * 4
    b[0] = A[0] * x + A[1] * y + A[2] * z + A[3] * u
    b[1] = A[4] * x + A[5] * y + A[6] * z + A[7] * u
    b[2] = A[8] * x + A[9] * y + A[10] * z + A[11] * u
    b[3] = A[12] * x + A[13] * y + A[14] * z + A[15] * u

    sol = [0.0] * 4

    res = LU(A, b, sol)
    if res is None:
        print("Error happened, no solution.")
    else:
        print("{:.1f}  {:.1f}  {:.1f}  {:.1f}".format(res[0], res[1], res[2], res[3]))

# Run the test
TestLU() #good!!!
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
coeff = [0.0] * 10
x = [0.0] * 10
y = [0.0] * 10
num = 0

for v in range(-4, 4, 1):
    v = v / 2.0
    x[num] = v
    y[num] = 2.2 + 5.5 * v + 6.6 * v * v + 9.9 * v * v * v
    num += 1

polyfit(num, x, y, 5, coeff)
for i in range(6):
    print(coeff[i])
#good! Error bound very small!

