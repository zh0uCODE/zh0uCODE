using System;

namespace work
{
    class Program
    {
       static void GaussSolve(int n, double[] A, double[] x, double[] b)
       {
        int i, j, k, r;
        double max;

        for (k = 0; k < n - 1; k++)
        {
            max = Math.Abs(A[k * n + k]);
            r = k;
            for (i = k + 1; i < n; i++)
            {
                if (max < Math.Abs(A[i * n + k]))
                {
                    max = Math.Abs(A[i * n + k]);
                    r = i;
                }
            }
            if (r != k)
            {
                for (i = 0; i < n; i++)
                {
                    max = A[k * n + i];
                    A[k * n + i] = A[r * n + i];
                    A[r * n + i] = max;
                }
                max = b[k];
                b[k] = b[r];
                b[r] = max;
            }

            for (i = k + 1; i < n; i++)
            {
                for (j = k + 1; j < n; j++)
                {
                    if (A[k * n + k] == 0)
                    {
                        Console.WriteLine("Error: Division by zero");
                        return;
                    }
                    A[i * n + j] -= A[i * n + k] * A[k * n + j] / A[k * n + k];
                }
                b[i] -= A[i * n + k] * b[k] / A[k * n + k];
            }
        }

        for (i = n - 1; i >= 0; i--)
        {
            x[i] = b[i];
            for (j = i + 1; j < n; j++)
            {
                x[i] -= A[i * n + j] * x[j];
            }
            x[i] /= A[i * n + i];
        }
       }

       static void Polyfit(int n, double[] x, double[] y, int poly_n, double[] p)
       {
        int i, j;
        double[] tempx = new double[n];
        double[] sumxx = new double[poly_n * 2 + 1];
        double[] tempy = new double[n];
        double[] sumxy = new double[poly_n + 1];
        double[] ata = new double[(poly_n + 1) * (poly_n + 1)];

        for (i = 0; i < n; i++)
        {
            tempx[i] = 1;
            tempy[i] = y[i];
        }

        for (i = 0; i < 2 * poly_n + 1; i++)
        {
            sumxx[i] = 0;
            for (j = 0; j < n; j++)
            {
                sumxx[i] += tempx[j];
                tempx[j] *= x[j];
            }
        }

        for (i = 0; i < poly_n + 1; i++)
        {
            sumxy[i] = 0;
            for (j = 0; j < n; j++)
            {
                sumxy[i] += tempy[j];
                tempy[j] *= x[j];
            }
        }

        for (i = 0; i < poly_n + 1; i++)
        {
            for (j = 0; j < poly_n + 1; j++)
            {
                ata[i * (poly_n + 1) + j] = sumxx[i + j];
            }
        }

        GaussSolve(poly_n + 1, ata, p, sumxy);
       }
       double[] LU(double[] A, double[] b, double[] solution)
       {
        const int NVAR = 4;
        int i, j, n;
        double temp = 0;

        double[][] mat = new double[NVAR][];
        for (i = 0; i < NVAR; i++)
        {
            mat[i] = new double[NVAR];
            Array.Clear(mat[i], 0, NVAR);
        }

        double[] x = solution;
        Array.Clear(x, 0, NVAR);

        i = 0;
        for (j = 0; j < NVAR; j++)
        {
            mat[i][j] = A[i * NVAR + j];
        }

        if (A[0] == 0)
        {
            Console.WriteLine("Error!, A[0] must not be 0!");
            return null;
        }

        j = 0;
        for (i = 1; i < NVAR; i++)
        {
            mat[i][j] = A[i * NVAR + j] / A[0];
        }

        for (i = 1; i < NVAR; i++)
        {
            for (j = 1; j < NVAR; j++)
            {
                if (i <= j)
                {
                    temp = 0;
                    for (n = 0; n < i; n++)
                    {
                        temp += mat[n][j] * mat[i][n];
                    }
                    mat[i][j] = A[i * NVAR + j] - temp;
                }
                else
                {
                    temp = 0;
                    for (n = 0; n < j; n++)
                    {
                        temp += mat[n][j] * mat[i][n];
                    }
                    if (mat[j][j] == 0)
                    {
                        Console.WriteLine($"Error! Matrix[{j}][{j}]==0, NO solution.");
                        return null;
                    }
                    mat[i][j] = (A[i * NVAR + j] - temp) / mat[j][j];
                }
            }
        }

        for (i = 0; i < NVAR; i++)
        {
            temp = 0;
            for (j = 0; j < i; j++)
            {
                temp += mat[i][j] * x[j];
            }
            x[i] = b[i] - temp;
        }

        for (i = NVAR - 1; i >= 0; i--)
        {
            temp = 0;
            for (j = NVAR - 1; j > i; j--)
            {
                temp += mat[i][j] * x[j];
            }
            if (mat[i][i] == 0)
            {
                Console.WriteLine($"Error! Matrix[{i}][{i}]==0, NO solution.");
                return null;
            }
            x[i] = (x[i] - temp) / mat[i][i];
        }

        return x;
       }
       static void TestLU(){
        double x = 10.0, y = 20.0, z = 30.0, u = 40.0;
        double[] A = {
            1, 2.0, 3.0, 4.0,
            1.0, 1.0, 3.0, 4.0,
            22 * 1.0, 2 * 2.0, 2 * 3.0, 2 * 4.0,
            1.0, 2.0, 3.0, -4.0
        };
        double[] b = new double[4];
        b[0] = A[0] * x + A[1] * y + A[2] * z + A[3] * u;
        b[1] = A[4] * x + A[5] * y + A[6] * z + A[7] * u;
        b[2] = A[8] * x + A[9] * y + A[10] * z + A[11] * u;
        b[3] = A[12] * x + A[13] * y + A[14] * z + A[15] * u;

        double[] sol = new double[4];
        Array.Clear(sol, 0, sol.Length);

        Program lu = new Program();
        double[] res = lu.LU(A, b, sol);
        if (res == null)
        {
            Console.WriteLine("Error happened, no solution.");
        }
        else
        {
            Console.WriteLine($"{sol[0]} {sol[1]} {sol[2]} {sol[3]}");
        }
       }

       static double PolyCalc(double temp, int order, double[] coeff){
            double sum = coeff[0];
            double x = temp;
            int i;
            for (i=0;i<order;i++) {
                sum += coeff[i+1] * x;
                x *= temp;
            }
            return sum;
       }

        static double Add(double x, double y)
        {
            return x+y;
        }
        static void Main(string[] args)
        {
            Console.WriteLine("Hello World!");
            double z = Add(5.0, 4.3);
            double [] coeff = new double[10];
            double [] x = new double[10];
            double [] y = new double[10];
            double v;
            int num = 0;
            for (v=-2.0; v<2.0; v+=0.5) {
                x[num] = v;
                y[num] = 1.0 + 3.3*v + 4.4*v*v + 5.5* v*v*v;
                num++;
            }
            Polyfit(num, x, y, 5, coeff);
            for (int i=0; i<=5; i++) {
                Console.WriteLine(coeff[i]);
            }
        }
    }
}
