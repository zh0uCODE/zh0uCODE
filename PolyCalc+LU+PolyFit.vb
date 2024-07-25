Imports System 'System is a Namespace  
Module Hello_Program

    Sub Main(args As String())

        Console.WriteLine("Hello, Welcome to the world of VB.NET")
        Console.WriteLine("Press any key to continue...")
        Dim coeff(9) As Double
        Dim x(9) As Double
        Dim y(9) As Double
        Dim v As Double
        Dim num As Integer = 0
        Dim c() As Double = {1.0, 2.0, 3.0, 4.0}
        Dim a As Double = PolyCalc(5.0, 2, c)
        Console.WriteLine(a)
    End Sub

    Public Function PolyCalc(ByVal temp As Double, ByVal order As Integer, ByVal coeff As Double()) As Double
        Dim sum As Double = coeff(0)
        Dim x As Double = temp
        Dim i As Integer

        For i = 0 To order - 1
            sum += coeff(i + 1) * x
            x *= temp
        Next

        Return sum
    End Function
    Sub TestLU()
        Dim x As Double = 25.0
        Dim y As Double = 35.0
        Dim z As Double = 45.0
        Dim u As Double = 55.0
        Dim A() As Double = {
        1, 2.0, 3.0, 4.0,
        1.0, 1.0, 3.0, 4.0,
        22 * 1.0, 2 * 2.0, 2 * 3.0, 2 * 4.0,
        1.0, 2.0, 3.0, -4.0
        }
        Dim b(3) As Double
        b(0) = A(0) * x + A(1) * y + A(2) * z + A(3) * u
        b(1) = A(4) * x + A(5) * y + A(6) * z + A(7) * u
        b(2) = A(8) * x + A(9) * y + A(10) * z + A(11) * u
        b(3) = A(12) * x + A(13) * y + A(14) * z + A(15) * u

        Dim sol(3) As Double
        Array.Clear(sol, 0, sol.Length)

        Dim res() As Double = LU(A, b, sol)
        If res Is Nothing Then
            Console.WriteLine("Error happened, no solution.")
        Else
            'Console.WriteLine("{sol(0)} {sol(1)} {sol(2)} {sol(3)}")
            Console.WriteLine("{0:F}  {1:F}  {2:F}  {3:F}", sol(0), sol(1), sol(2), sol(3))
        End If
    End Sub
    Function LU(ByVal A() As Double, ByVal b() As Double, ByVal solution() As Double) As Double()
        Const NVAR As Integer = 4
        Dim i, j, n As Integer
        Dim temp As Double = 0

        Dim mat(NVAR - 1)() As Double
        For i = 0 To NVAR - 1
            mat(i) = New Double(NVAR - 1) {}
            Array.Clear(mat(i), 0, NVAR)
        Next

        Dim x() As Double = solution
        Array.Clear(x, 0, NVAR)

        i = 0
        For j = 0 To NVAR - 1
            mat(i)(j) = A(i * NVAR + j)
        Next

        If A(0) = 0 Then
            Console.WriteLine("Error!, A[0] must not be 0!")
            Return Nothing
        End If

        j = 0
        For i = 1 To NVAR - 1
            mat(i)(j) = A(i * NVAR + j) / A(0)
        Next

        For i = 1 To NVAR - 1
            For j = 1 To NVAR - 1
                If i <= j Then
                    temp = 0
                    For n = 0 To i - 1
                        temp += mat(n)(j) * mat(i)(n)
                    Next
                    mat(i)(j) = A(i * NVAR + j) - temp
                Else
                    temp = 0
                    For n = 0 To j - 1
                        temp += mat(n)(j) * mat(i)(n)
                    Next
                    If mat(j)(j) = 0 Then
                        Console.WriteLine("Error! Matrix[{j}][{j}]==0, NO solution.")
                        Return Nothing
                    End If
                    mat(i)(j) = (A(i * NVAR + j) - temp) / mat(j)(j)
                End If
            Next
        Next

        For i = 0 To NVAR - 1
            temp = 0
            For j = 0 To i - 1
                temp += mat(i)(j) * x(j)
            Next
            x(i) = b(i) - temp
        Next

        For i = NVAR - 1 To 0 Step -1
            temp = 0
            For j = NVAR - 1 To i + 1 Step -1
                temp += mat(i)(j) * x(j)
            Next
            If mat(i)(i) = 0 Then
                Console.WriteLine("Error! Matrix[{i}][{i}]==0, NO solution.")
                Return Nothing
            End If
            x(i) = (x(i) - temp) / mat(i)(i)
        Next

        Return x
    End Function
    Sub Polyfit(ByVal n As Integer, ByVal x() As Double, ByVal y() As Double, ByVal poly_n As Integer, ByVal p() As Double)
        Dim i, j As Integer
        Dim tempx(n - 1) As Double
        Dim sumxx(poly_n * 2) As Double
        Dim tempy(n - 1) As Double
        Dim sumxy(poly_n) As Double
        Dim ata((poly_n + 1) * (poly_n + 1) - 1) As Double

        For i = 0 To n - 1
            tempx(i) = 1
            tempy(i) = y(i)
        Next

        For i = 0 To 2 * poly_n
            sumxx(i) = 0
            For j = 0 To n - 1
                sumxx(i) += tempx(j)
                tempx(j) *= x(j)
            Next
        Next

        For i = 0 To poly_n
            sumxy(i) = 0
            For j = 0 To n - 1
                sumxy(i) += tempy(j)
                tempy(j) *= x(j)
            Next
        Next

        For i = 0 To poly_n
            For j = 0 To poly_n
                ata(i * (poly_n + 1) + j) = sumxx(i + j)
            Next
        Next

        GaussSolve(poly_n + 1, ata, p, sumxy)
    End Sub

    Sub GaussSolve(ByVal n As Integer, ByVal a() As Double, ByVal x() As Double, ByVal b() As Double)
        Dim i, j, k As Integer
        Dim max As Double
        Dim temp As Double
        Dim index(n - 1) As Integer
        Dim c(n - 1) As Double

        For i = 0 To n - 1
            index(i) = i
        Next

        For i = 0 To n - 1
            c(i) = 0
            For j = 0 To n - 1
                If Math.Abs(a(i * n + j)) > c(i) Then
                    c(i) = Math.Abs(a(i * n + j))
                End If
            Next
        Next

        For j = 0 To n - 1
            max = 0
            For i = j To n - 1
                temp = Math.Abs(a(index(i) * n + j)) / c(index(i))
                If temp > max Then
                    max = temp
                    k = i
                End If
            Next

            temp = index(j)
            index(j) = index(k)
            index(k) = temp

            For i = j + 1 To n - 1
                temp = a(index(i) * n + j) / a(index(j) * n + j)
                a(index(i) * n + j) = temp
                For k = j + 1 To n - 1
                    a(index(i) * n + k) -= temp * a(index(j) * n + k)
                Next
            Next
        Next

        For i = 0 To n - 1
            x(i) = b(index(i))
        Next

        For i = 1 To n - 1
            For j = 0 To i - 1
                x(i) -= a(index(i) * n + j) * x(j)
            Next
        Next

        x(n - 1) = x(n - 1) / a(index(n - 1) * n + (n - 1))
        For i = n - 2 To 0 Step -1
            For j = i + 1 To n - 1
                x(i) -= a(index(i) * n + j) * x(j)
            Next
            x(i) = x(i) / a(index(i) * n + i)
        Next
    End Sub

    End Module


