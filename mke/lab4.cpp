#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

/* Константы */
// Ограничния
const double X0 = 0.0;
const double XN = 90.0;

// Граничные условия
const double BC0 = 0.0;
const double BCN = 2.0;

const int N_CUBIC = 20;

/* Решение СЛАУ с трехдиагональной матрицей методом Томаса */
void solve_thomas(double* a, double* b, double* c, double* d, int n) 
{
    /*
    // n is the number of unknowns

    |b0 c0 0 ||x0| |d0|
    |a1 b1 c1||x1|=|d1|
    |0  a2 b2||x2| |d2|

    1st iteration: b0x0 + c0x1 = d0 -> x0 + (c0/b0)x1 = d0/b0 ->

        x0 + g0x1 = r0               where g0 = c0/b0        , r0 = d0/b0

    2nd iteration:     | a1x0 + b1x1   + c1x2 = d1
        from 1st it.: -| a1x0 + a1g0x1        = a1r0
                    -----------------------------
                          (b1 - a1g0)x1 + c1x2 = d1 - a1r0

        x1 + g1x2 = r1               where g1=c1/(b1 - a1g0) , r1 = (d1 - a1r0)/(b1 - a1g0)

    3rd iteration:      | a2x1 + b2x2   = d2
        from 2nd it. : -| a2x1 + a2g1x2 = a2r2
                       -----------------------
                       (b2 - a2g1)x2 = d2 - a2r2
        x2 = r2                      where                     r2 = (d2 - a2r2)/(b2 - a2g1)
    Finally we have a triangular matrix:
    |1  g0 0 ||x0| |r0|
    |0  1  g1||x1|=|r1|
    |0  0  1 ||x2| |r2|

    Condition: ||bi|| > ||ai|| + ||ci||

    in this version the c matrix reused instead of g
    and             the d matrix reused instead of r and x matrices to report results
    Written by Keivan Moradi, 2014
    */
    n--; // since we start from x0 (not x1)
    c[0] /= b[0];
    d[0] /= b[0];

    for (int i = 1; i < n; i++) {
        c[i] /= b[i] - a[i]*c[i-1];
        d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
    }

    d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);

    for (int i = n; i-- > 0;) {
        d[i] -= c[i]*d[i+1];
    }
}

/* Аналитическое решение */
// Решение ДУ
double u_analytical(double x)
{
    const double C1 = -(2.0*(2.0*exp(270.0*sqrt(2.0/71.0)) - 5.0)) / (3.0*(exp(540.0*sqrt(2.0/71.0)) - 1.0));
    const double C2 = -(2.0*exp(270.0*sqrt(2.0/71.0))*(5.0*exp(270.0*sqrt(2.0/71.0)) - 2.0)) / (3.0*(exp(540.0*sqrt(2.0/71.0)) - 1.0));

    return (C1*exp(3.0*sqrt(2.0/71.0)*x) + C2*exp(-3.0*sqrt(2.0/71.0)*x) + 10.0/3.0);
}

int solve_analytical(int n, double* &x, double* &u)
{
    u = new double[n];

    for (int i = 0; i < n; i++)
    {
        u[i] = u_analytical(x[i]);
    } 

    return n;
}

/* Задание ГУ */
void border_condition(int len, double value, bool left, int type, double *a, double *b, double *c, double *d, double der_b[2])
{
    if (type == 1)
    {
        int i;
        if (left)
        {
            i = 0;
        }
        else
        {
            i = len - 1;
        }

        a[i] = 0.0;
        b[i] = 1.0;
        c[i] = 0.0;

        d[i] = value;
    }
    else if (type == 2)
    {
        int i;
        double coef;
        if (left)
        {
            i = 0;
            coef = der_b[0];
        }
        else
        {
            i = len - 1;
            coef = der_b[1];
        }

        d[i] += coef * value;
    }
}

/* Решение линейными КЭ */
// ММ конечного элемента
void finite_element_linear(double a[2][2], double b[2], double der_b[2], double l)
{
    a[0][0] = 71.0/l + 18.0/3.0 * l; 
    a[0][1] = -71.0/l + 18.0/6.0 * l;

    a[1][0] = -71.0/l + 18.0/6.0 * l; 
    a[1][1] = 71.0/l + 18.0/3.0 * l;

    b[0] = 60.0 * l/2.0;
    b[1] = 60.0 * l/2.0;

    der_b[0] = -71.0;
    der_b[1] = 71.0;
}

int solve_linear(int n, double* &x, double* &u)
{
    int len = n + 1;

    double l = (XN - X0) / n;

    // ММ конечного элемента (без производных)
    double fe_a[2][2], fe_b[2], fe_der_b[2];
    finite_element_linear(fe_a, fe_b, fe_der_b, l);

    // Ансамблирование
    double *a, *b, *c, *d; // Матрица abc и вектор d
    a = new double[len];
    b = new double[len];
    c = new double[len];
    d = new double[len];
    for (int i = 0; i < len; i++)
    {
        a[i] = b[i] = c[i] = d[i] = 0.0;
    }

    for (int i = 0; i < n; i++)
    {
        b[i] += fe_a[0][0]; c[i] += fe_a[0][1];
        a[i + 1] += fe_a[1][0]; b[i + 1] += fe_a[1][1];

        d[i] += fe_b[0];
        d[i + 1] += fe_b[1];
    }
    // ГУ
    border_condition(len, BC0, true, 1, a, b, c, d, fe_der_b);
    border_condition(len, BCN, false, 1, a, b, c, d, fe_der_b);

    // Решение СЛАУ
    solve_thomas(a, b, c, d, len);

    // Сохранение результатов
    x = new double[len];
    u = new double[len];

    double curr_x = X0;
    for (int i = 0; i < len; i++)
    {
        x[i] = curr_x;
        u[i] = d[i];

        curr_x += l;
    }

    // Освобождение памяти
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;

    return len;
}

/* Решение кубическим КЭ */
// ММ конечного элемента (без производных)
void finite_element_cubic(double a[4][4], double b[4], double der_b[2], double l)
{
    // Полный КЭ
    a[0][0] = 71.0 * 37.0/(30.0*l) - 18.0 * (-8.0*l)/35.0; 
    a[0][1] = 71.0 * (-63.0)/(40.0*l) - 18.0 * (-99.0*l)/560.0; 
    a[0][2] = 71.0 * 9.0/(20.0*l) - 18.0 * 9.0*l/140.0; 
    a[0][3] = 71.0 * (-13.0)/(120.0*l) - 18.0 * (-19.0*l)/560.0;

    a[1][0] = a[0][1]; 
    a[1][1] = 71.0 * 18.0/(5.0*l) - 18.0 * (-81.0*l)/70.0; 
    a[1][2] = 71.0 * (-99.0)/(40.0*l) - 18.0 * 81.0*l/560.0; 
    a[1][3] = a[0][2]; 

    a[2][0] = a[0][2]; 
    a[2][1] = a[1][2]; 
    a[2][2] = a[1][1]; 
    a[2][3] = a[0][1];

    a[3][0] = a[0][3]; 
    a[3][1] = a[1][3]; 
    a[3][2] = a[2][3]; 
    a[3][3] = a[0][0];

    b[0] = 60.0 * 3.0*l/8.0;
    b[1] = 60.0 * 9.0*l/8.0;
    b[2] = b[1];
    b[3] = b[0];

    der_b[0] = -71.0;
    der_b[1] = 71.0;

    // Упрощение методом Гаусса
    int index[4] = {1, 2, 0, 3};
    for (int* k = index; k < index + 2; k++) // 2 итерации
    {
        for (int* i = k + 1; i < index + 4; i++)
        {
            double coef = a[*i][*k] / a[*k][*k];
            a[*i][*k] = 0.0;
            for (int* j = k + 1; j < index + 4; j++)
            {
                a[*i][*j] -= coef * a[*k][*j];
            }
            b[*i] -= coef * b[*k];
        }
    }
}

void restore_cubic(double *x, double *u, double a[4][4], double b[4], double l)
{
    u[2] = (b[2] - a[2][0]*u[0] - a[2][3]*u[3]) / a[2][2];
    u[1] = (b[1] - a[1][0]*u[0] - a[1][2]*u[2] - a[1][3]*u[3]) / a[1][1];

    x[1] = x[0] + l;
    x[2] = x[0] + 2.0*l; 
}

int solve_cubic(int n, double* &x, double* &u)
{
    int size = n + 1;
    int len = 3*n + 1;

    double l = ((XN - X0) / n) / 3.0;

    // ММ конечного элемента (без производных)
    double fe_a[4][4], fe_b[4], fe_der_b[2];
    finite_element_cubic(fe_a, fe_b, fe_der_b, l);

    // Ассамблирование
    double *a, *b, *c, *d; // Матрица abc и вектор d
    a = new double[size];
    b = new double[size];
    c = new double[size];
    d = new double[size];
    for (int i = 0; i < size; i++)
    {
        a[i] = b[i] = c[i] = d[i] = 0.0;
    }

    // Между границами
    for (int i = 0; i < n; i++)
    {
        b[i] += fe_a[0][0]; c[i] += fe_a[0][3];
        a[i + 1] += fe_a[3][0]; b[i + 1] += fe_a[3][3];

        d[i] += fe_b[0];
        d[i + 1] += fe_b[3];
    }
    // ГУ
    border_condition(size, BC0, true, 1, a, b, c, d, fe_der_b);
    border_condition(size, BCN, false, 1, a, b, c, d, fe_der_b);

    // Решение СЛАУ
    solve_thomas(a, b, c, d, size);

    // Сохранение результатов
    x = new double[len];
    u = new double[len];

    double curr_x = X0;

    x[0] = curr_x;
    u[0] = d[0];

    curr_x += 3.0*l;

    for (int i = 1; i < size; i++)
    {
        x[3*i] = curr_x;
        u[3*i] = d[i];

        restore_cubic(x + 3*(i - 1), u + 3*(i - 1), fe_a, fe_b, l);

        curr_x += 3.0*l;
    }

    // Освобождение памяти
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;

    return len;
}

/* Вывод данных в файл */
void output(const char *filename, int n, double *x, double *u)
{
    ofstream file(filename);

    for (int i = 0; i < n; i++)
    {
        file << x[i] << ' ' << u[i] << '\n';
    }
}

void output_error(const char *filename, int n, double *x, double *u_analitical, double *u_fe)
{
    ofstream file(filename);

    for (int i = 0; i < n; i++)
    {
        file << x[i] << ' ' << u_analitical[i] << ' ' << u_fe[i] << ' ' << abs(u_fe[i] - u_analitical[i]) << '\n';
    }
}

/* Нахождение погрешности */
double max_error(int n, double *x, double *u)
{
    double max = 0.0;
    for (int i = 0; i < n; i++)
    {
        double curr_error = abs(u[i] - u_analytical(x[i]));
        if (curr_error > max)
        {
            max = curr_error;
        }
    }

    return max;
}

/* main */
int main()
{
    // Кубические КЭ
    int len_cubic;
    double *x_cubic, *u_cubic;
    double *u_analytical_for_cubic;
    double error_cubic;
    len_cubic = solve_cubic(N_CUBIC, x_cubic, u_cubic);
    solve_analytical(len_cubic, x_cubic, u_analytical_for_cubic);
    error_cubic = max_error(len_cubic, x_cubic, u_cubic);

    // Линейные КЭ
    int len_linear;
    double *x_linear, *u_linear;
    double *u_analytical_for_linear;
    double error_linear;
    double delta_error_prev, delta_error_curr = 1000.0;
    
    int n = 3 * N_CUBIC;
    do
    {
        n++;
        len_linear = solve_linear(n, x_linear, u_linear);

        error_linear = max_error(len_linear, x_linear, u_linear);
        delta_error_prev = delta_error_curr;
        delta_error_curr = abs(error_cubic - error_linear);

        delete[] x_linear;
        delete[] u_linear;

    } while (delta_error_curr < delta_error_prev);

    n--;
    len_linear = solve_linear(n, x_linear, u_linear);
    solve_analytical(len_linear, x_linear, u_analytical_for_linear);
    error_linear = max_error(len_linear, x_linear, u_linear);

    // Вывод результатов
    cout << "Погрешность для " << N_CUBIC << " кубических КЭ: " << error_cubic << '\n';
    output("analytical_for_cubic.txt", len_cubic, x_cubic, u_analytical_for_cubic);
    output("cubic.txt", len_cubic, x_cubic, u_cubic);
    output_error("cubic_error.txt", len_cubic, x_cubic, u_analytical_for_cubic, u_cubic);

    cout << "Погрешность для " << n << " линейных КЭ: " << error_linear << '\n';
    output("analytical_for_linear.txt", len_linear, x_linear, u_analytical_for_linear);
    output("linear.txt", len_linear, x_linear, u_linear);
    output_error("linear_error.txt", len_linear, x_linear, u_analytical_for_linear, u_linear);
    
    // Освобождение памяти
    delete[] x_cubic;
    delete[] u_cubic;
    delete[] u_analytical_for_cubic;

    delete[] x_linear;
    delete[] u_linear;
    delete[] u_analytical_for_linear;

    return 0;
}
