#include <iostream>
#include "sparse_lu.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>

// Хранимые моменты времени
enum TimePos
{
    CURR_TIME,
    PREV_TIME,

    TIME_POS_COUNT_REACTIVE, // Количество для значений в реактивных элементах

    PREV_PREV_TIME = TIME_POS_COUNT_REACTIVE,

    TIME_POS_COUNT_BASIC // Количество для базисных хначений
};

/* Параметры моделирования */
const double TIME = 1e-3;
const double DELTA_TIME_MIN = 1e-9;
const double DELTA_TIME_PRINT = 0.0;

const int NEWTON_MAX_ITER = 100;
const double NEWTON_EPSILON = 1e-4;

const double ERROR_DELTA_1 = 1e-4;
const double ERROR_DELTA_2 = 1e-2;

// Константы модели
const double L1 = 0.001;
const double L2 = 0.001;
const double C1 = 1e-6;
const double R1 = 1000.0;
const double R2 = 1000.0;
// Источник
const double A = 10.0;
const double P = 1e-4;
// Диод
const double I_T = 1e-12;
const double C_B = 2e-12;
const double MF_T = 0.026;
const double TAU = 5e-9;
const double R_U = 1000000.0;
const double R_B = 20.0;

// Базисные переменные
enum BasicValueName
{
    // Переменные интегрирования

    INT_VALUE_NAMES_COUNT, // Количество переменных интегрирования

    // Остальные переменные
    I_L1 = INT_VALUE_NAMES_COUNT,
    I_L2,
    FI1,
    FI2,
    FI3,
    FI4,
    I_E,

    BASIC_VALUE_NAMES_COUNT // Количество переменных
};

// Базисные переменные на вывод
const BasicValueName basicValuesToPrint[] = { FI2, FI4, I_L1, I_L2 };
const int BASIC_VALUES_TO_PRINT_COUNT = sizeof(basicValuesToPrint) / sizeof(BasicValueName);

// Переменные реактивных элементов (не входящие в базис)
enum ReactiveValueName
{
    U_C1,
    U_C,

    REACTIVE_VALUE_NAMES_COUNT // Количество переменных
};

// Расчет переменных реактивных элементов
double *getReactiveValues(
    double *reactiveValuesCurr, const double *reactiveValuesPrev,
    const double *basicValuesCurr, const double *basicValuesPrev,
    double delta_time, double currTime
)
{
    reactiveValuesCurr[U_C1] = basicValuesCurr[FI2] - basicValuesCurr[FI1];
    reactiveValuesCurr[U_C] = basicValuesCurr[FI3] - basicValuesCurr[FI2];

    return reactiveValuesCurr;
}

// Переменных реактивных элементов на вывод
const ReactiveValueName reactiveValuesToPrint[] = {  };
const int REACTIVE_VALUES_TO_PRINT_COUNT = sizeof(reactiveValuesToPrint) / sizeof(ReactiveValueName);

// Установка начальных значений
void setInitialValues(double* const *basicValues, double* const *reactiveValues, double *timePos)
{
    for (int i = 0; i < TIME_POS_COUNT_BASIC; i++)
    {
        // Базисные переменные
        for (int value = 0; value < BASIC_VALUE_NAMES_COUNT; value++)
        {
            basicValues[i][value] = 0.0;
        }
        // Моменты времени
        timePos[i] = 0.0;
    }
    for (int i = 0; i < TIME_POS_COUNT_REACTIVE; i++)
    {
        // Переменные в реактивных элементах
        for (int value = 0; value < REACTIVE_VALUE_NAMES_COUNT; value++)
        {
            reactiveValues[i][value] = 0.0;
        }
    }
}

// Получение математической модели для метода Ньютона
void getNewtonModel(
    Sparse &jacobi, double *minus_f, 
    const double *basicValuesCurr, const double *basicValuesPrev, 
    const double *reactiveValuesPrev, 
    double delta_time, double currTime
)
{
    // Вычисляемые значения
    const double _E = A*std::sin(2.0*M_PI*currTime / P);
    const double _C = C_B + TAU/MF_T * I_T*std::exp((basicValuesCurr[FI3] - basicValuesCurr[FI2]) / MF_T);
    const double _I = I_T*(std::exp((basicValuesCurr[FI3] - basicValuesCurr[FI2]) / MF_T) - 1.0);

    const double _A = 1.0/delta_time * (_C + ((basicValuesCurr[FI3] - basicValuesCurr[FI2]) - reactiveValuesPrev[U_C]) * TAU / (MF_T*MF_T) * I_T*std::exp((basicValuesCurr[FI3] - basicValuesCurr[FI2]) / MF_T));
    const double _B = I_T/MF_T * std::exp((basicValuesCurr[FI3] - basicValuesCurr[FI2]) / MF_T);

    const double _I_R1 = (basicValuesCurr[FI2] - basicValuesCurr[FI1]) / R1;
    const double _I_C1 = C1/delta_time * ((basicValuesCurr[FI2] - basicValuesCurr[FI1]) - reactiveValuesPrev[U_C1]);
    const double _I_C = _C/delta_time * ((basicValuesCurr[FI3] - basicValuesCurr[FI2]) - reactiveValuesPrev[U_C]);
    const double _I_R_U = (basicValuesCurr[FI3] - basicValuesCurr[FI2]) / R_U;
    const double _I_R_B = (basicValuesCurr[FI4] - basicValuesCurr[FI3]) / R_B;
    const double _I_R2 = basicValuesCurr[FI4] / R2;

    // Формирование матрицы Якоби, хранящейся в формате CSR
    static int staticLineEnds[] = { 3, 5, 9, 13, 16, 19, 20 };
    static int staticIndices[] = {
        0, 2, 3,
        1, 5,
        0, 2, 3, 6,
        0, 2, 3, 4,
        3, 4, 5,
        1, 4, 5,
        2

    };
    static double staticValues[] = {
        NAN, -1.0, 1.0,
        NAN, 1.0,
        -1.0, NAN, NAN, 1.0,
        1.0, NAN, NAN, NAN,
        NAN, NAN, -1.0/R_B,
        1.0, -1.0/R_B, (1.0/R_B + 1.0/R2),
        1.0
    };

    staticValues[0] = -L1/delta_time;

    staticValues[staticLineEnds[0] + 0] = -L2/delta_time;

    staticValues[staticLineEnds[1] + 1] = 1.0/R1 + C1/delta_time;
    staticValues[staticLineEnds[1] + 2] = -1.0/R1 - C1/delta_time;

    staticValues[staticLineEnds[2] + 1] = -1.0/R1 - C1/delta_time;
    staticValues[staticLineEnds[2] + 2] = 1.0/R1 + C1/delta_time + _A + _B + 1.0/R_U;
    staticValues[staticLineEnds[2] + 3] = -_A - _B - 1.0/R_U;

    staticValues[staticLineEnds[3] + 0] = -_A - _B - 1.0/R_U;
    staticValues[staticLineEnds[3] + 1] = _A + _B + 1.0/R_U + 1.0/R_B;

    jacobi.type = CSR;
    jacobi.size = BASIC_VALUE_NAMES_COUNT;
    jacobi.line_end = staticLineEnds;
    jacobi.index = staticIndices;
    jacobi.value = staticValues;

    // Формирование вектора значений
    minus_f[0] = -((basicValuesCurr[FI2] - basicValuesCurr[FI1]) - L1/delta_time * (basicValuesCurr[I_L1] - basicValuesPrev[I_L1]));
    minus_f[1] = -(basicValuesCurr[FI4] - L2/delta_time * (basicValuesCurr[I_L2] - basicValuesPrev[I_L2]));
    minus_f[2] = -(basicValuesCurr[I_E] - _I_R1 - _I_C1 - basicValuesCurr[I_L1]);
    minus_f[3] = -(_I_R1 + _I_C1 + basicValuesCurr[I_L1] - _I_C - _I - _I_R_U);
    minus_f[4] = -(_I_C + _I + _I_R_U - _I_R_B);
    minus_f[5] = -(_I_R_B + basicValuesCurr[I_L2] + _I_R2);
    minus_f[6] = -(basicValuesCurr[FI1] - _E);
}

/* Остальное */
// Операции с векторами
double *vectorCopy(double *vector1, const double *vector2, int size)
{
    for (int i = 0; i < size; i++)
    {
        vector1[i] = vector2[i];
    }

    return vector1;
}

double *vectorSub(double *vector1, const double *vector2, int size)
{
    for (int i = 0; i < size; i++)
    {
        vector1[i] -= vector2[i];
    }

    return vector1;
}

double *vectorAdd(double *vector1, const double *vector2, int size)
{
    for (int i = 0; i < size; i++)
    {
        vector1[i] += vector2[i];
    }

    return vector1;
}

double *vectorMul(double *vector, double number, int size)
{
    for (int i = 0; i < size; i++)
    {
        vector[i] *= number;
    }

    return vector;
}

double vectorNorm(const double *vector, int size)
{
    double max = 0;
    for (int i = 0; i < size; i++)
    {
        double abs = std::abs(vector[i]);
        if (abs > max)
        {
            max = abs;
        }
    }

    return max;
}

// Оценка ошибки моделирования
double estimateError(const double* const *values, const double *timePos, double delta_time)
{
    double max = 0.0;

    const double delta_timeLeft = timePos[PREV_TIME] - timePos[PREV_PREV_TIME];
    const double delta_timeRight = timePos[CURR_TIME] - timePos[PREV_TIME];
    for (int i = INT_VALUE_NAMES_COUNT; i < BASIC_VALUE_NAMES_COUNT; i++)
    {
        double delta_local = std::abs(delta_time*delta_time / 2.0 * (
            (values[CURR_TIME][i] - values[PREV_TIME][i]) / delta_timeRight -
            (values[PREV_TIME][i] - values[PREV_PREV_TIME][i]) / delta_timeLeft 
        ) / delta_timeLeft);

        if (delta_local > max)
        {
            max = delta_local;
        }
    }

    return max;
}

// Итерация метода Ньютона
double newtonIter(
    double *basicValuesCurr, const double *basicValuesPrev, 
    const double *reactiveValuesPrev, 
    double delta_time, double currTime
)
{
    Sparse jacobi;
    double minus_f[BASIC_VALUE_NAMES_COUNT];

    Sparse l, u;

    double delta_x[BASIC_VALUE_NAMES_COUNT];

    getNewtonModel(jacobi, minus_f, basicValuesCurr, basicValuesPrev, reactiveValuesPrev, delta_time, currTime);
    sparse_lu(&jacobi, &l, &u);
    solve_sparse_lu(&l, &u, minus_f, delta_x);

    vectorAdd(basicValuesCurr, delta_x, BASIC_VALUE_NAMES_COUNT);

    free_sparse(&l);
    free_sparse(&u);

    return vectorNorm(delta_x, BASIC_VALUE_NAMES_COUNT);
}

int main()
{
    double *basicValues[TIME_POS_COUNT_BASIC];
    double *reactiveValues[TIME_POS_COUNT_REACTIVE];
    double timePos[TIME_POS_COUNT_BASIC];
    for (int i = 0; i < TIME_POS_COUNT_BASIC; i++)
    {
        basicValues[i] = new double[BASIC_VALUE_NAMES_COUNT];
    }
    for (int i = 0; i < TIME_POS_COUNT_REACTIVE; i++)
    {
        reactiveValues[i] = new double[REACTIVE_VALUE_NAMES_COUNT];
    }

    double currTime = 0.0;
    double delta_time = DELTA_TIME_MIN;
    double timePrint = -1.0;

    int step = 0;
    bool isNextStep = false;

    while (currTime < TIME)
    {
        // Следующий момент
        if (isNextStep)
        {
            std::rotate(basicValues, basicValues + PREV_PREV_TIME, basicValues + TIME_POS_COUNT_BASIC);
            std::rotate(reactiveValues, reactiveValues + PREV_TIME, reactiveValues + TIME_POS_COUNT_REACTIVE);
            std::rotate(timePos, timePos + PREV_PREV_TIME, timePos + TIME_POS_COUNT_BASIC);
        }

        // Текущий момент
        timePos[CURR_TIME] = currTime;

        // Начальное приближение
        if (step >= PREV_PREV_TIME) // Если рассчитан предыдущий момент
        {
            const double timeCoef = delta_time / (timePos[PREV_TIME] - timePos[PREV_PREV_TIME]);
            for (int i = 0; i < BASIC_VALUE_NAMES_COUNT; i++)
            {
                basicValues[CURR_TIME][i] = basicValues[PREV_TIME][i] + (basicValues[PREV_TIME][i] - basicValues[PREV_PREV_TIME][i]) * timeCoef;
            }
        }
        else if (step == PREV_TIME) // Если рассчитан только предыдущий момент
        {
            vectorCopy(basicValues[CURR_TIME], basicValues[PREV_TIME], BASIC_VALUE_NAMES_COUNT);
        }
        else // Если на самом первом шаге
        {
            setInitialValues(basicValues, reactiveValues, timePos);
        }

        // Метод Ньютона
        try
        {
            for (int iter = 0; iter < NEWTON_MAX_ITER; iter++)
            {
                double delta_xNorm = newtonIter(basicValues[CURR_TIME], basicValues[PREV_TIME], reactiveValues[PREV_TIME], delta_time, currTime);

                if (delta_xNorm < NEWTON_EPSILON)
                {
                    throw true;
                }
            }

            // Не достигнута требуемая точность
            isNextStep = false;

            delta_time /= 2.0;
            if (delta_time < DELTA_TIME_MIN)
            {
                std::cerr << "*** Не удалось достигнуть требуемой точности метода Ньютона ***\n";
                return 1;
            }

            continue;
        }
        catch (bool ok) {} // Продолжение

        // Проверка на ошибку в расчетах
        for (int i = 0; i < BASIC_VALUE_NAMES_COUNT; i++)
        {
            if (std::isnan(basicValues[CURR_TIME][i]) || std::isinf(basicValues[CURR_TIME][i]))
            {
                std::cerr << "*** Ошибка в расчетах ***\n";
                return 3;
            }
        }

        // Оценка погрешности и подсчет временного шага
        if (step >= PREV_PREV_TIME) // Если рассчитан предыдущий момент предыдущего
        {
            const double delta_local = estimateError(basicValues, timePos, delta_time);
            if (delta_local < ERROR_DELTA_1)
            {
                currTime += delta_time;
                isNextStep = true;
                step++;

                delta_time *= 2.0;
            }
            else if (delta_local < ERROR_DELTA_2)
            {
                currTime += delta_time;
                isNextStep = true;
                step++;
            }
            else
            {
                isNextStep = false;

                delta_time /= 2.0;
                if (delta_time < DELTA_TIME_MIN)
                {
                    std::cerr << "*** Не удалось достигнуть требуемой точности моделирования ***\n";
                    return 2;
                }

                continue;
            }
        }
        else
        {
            currTime += delta_time;
            isNextStep = true;
            step++;
        }
        
        // Получение значений в реактивных элементах
        getReactiveValues(reactiveValues[CURR_TIME], reactiveValues[PREV_TIME], basicValues[CURR_TIME], basicValues[PREV_TIME], delta_time, timePos[CURR_TIME]);

        // Печать результатов
        if (timePos[CURR_TIME] > timePrint)
        {
            std::cout << timePos[CURR_TIME] << "    "; 
            for (int i = 0; i < BASIC_VALUES_TO_PRINT_COUNT; i++)
            {
                std::cout << basicValues[CURR_TIME][basicValuesToPrint[i]] << "    ";
            }
            for (int i = 0; i < REACTIVE_VALUES_TO_PRINT_COUNT; i++)
            {
                std::cout << reactiveValues[CURR_TIME][reactiveValuesToPrint[i]] << "    ";
            }
            std::cout << "\n";

            timePrint = timePos[CURR_TIME] + DELTA_TIME_PRINT;
        }
    }

    // Очистка памяти
    for (int i = 0; i < TIME_POS_COUNT_BASIC; i++)
    {
        delete[] basicValues[i];
    }
    for (int i = 0; i < TIME_POS_COUNT_REACTIVE; i++)
    {
        delete[] reactiveValues[i];
    }

    return 0;
}
