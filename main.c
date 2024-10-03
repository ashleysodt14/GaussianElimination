#define _GNU_SOURCE

#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <fenv.h>
#include <setjmp.h>
#include <Python.h>  // Include Python API

#include "gauss_solve.h"
#include "helpers.h"

/* Size of the matrix */
#define N  3

void test_gauss_solve()
{
    printf("Entering function: %s\n", __func__);
  
    const double A0[N][N] = {
        {2, 3, -1},
        {4, 1, 2},
        {-2, 7, 2}
    };

    const double b0[N] = {5, 6, 3};
    double A[N][N], b[N], x[N], y[N];

    /* Create copies of the matrices. NOTE: the copies will get destroyed. */
    memcpy(A, A0, sizeof(A0));
    memcpy(b, b0, sizeof(b0));  

    gauss_solve_in_place(N, A, b);
    memcpy(x, b, sizeof(b0));
    matrix_times_vector(N, A0, x, y);

    double eps = 1e-6, dist = norm_dist(N, b0, y);
    assert(dist < eps);

    /* Print x */
    puts("x:\n");
    print_vector(N, x);

    /* Print U */
    puts("U:\n");
    print_matrix(N, A, FLAG_UPPER_PART);
  
    /* Print L */
    puts("L:\n");
    print_matrix(N, A, FLAG_LOWER_PART);
}

jmp_buf env;  // Buffer to store the state for setjmp/longjmp

void test_gauss_solve_with_zero_pivot()
{
    printf("Entering function: %s\n", __func__);
  
    double A[N][N] = {
        {0, 3, -1},
        {4, 1, 2},
        {-2, 7, 2}
    };

    double b[N] = {5, 6, 3};

    // Save the program state with setjmp
    if (setjmp(env) == 0) {
        gauss_solve_in_place(N, A, b);
        print_matrix(N, A, FLAG_LOWER_PART);
    } else {
        // This block is executed when longjmp is called
        printf("Returned to main program flow after exception\n");
    }
}

void test_lu_in_place()
{
    printf("Entering function: %s\n", __func__);

    const double A0[N][N] = {
        {2, 3, -1},
        {4, 1, 2},
        {-2, 7, 2}
    };

    const double b0[N] = {5, 6, 3};
    double A[N][N];

    memcpy(A, A0, sizeof(A0));
    lu_in_place(N, A);

    /* Print U */
    puts("U:\n");
    print_matrix(N, A, FLAG_UPPER_PART);
  
    /* Print L */
    puts("L:\n");
    print_matrix(N, A, FLAG_LOWER_PART);

    lu_in_place_reconstruct(N, A);

    /* Print U */
    puts("Reconstructed A:\n");
    print_matrix(N, A, FLAG_WHOLE);

    memcpy(A, A0, sizeof(A0));
    puts("Original A:\n");
    print_matrix(N, A, FLAG_WHOLE);

    double eps = 1e-6;
    assert(frobenius_norm_dist(N, A, A0) < eps);
}

void fpe_handler(int sig) {
    printf("Entering %s...\n", __func__);
    if(sig == SIGFPE) {
        printf("Floating point exception occurred, ignoring...\n");
        longjmp(env, 1);  // Jump back to where setjmp was called
    }
}

int main() {
    // Initialize the Python interpreter
    Py_Initialize();

    // Set the path to find your Python modules if needed
    PyRun_SimpleString("import sys; sys.path.append('.')");

    // Import the Python module and function
    PyObject *pName = PyUnicode_DecodeFSDefault("python_module");  // Replace "python_module" with your module name
    PyObject *pModule = PyImport_Import(pName);
    Py_XDECREF(pName);

    if (pModule != NULL) {
        // Get the function from the module
        PyObject *pFunc = PyObject_GetAttrString(pModule, "python_function");  // Replace "python_function" with your function name
        
        if (PyCallable_Check(pFunc)) {
            // Prepare arguments to pass to the Python function
            PyObject *pArgs = PyTuple_New(1);  // For example, a single argument
            PyObject *pValue = PyLong_FromLong(10);  // Example argument
            PyTuple_SetItem(pArgs, 0, pValue);

            // Call the Python function
            PyObject *pResult = PyObject_CallObject(pFunc, pArgs);

            // Check the result and handle it (assuming the function returns an integer)
            if (pResult != NULL) {
                printf("Result from Python: %ld\n", PyLong_AsLong(pResult));
                Py_XDECREF(pResult);
            } else {
                PyErr_Print();
            }

            Py_XDECREF(pArgs);
        } else {
            PyErr_Print();
        }

        Py_XDECREF(pFunc);
        Py_XDECREF(pModule);
    } else {
        PyErr_Print();
    }

    // Finalize the Python interpreter
    Py_Finalize();

    return 0;
}
