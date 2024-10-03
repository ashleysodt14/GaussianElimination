#include <stdio.h>
#include <signal.h>

#ifdef __linux__
#include <fenv.h>
#endif

// Floating-point exception handler prototype
void fpe_handler(int sig);

int main() {
    // Enable floating-point exceptions on Linux
    #ifdef __linux__
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    #endif

    // Set signal handler for SIGFPE
    void (*old_handler)(int) = signal(SIGFPE, fpe_handler);

    printf("Starting the main program...\n");

    // Example: trigger floating-point exception (division by zero)
    double x = 1.0 / 0.0;

    return 0;
}

// Floating-point exception handler function
void fpe_handler(int sig) {
    printf("Floating-point exception occurred! Signal: %d\n", sig);
}
