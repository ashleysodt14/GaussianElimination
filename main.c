#include <stdio.h>
#include <signal.h>
#include <fenv.h>

// Prototype for the floating-point exception handler
void fpe_handler(int sig);

int main() {
    // Enable floating-point exceptions on Linux-based systems
    #ifdef __linux__
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    #endif

    // Set up signal handler for SIGFPE (Floating-point exception)
    void (*old_handler)(int) = signal(SIGFPE, fpe_handler);

    // Your main logic here, e.g., test cases or benchmarks
    printf("Starting the main program...\n");

    // Example floating-point exception
    double x = 1.0 / 0.0;  // Division by zero, just for testing

    return 0;
}

// Example floating-point exception handler
void fpe_handler(int sig) {
    printf("Floating-point exception occurred! Signal: %d\n", sig);
}
