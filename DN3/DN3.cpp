#include <iostream>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <iomanip>

// Funkcija za izracun arctan(x) z uporabo Taylorjeve vrste
// x: kazalec na vhodno vrednost
// N_steps: kazalec na stevilo clenov v Taylorjevi vrsti
// Vrne izracunano vrednost arctan(x)
double izracunArctan(double* x, int* N_steps) {
    double rezultat = 0.0;
    double clen = *x; // Zacetek z x (n=0)
    for (int n = 0; n < *N_steps; ++n) {
        if (n > 0) {
            clen *= -(*x) * (*x); // Pomnozi z -x^2 za naslednji clen
        }
        rezultat += clen / (2 * n + 1);
    }
    // Klic funkcije arctan za testiranje z znanim kotom
    if (std::abs(*x - std::tan(M_PI / 4)) < 1e-9) {
        std::cout << "Klic funkcije arctan za znan kot x = tan(?/4): " << rezultat << std::endl;
    }
    return rezultat;
}

// Funkcija za izracun integrala z metodo trapezov
// a: spodnja meja integracije
// b: zgornja meja integracije
// n: stevilo podintervalov
// N_steps: stevilo clenov za Taylorjevo vrsto arctan
// Vrne izracunan integral
double trapeznaMetoda(double a, double b, int n, int N_steps) {
    double deltaX = (b - a) / n;
    double vsota = 0.0;

    // Zanka cez vse podintervale
    for (int i = 0; i <= n; ++i) {
        double x = a + i * deltaX;
        double polovicaX = x / 2;
        double fx = std::exp(3 * x) * izracunArctan(&polovicaX, &N_steps);

        // Teza funkcijskih vrednosti
        if (i == 0 || i == n) {
            vsota += fx; // Robne tocke
        }
        else {
            vsota += 2 * fx; // Notranje tocke
        }
    }

    // Pomnozi z deltaX/2 za koncno vrednost integrala
    return (deltaX / 2) * vsota;
}

int main() {
    // Parametri integracije
    double a = 0.0; // Spodnja meja
    double b = M_PI / 4; // Zgornja meja
    int n = 1000; // Stevilo podintervalov
    int N_steps = 10; // Stevilo clenov v Taylorjevi vrsti za arctan

    // Izracun integrala
    double rezultat = trapeznaMetoda(a, b, n, N_steps);

    // Izpis rezultata
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Izracunan integral je: " << rezultat << std::endl;

    return 0;
}
