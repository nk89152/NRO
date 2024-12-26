#define _USE_MATH_DEFINES
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <omp.h>

using namespace std;

int main() {
    // Inicializiramo matriko A in vektor b
    vector<vector<double>> A;
    vector<double> b;

    // Podamo ime datoteke
    string filename = "./datoteka_A_b.txt";

    // Preberemo datoteko
    ifstream infile(filename);
    if (!infile) {
        cerr << "Napaka pri odpiranju datoteke." << endl;
        return 1;
    }

    // Preberemo velikost matrike A
    string string_first_line;
    getline(infile, string_first_line);
    replace(string_first_line.begin(), string_first_line.end(), '=', ' ');
    istringstream iss(string_first_line);
    string nepomemben_del1, nepomemben_del2;
    int n;
    iss >> nepomemben_del1 >> nepomemben_del2 >> n;

    cout << "Velikost matrike A: " << n << "x" << n << endl;

    // Branje matrike A
    for (int iiA = 0; iiA < n; iiA++) {
        string line;
        getline(infile, line);
        replace(line.begin(), line.end(), ';', ' ');
        istringstream iss_column(line);
        vector<double> row(n);
        for (int column = 0; column < n; column++) {
            iss_column >> row[column];
        }
        A.push_back(row);
    }

    // Presko�imo prazno vrstico
    string empty_line;
    getline(infile, empty_line);

    // Preberemo vektor b
    string string_line_b;
    getline(infile, string_line_b);
    replace(string_line_b.begin(), string_line_b.end(), '>', ' ');
    istringstream iss_b(string_line_b);
    int n_b;
    iss_b >> nepomemben_del1 >> nepomemben_del2 >> n_b;

    cout << "Velikost vektorja b: " << n_b << endl;

    for (int iib = 0; iib < n_b; iib++) {
        string line_b_element;
        getline(infile, line_b_element);
        istringstream iss_b_element(line_b_element);
        double b_element = 0;
        iss_b_element >> b_element;
        b.push_back(b_element);
    }

    // Inicializiramo vektor re�itve T in SOR parameter
    vector<double> T(n, 100.0);  // Vsi elementi vektorja T so na za�etku 100
    vector<double> T_old(n, 100.0);  // Kopija za preverjanje konvergence
    double max_error = 1e-10;  // Najve�ja dovoljena napaka (kriterij za konvergenco)
    double omega = 1.5;  // Relaksacijski faktor za SOR metodo (1.0 pomeni navadno Gauss-Seidel metodo)
    bool converged = false;  // Indikator, ali je metoda dosegla konvergenco

    // Za�etek merjenja �asa
    auto start_time = chrono::high_resolution_clock::now();

    // Iterativni postopek za re�evanje sistema ena�b
    for (int iter = 0; iter < 2000 && !converged; iter++) {
        converged = true;

        // Gauss-Seidel - rde�i vozli (parni indeksi)
        // Ker Gauss-Seidel metoda uporablja rezultate iz iste iteracije, je te�ko paralelizirati.
        // Da to omogo�imo, lo�imo vozli��a na parna (rde�a) in neparna (�rna).
#pragma omp parallel for schedule(dynamic) shared(A, b, T, n)
        for (int i = 0; i < n; i++) {
            if (i % 2 == 0) {  // Rde�i vozli
                double sum = b[i];
                for (int j = 0; j < n; j++) {
                    if (j != i) sum -= A[i][j] * T[j];  // Izra�un prispevka drugih �lenov
                }
                double T_new = sum / A[i][i];  // Delitev z diagonalo
                T[i] = omega * T_new + (1 - omega) * T[i];  // SOR posodobitev
            }
        }

        // Gauss-Seidel - �rni vozli (neparni indeksi)
#pragma omp parallel for schedule(dynamic) shared(A, b, T, n)
        for (int i = 0; i < n; i++) {
            if (i % 2 != 0) {  // �rni vozli
                double sum = b[i];
                for (int j = 0; j < n; j++) {
                    if (j != i) sum -= A[i][j] * T[j];  // Izra�un prispevka drugih �lenov
                }
                double T_new = sum / A[i][i];  // Delitev z diagonalo
                T[i] = omega * T_new + (1 - omega) * T[i];  // SOR posodobitev
            }
        }

        // Preverjanje konvergence
        // Primerjamo novo in staro vrednost vektorja T, da ugotovimo, ali smo dosegli �eleno natan�nost
#pragma omp parallel for schedule(dynamic) shared(T, T_old, converged)
        for (int i = 0; i < n; i++) {
            if (abs(T[i] - T_old[i]) > max_error) {
                converged = false;  // �e razlika presega dovoljeno, �e nismo konvergirali
            }
            T_old[i] = T[i];  // Posodobimo staro vrednost za naslednjo iteracijo
        }
    }

    // Konec merjenja �asa
    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> time_duration = end_time - start_time;

    // Izpis �asa izvajanja
    cout << "�as izvajanja Gauss-Seidelove metode: " << time_duration.count() << " sekund" << endl;

    // Najvi�ja temperatura
    double max_T = *max_element(T.begin(), T.end());
    cout << "Najvisja temperatura: " << max_T << " �C" << endl;

    return 0;
}
// Gauss-Seidel metoda ni naravno paralelizirana, saj nove vrednosti vektorja T 
// temeljijo na �e posodobljenih vrednostih iz iste iteracije. To povzro�a serijske
// odvisnosti, ki prepre�ujejo enostavno vzporedno izvajanje.

// Za omogo�anje paralelizacije smo metodo prilagodili z razdelitvijo vozli�� na
// rde�e (parni indeksi) in �rne (neparni indeksi) vozle. Znotraj vsake skupine
// (rde�e ali �rne) so izra�uni neodvisni, kar omogo�a, da se vozli v isti skupini
// izra�unavajo vzporedno z uporabo OpenMP.

// Uvedli smo tudi metodo Successive Over-Relaxation (SOR) za pospe�itev konvergence.
// Relaksacijski faktor omega (? > 1) omogo�a hitrej�o konvergenco, vendar je optimalna
// vrednost ? dolo�ena eksperimentalno, saj mora zagotavljati stabilnost in u�inkovitost metode.
