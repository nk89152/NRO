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

    // Preskoèimo prazno vrstico
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

    // Inicializiramo vektor rešitve T in SOR parameter
    vector<double> T(n, 100.0);  // Vsi elementi vektorja T so na zaèetku 100
    vector<double> T_old(n, 100.0);  // Kopija za preverjanje konvergence
    double max_error = 1e-10;  // Najveèja dovoljena napaka (kriterij za konvergenco)
    double omega = 1.5;  // Relaksacijski faktor za SOR metodo (1.0 pomeni navadno Gauss-Seidel metodo)
    bool converged = false;  // Indikator, ali je metoda dosegla konvergenco

    // Zaèetek merjenja èasa
    auto start_time = chrono::high_resolution_clock::now();

    // Iterativni postopek za reševanje sistema enaèb
    for (int iter = 0; iter < 2000 && !converged; iter++) {
        converged = true;

        // Gauss-Seidel - rdeèi vozli (parni indeksi)
        // Ker Gauss-Seidel metoda uporablja rezultate iz iste iteracije, je težko paralelizirati.
        // Da to omogoèimo, loèimo vozlišèa na parna (rdeèa) in neparna (èrna).
#pragma omp parallel for schedule(dynamic) shared(A, b, T, n)
        for (int i = 0; i < n; i++) {
            if (i % 2 == 0) {  // Rdeèi vozli
                double sum = b[i];
                for (int j = 0; j < n; j++) {
                    if (j != i) sum -= A[i][j] * T[j];  // Izraèun prispevka drugih èlenov
                }
                double T_new = sum / A[i][i];  // Delitev z diagonalo
                T[i] = omega * T_new + (1 - omega) * T[i];  // SOR posodobitev
            }
        }

        // Gauss-Seidel - èrni vozli (neparni indeksi)
#pragma omp parallel for schedule(dynamic) shared(A, b, T, n)
        for (int i = 0; i < n; i++) {
            if (i % 2 != 0) {  // Èrni vozli
                double sum = b[i];
                for (int j = 0; j < n; j++) {
                    if (j != i) sum -= A[i][j] * T[j];  // Izraèun prispevka drugih èlenov
                }
                double T_new = sum / A[i][i];  // Delitev z diagonalo
                T[i] = omega * T_new + (1 - omega) * T[i];  // SOR posodobitev
            }
        }

        // Preverjanje konvergence
        // Primerjamo novo in staro vrednost vektorja T, da ugotovimo, ali smo dosegli želeno natanènost
#pragma omp parallel for schedule(dynamic) shared(T, T_old, converged)
        for (int i = 0; i < n; i++) {
            if (abs(T[i] - T_old[i]) > max_error) {
                converged = false;  // Èe razlika presega dovoljeno, še nismo konvergirali
            }
            T_old[i] = T[i];  // Posodobimo staro vrednost za naslednjo iteracijo
        }
    }

    // Konec merjenja èasa
    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> time_duration = end_time - start_time;

    // Izpis èasa izvajanja
    cout << "Èas izvajanja Gauss-Seidelove metode: " << time_duration.count() << " sekund" << endl;

    // Najvišja temperatura
    double max_T = *max_element(T.begin(), T.end());
    cout << "Najvisja temperatura: " << max_T << " °C" << endl;

    return 0;
}
// Gauss-Seidel metoda ni naravno paralelizirana, saj nove vrednosti vektorja T 
// temeljijo na že posodobljenih vrednostih iz iste iteracije. To povzroèa serijske
// odvisnosti, ki prepreèujejo enostavno vzporedno izvajanje.

// Za omogoèanje paralelizacije smo metodo prilagodili z razdelitvijo vozlišè na
// rdeèe (parni indeksi) in èrne (neparni indeksi) vozle. Znotraj vsake skupine
// (rdeèe ali èrne) so izraèuni neodvisni, kar omogoèa, da se vozli v isti skupini
// izraèunavajo vzporedno z uporabo OpenMP.

// Uvedli smo tudi metodo Successive Over-Relaxation (SOR) za pospešitev konvergence.
// Relaksacijski faktor omega (? > 1) omogoèa hitrejšo konvergenco, vendar je optimalna
// vrednost ? doloèena eksperimentalno, saj mora zagotavljati stabilnost in uèinkovitost metode.
