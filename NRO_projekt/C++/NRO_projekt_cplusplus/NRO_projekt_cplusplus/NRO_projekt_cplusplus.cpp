// Projektna naloga pri predmetu Napredna raèunalniška orodja

// Uvoz knjižnic
#include <iostream> // Za vhodno/izhodne operacije
#include <fstream>  // Za delo z datotekami
#include <sstream>  // Za pretvorbo nizov v druge tipe
#include <vector>   // Za uporabo dinamiènih tabel
#include <cmath>    // Za matematiène funkcije
#include <algorithm> // Za algoritme, kot je std::replace
#include <chrono>   // Za merjenje èasa
#include <omp.h>    // Za uporabo OpenMP za paralelizacijo

using namespace std;

void printVector(const std::vector<int>& vec) {
    for (const int& val : vec) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}

void printVector(const std::vector<double>& vec) {
    for (const double& val : vec) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}

void print2DVector(const std::vector<std::vector<int>>& vec) {
    for (const auto& innerVec : vec) {
        for (const int& val : innerVec) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
}

int main() {
    // Definiranje spremenljivk
    std::string mreza = "./primer4mreza.txt"; // Pot do datoteke z mrežo

    // Vektorji za shranjevanje podatkov
    std::vector<double> X; // X koordinate vozlišè
    std::vector<double> Y; // Y koordinate vozlišè
    std::vector<vector<int>> celice; // Celice, ki jih sestavljajo vozlišèa
    std::vector<vector<int>> vozlisca_robnih_pogojev; // Vozlišèa robnih pogojev
    std::vector<int> tipi_robnih_pogojev; // Tipi robnih pogojev
    std::vector<double> vrednosti_robnih_pogojev; // Vrednosti robnih pogojev
    std::vector<double> vrednosti_prestopa_toplote; // Vrednosti toplotnega prestopa

    // Odpri datoteko
    std::ifstream file;
    file.open(mreza);

    // Preberemo prvo vrstico, ki vsebuje število vozlišè
    std::string string_first_line;
    std::getline(file, string_first_line);
    std::istringstream iss(string_first_line);
    std::string del1;

    // Branje števila vozlišè
    int n_vozlisc; // Število vozlišè
    iss >> del1; // Preskoèi prvi del
    iss >> n_vozlisc; // Preberemo število vozlišè

    // Branje vozlišè
    for (int i = 0; i < n_vozlisc; i++) {
        std::string s;
        std::getline(file, s); // Preberemo vrstico
        std::replace(s.begin(), s.end(), ';', ' '); // Zamenjaj ';' z ' '
        std::replace(s.begin(), s.end(), ',', ' '); // Zamenjaj ',' z ' '
        std::istringstream iss(s);
        int node_id;
        double x;
        double y;
        iss >> node_id >> x >> y; // Preberemo ID vozlišèa in njegove koordinate

        // Dodamo koordinate v vektorje
        X.push_back(x);
        Y.push_back(y);
    }

    // Preberemo prazno vrstico
    std::string prazna_vrstica;
    std::getline(file, prazna_vrstica);

    // Branje celic
    std::string cell_line;
    std::getline(file, cell_line); // Preberemo vrstico s celicami

    std::vector<std::string> cells_string;
    std::string cell;
    std::istringstream iss2(cell_line);
    while (iss2 >> cell) {
        cells_string.push_back(cell); // Shranimo celice v vektor
    }
    int n_celice = std::stoi(cells_string[1]); // Preberemo število celic

    // Branje celic
    for (int i = 0; i < n_celice; i++) {
        std::string s;
        std::getline(file, s); // Preberemo vrstico
        std::replace(s.begin(), s.end(), ';', ' '); // Zamenjaj ';' z ' '
        std::replace(s.begin(), s.end(), ',', ' '); // Zamenjaj ',' z ' '
        std::istringstream iss(s);
        int cell_id;
        int node1_id;
        int node2_id;
        int node3_id;
        int node4_id;
        iss >> cell_id >> node1_id >> node2_id >> node3_id >> node4_id; // Preberemo ID celice in njene vozlišèa
        vector<int> celica;

        // Dodamo vozlišèa celice v vektor
        celica.push_back(node1_id);
        celica.push_back(node2_id);
        celica.push_back(node3_id);
        celica.push_back(node4_id);

        celice.push_back(celica); // Dodamo celico v vektor celic
    }

    // Branje robnih pogojev
    std::string prazna_vrstica2;
    std::getline(file, prazna_vrstica2); // Preberemo prazno vrstico

    std::string robni_pogoji_line;
    std::getline(file, robni_pogoji_line); // Preberemo vrstico z robnimi pogoji

    std::vector<std::string> robni_pogoji_string;
    std::string rp;
    std::istringstream iss3(robni_pogoji_line);

    int n_pogoji = 0; // Število robnih pogojev
    std::string del2;
    std::string del3;

    iss3 >> del2; // Preskoèi prvi del
    iss3 >> del3; // Preskoèi drugi del
    iss3 >> n_pogoji; // Preberemo število robnih pogojev

    // Branje robnih pogojev
    for (int i = 0; i < n_pogoji; i++) {
        std::string s;
        string tip_pogoja;
        std::getline(file, s); // Preberemo vrstico
        std::istringstream iss4(s);
        std::string npd1;
        std::string npd2;
        iss4 >> npd1; // Preskoèi prvi del
        iss4 >> npd2; // Preskoèi drugi del
        iss4 >> tip_pogoja; // Preberemo tip robnega pogoja

        // Branje datoteke, èe imamo opravka z znano temperaturo na robu
        if (tip_pogoja == "temperatura") {
            tipi_robnih_pogojev.push_back(0); // Dodamo tip robnega pogoja
            std::string rp_temp;
            std::getline(file, rp_temp); // Preberemo temperaturo
            std::istringstream issrp1(rp_temp);
            double temperatura;
            std::string del4;

            issrp1 >> del4; // Preskoèi del
            issrp1 >> temperatura; // Preberemo temperaturo

            vrednosti_robnih_pogojev.push_back(temperatura); // Dodamo temperaturo v vektor
            vrednosti_prestopa_toplote.push_back(-1); // Dodamo placeholder za toplotni prestop
        }
        // Branje datoteke, èe imamo opravka z toplotnim tokom na robu
        else if (tip_pogoja == "toplotni") {
            tipi_robnih_pogojev.push_back(1); // Dodamo tip robnega pogoja
            std::string rp_q;
            std::getline(file, rp_q); // Preberemo toplotni tok
            std::istringstream issrp2(rp_q);
            double toplotni_tok = 0;

            std::string del4;
            issrp2 >> del4; // Preskoèi del
            issrp2 >> toplotni_tok; // Preberemo toplotni tok

            vrednosti_robnih_pogojev.push_back(toplotni_tok); // Dodamo toplotni tok v vektor
        }

        int st_vozlisc_v_rp = 0; // Število vozlišè v robnem pogoju
        std::string st_vozlisc_rp;
        std::getline(file, st_vozlisc_rp); // Preberemo število vozlišè v robnem pogoju
        std::istringstream strp(st_vozlisc_rp);
        strp >> st_vozlisc_v_rp; // Preberemo število vozlišè

        vector<int> vozlisca_v_robnem_pogoju; // Vektor za shranjevanje vozlišè robnega pogoja

        for (int i = 0; i < st_vozlisc_v_rp; i++) {
            std::string id_voz_rp;
            std::getline(file, id_voz_rp); // Preberemo ID vozlišèa
            std::istringstream idvozrp(id_voz_rp);
            int id_vozlisca = 0;
            idvozrp >> id_vozlisca; // Preberemo ID vozlišèa

            vozlisca_v_robnem_pogoju.push_back(id_vozlisca); // Dodamo vozlišèe v vektor
        }

        vozlisca_robnih_pogojev.push_back(vozlisca_v_robnem_pogoju); // Dodamo vozlišèa robnega pogoja v vektor

        std::string prazna_vrstica2;
        std::getline(file, prazna_vrstica2); // Preberemo prazno vrstico
    }

    // Izpis vseh vektorjev
    std::cout << "Tipi robnih pogojev: ";
    printVector(tipi_robnih_pogojev);

    std::cout << "Vrednosti robnih pogojev: ";
    printVector(vrednosti_robnih_pogojev);

    std::cout << "Vrednosti prestopa toplote: ";
    printVector(vrednosti_prestopa_toplote);

    std::cout << "Vozlišèa robnih pogojev: " << std::endl;
    print2DVector(vozlisca_robnih_pogojev);

    // Konec branja datoteke

    // Sledi blok kode, ki preveri sosednjost vozlišè in to zapiše v obliki vektorjev
    double deltaX = 1.0; // Razdalja med vozlišèi v X smeri
    double deltaY = 1.0; // Razdalja med vozlišèi v Y smeri
    double k = 24.0; // Koeficient toplotne prevodnosti

    vector<vector<int>> sosednja_vozlisca; // Vektor za shranjevanje sosednjih vozlišè

    for (int node_id = 0; node_id < n_vozlisc; node_id++) {
        vector<int> node_i_neighbours = { -1,-1,-1,-1 }; // Inicializacija sosednjih vozlišè

        for (int nd = 0; nd < n_celice; nd++) {
            vector<int> trenutna_celica = celice[nd]; // Pridobimo trenutne celice

            int vozlisce1 = trenutna_celica[0];
            int vozlisce2 = trenutna_celica[1];
            int vozlisce3 = trenutna_celica[2];
            int vozlisce4 = trenutna_celica[3];

            // Preverimo, ali je vozlišèe del trenutne celice
            if (node_id == vozlisce1 || node_id == vozlisce2 || node_id == vozlisce3 || node_id == vozlisce4) {
                for (int vozl = 0; vozl < 4; vozl++) {
                    int sosednje_vozlisce = trenutna_celica[vozl]; // Pridobimo sosednje vozlišèe

                    int pozicija = 0; // Inicializacija pozicije

                    // Preverimo, ali je sosednje vozlišèe razliène od trenutnega vozlišèa
                    if (sosednje_vozlisce != node_id) {
                        double x_obravnavano_vozl = X[node_id]; // X koordinata trenutnega vozlišèa
                        double y_obravnavano_vozl = Y[node_id]; // Y koordinata trenutnega vozlišèa

                        double x_sosed = X[sosednje_vozlisce]; // X koordinata sosednjega vozlišèa
                        double y_sosed = Y[sosednje_vozlisce]; // Y koordinata sosednjega vozlišèa

                        // Doloèimo pozicijo sosednjega vozlišèa glede na trenutnega
                        if ((x_obravnavano_vozl - x_sosed) < 1e-9 && (x_obravnavano_vozl - x_sosed) > -(1e-9)) {
                            if ((y_obravnavano_vozl - y_sosed) > 0.0) {
                                pozicija = 1; // Sosednje vozlišèe je zgoraj
                            }
                            else {
                                pozicija = 3; // Sosednje vozlišèe je spodaj
                            }
                        }
                        else if ((y_obravnavano_vozl - y_sosed) < 1e-9 && (y_obravnavano_vozl - y_sosed) > -(1e-9)) {
                            if ((x_obravnavano_vozl - x_sosed) > 0.0) {
                                pozicija = 0; // Sosednje vozlišèe je levo
                            }
                            else {
                                pozicija = 2; // Sosednje vozlišèe je desno
                            }
                        }
                        else {
                            pozicija = -1; // Sosednje vozlišèe ni v pravilni poziciji
                        }

                        // Dodamo sosednje vozlišèe v ustrezno pozicijo
                        if (pozicija != -1) {
                            node_i_neighbours[pozicija] = sosednje_vozlisce;
                        }
                    }
                }
            }
        }
        sosednja_vozlisca.push_back(node_i_neighbours); // Dodamo sosednja vozlišèa v vektor
    }

    std::cout << "Vektor sosednjih vozlišè: " << std::endl;
    print2DVector(sosednja_vozlisca);

    int n = n_vozlisc; // Število vozlišè

    // Ustvarimo matriko koeficientov A in vektor neznank b
    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0)); // Matriko A inicializiramo z nièlami
    std::vector<double> b(n, 0.0); // Vektor b inicializiramo z nièlami

    // Izpolnimo matriko A in vektor b
    for (int node_id = 0; node_id < n_vozlisc; node_id++) {
        std::vector<int> sosedi = sosednja_vozlisca[node_id]; // Pridobimo sosednja vozlišèa
        int levi_sosed = sosedi[0];
        int spodnji_sosed = sosedi[1];
        int desni_sosed = sosedi[2];
        int zgornji_sosed = sosedi[3];

        // Preverimo, ali je vozlišèe v notranjosti
        if (levi_sosed != -1 && spodnji_sosed != -1 && desni_sosed != -1 && zgornji_sosed != -1) {
            A[node_id][levi_sosed] = 1.0; // Doloèimo koeficiente za sosednje vozlišèa
            A[node_id][spodnji_sosed] = 1.0;
            A[node_id][desni_sosed] = 1.0;
            A[node_id][zgornji_sosed] = 1.0;
            A[node_id][node_id] = -4.0; // Doloèimo koeficient za trenutno vozlišèe
        }
        else { // Vozlišèa na robu
            int tip_robnega_pogoja = 0; // Inicializiramo tip robnega pogoja
            double vrednost = 0; // Inicializiramo vrednost

            // Preverimo robne pogoje
            for (int robni_pogoj_id = 0; robni_pogoj_id < 5; robni_pogoj_id++) {
                std::vector<int> vozlisca_v_trenutnem_rp = vozlisca_robnih_pogojev[robni_pogoj_id];

                for (int id_vozlisce_rp = 0; id_vozlisce_rp < vozlisca_v_trenutnem_rp.size(); id_vozlisce_rp++) {
                    int vozlisce_v_trenutnem_rp = vozlisca_v_trenutnem_rp[id_vozlisce_rp];

                    if (node_id == vozlisce_v_trenutnem_rp) {
                        tip_robnega_pogoja = tipi_robnih_pogojev[robni_pogoj_id]; // Doloèimo tip robnega pogoja
                        vrednost = vrednosti_robnih_pogojev[robni_pogoj_id]; // Doloèimo vrednost robnega pogoja
                    }
                }
            }

            // Obdelava Dirichletovega robnega pogoja
            if (tip_robnega_pogoja == 0) {
                A[node_id][node_id] = 1.0; // Doloèimo koeficient za robni pogoj
                b[node_id] = vrednost; // Doloèimo vrednost vektorja b
            }
            // Obdelava Neumannovega robnega pogoja
            else if (tip_robnega_pogoja == 1) {
                int stevilo_sosedov = 0; // Število sosedov

                for (int st = 0; st < 4; st++) {
                    if (sosedi[st] != -1) {
                        stevilo_sosedov++; // Preštejemo sosednje vozlišèa
                    }
                }

                // Doloèimo koeficiente za robne pogoje
                if (stevilo_sosedov == 3) {
                    if (levi_sosed == -1) {
                        A[node_id][node_id] -= 4.0;
                        A[node_id][desni_sosed] += 2.0;
                        A[node_id][spodnji_sosed] += 1.0;
                        A[node_id][zgornji_sosed] += 1.0;
                        b[node_id] = -2.0 * (vrednost * deltaX / k);
                    }
                    if (desni_sosed == -1) {
                        A[node_id][node_id] -= 4.0;
                        A[node_id][levi_sosed] += 2.0;
                        A[node_id][spodnji_sosed] += 1.0;
                        A[node_id][zgornji_sosed] += 1.0;
                        b[node_id] = -2.0 * (vrednost * deltaX / k);
                    }
                    if (spodnji_sosed == -1) {
                        A[node_id][node_id] -= 4.0;
                        A[node_id][levi_sosed] += 1.0;
                        A[node_id][desni_sosed] += 1.0;
                        A[node_id][zgornji_sosed] += 2.0;
                        b[node_id] = -2.0 * (vrednost * deltaX / k);
                    }
                    if (zgornji_sosed == -1) {
                        A[node_id][node_id] -= 4.0;
                        A[node_id][levi_sosed] += 1.0;
                        A[node_id][desni_sosed] += 1.0;
                        A[node_id][spodnji_sosed] += 2.0;
                        b[node_id] = -2.0 * (vrednost * deltaX / k);
                    }
                }
            }
        }
    }

    // Shranjevanje matrike A in vektorja b za uporabo v MATLABu
    std::ofstream fileIDA;
    fileIDA.open("A.csv"); // Odpri datoteko za shranjevanje matrike A

    for (int i = 0; i < A.size(); i++) {
        fileIDA << A[i][0]; // Zapiši prvi element vrstice

        for (int j = 1; j < A[i].size(); j++) {
            fileIDA << ", " << A[i][j]; // Zapiši preostale elemente vrstice
        }

        fileIDA << std::endl; // Prehod na novo vrstico
    }
    fileIDA.close(); // Zapri datoteko

    std::ofstream fileIDb;
    fileIDb.open("b.csv"); // Odpri datoteko za shranjevanje vektorja b

    for (int i = 0; i < b.size(); i++) {
        fileIDb << b[i]; // Zapiši element vektorja b
        if (i != b.size() - 1) {
            fileIDb << ","; // Dodaj vejico, razen pri zadnjem elementu
        }
        fileIDb << std::endl; // Prehod na novo vrstico
    }
    fileIDb.close(); // Zapri datoteko

    // Inicializacija vektorja T za shranjevanje rezultatov
    vector<double> T(n, 100); // Inicializiraj vse temperature na 100
    vector<double> T_old = T; // Ustvari kopijo za preverjanje konvergence

    double max_error = 1e-10; // Toleranca za preverjanje konvergence
    double omega = 1.5; // Successive Over-Relaxation (SOR) parameter
    bool converged = false; // Spremenljivka za spremljanje konvergence

    int st_iter = 1000; // Najveèje število iteracij

    // Nastavitev števila niti
    omp_set_num_threads(16); // Nastavi število niti za delovanje
    omp_set_dynamic(1); // Omogoèi dinamièno dodeljevanje niti

    // Zaèetek merjenja èasa sistema enaèb
    auto start_time = std::chrono::high_resolution_clock::now();

    // Gauss-Seidel metoda s SOR
    for (int iter = 0; iter < st_iter && !converged; iter++) {
        converged = true;

    //Successive Over-Relaxation (SOR) je iterativna metoda za reševanje linearnih sistemov, 
    //ki izboljša hitrost konvergence Gauss-Seidelove metode z uporabo relaksacijskega faktorja 
    //omega, ki prilagaja posodobitve vrednosti, da se pospeši konvergenca in zmanjša število potrebnih iteracij

        // Rdeèi vozli (parni indeksi)
#pragma omp parallel for schedule(dynamic) shared(A, b, T, n) reduction(&&: converged)
        for (int jj = 0; jj < n; jj++) {
            if (jj % 2 == 0) { // Rdeèi vozli
                double d = b[jj];
                for (int ii = 0; ii < n; ii++) {
                    if (jj != ii) {
                        d -= A[jj][ii] * T[ii]; // Izraèunaj novo vrednost
                    }
                }
                double T_new = d / A[jj][jj]; // Izraèunaj novo temperaturo
                double T_updated = omega * T_new + (1 - omega) * T[jj]; // Uporabi SOR
                if (abs(T_updated - T_old[jj]) > max_error) {
                    converged = false;
                }
                T[jj] = T_updated;
            }
        }

        // Èrni vozli (neparni indeksi)
#pragma omp parallel for schedule(dynamic) shared(A, b, T, n) reduction(&&: converged)
        for (int jj = 0; jj < n; jj++) {
            if (jj % 2 != 0) { // Èrni vozli
                double d = b[jj];
                for (int ii = 0; ii < n; ii++) {
                    if (jj != ii) {
                        d -= A[jj][ii] * T[ii]; // Izraèunaj novo vrednost
                    }
                }
                double T_new = d / A[jj][jj]; // Izraèunaj novo temperaturo
                double T_updated = omega * T_new + (1 - omega) * T[jj]; // Uporabi SOR
                if (abs(T_updated - T_old[jj]) > max_error) {
                    converged = false;
                }
                T[jj] = T_updated;
            }
        }

        // Posodobi stari vektor za naslednjo iteracijo
        T_old = T;
    }

    // Konec merjenja èasa
    auto end_time = std::chrono::high_resolution_clock::now();
    
    // Izpis trajanja
    std::chrono::duration<double> time_duration = end_time - start_time;
    std::cout << "Èas Gauss-Seidl metode: " << time_duration.count() << " sekund" << std::endl;

    // Zapis rezultatov v .vtk datoteki
    std::ofstream fileID;
    fileID.open("rezultat.vtk"); // Odpri datoteko za shranjevanje rezultatov

    fileID << "# vtk DataFile Version 3.0\n"; // Zapis verzije VTK datoteke
    fileID << "Mesh_1\n"; // Ime mreže
    fileID << "ASCII\n"; // Oblika datoteke
    fileID << "DATASET UNSTRUCTURED_GRID\n"; // Tip podatkov
    fileID << "POINTS " << n_vozlisc << " float\n"; // Zapis števila vozlišè
    for (int koordinata_id = 0; koordinata_id < n_vozlisc; koordinata_id++) {
        fileID << X[koordinata_id] << " " << Y[koordinata_id] << " 0\n"; // Zapis koordinat vozlišè
    }
    fileID << "\n";
    fileID << "CELLS " << n_celice << " " << n_celice * 5 << "\n"; // Zapis števila celic
    for (int celica_id = 0; celica_id < n_celice; celica_id++) {
        int vozl_id1 = celice[celica_id][0];
        int vozl_id2 = celice[celica_id][1];
        int vozl_id3 = celice[celica_id][2];
        int vozl_id4 = celice[celica_id][3];
        fileID << "4 " << vozl_id1 << " " << vozl_id2 << " " << vozl_id3 << " " << vozl_id4 << "\n"; // Zapis celic
    }
    fileID << "\n";
    fileID << "CELL_TYPES " << n_celice << "\n"; // Zapis tipov celic
    for (int celica_id = 0; celica_id < n_celice; celica_id++) {
        fileID << "9\n"; // Tip celice (kvadrat)
    }
    fileID << "\n";
    fileID << "POINT_DATA " << n_vozlisc << "\n"; // Zapis podatkov o vozlišèih
    fileID << "SCALARS Temperature float 1\n"; // Zapis temperature kot skalarne vrednosti
    fileID << "LOOKUP_TABLE default\n"; // Doloèitev privzete tabele
    for (int koordinata_id = 0; koordinata_id < n_vozlisc; koordinata_id++) {
        fileID << T[koordinata_id] << "\n"; // Zapis temperature za vsako vozlišèe
    }

    fileID.close(); // Zapri datoteko

    return 0; // Konec programa
}