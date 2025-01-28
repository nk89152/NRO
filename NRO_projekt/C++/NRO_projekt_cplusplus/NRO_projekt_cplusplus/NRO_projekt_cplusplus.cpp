// Projektna naloga pri predmetu Napredna ra�unalni�ka orodja

// Uvoz knji�nic
#include <iostream> // Za vhodno/izhodne operacije
#include <fstream>  // Za delo z datotekami
#include <sstream>  // Za pretvorbo nizov v druge tipe
#include <vector>   // Za uporabo dinami�nih tabel
#include <cmath>    // Za matemati�ne funkcije
#include <algorithm> // Za algoritme, kot je std::replace
#include <chrono>   // Za merjenje �asa
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
    std::string mreza = "./primer4mreza.txt"; // Pot do datoteke z mre�o

    // Vektorji za shranjevanje podatkov
    std::vector<double> X; // X koordinate vozli��
    std::vector<double> Y; // Y koordinate vozli��
    std::vector<vector<int>> celice; // Celice, ki jih sestavljajo vozli��a
    std::vector<vector<int>> vozlisca_robnih_pogojev; // Vozli��a robnih pogojev
    std::vector<int> tipi_robnih_pogojev; // Tipi robnih pogojev
    std::vector<double> vrednosti_robnih_pogojev; // Vrednosti robnih pogojev
    std::vector<double> vrednosti_prestopa_toplote; // Vrednosti toplotnega prestopa

    // Odpri datoteko
    std::ifstream file;
    file.open(mreza);

    // Preberemo prvo vrstico, ki vsebuje �tevilo vozli��
    std::string string_first_line;
    std::getline(file, string_first_line);
    std::istringstream iss(string_first_line);
    std::string del1;

    // Branje �tevila vozli��
    int n_vozlisc; // �tevilo vozli��
    iss >> del1; // Presko�i prvi del
    iss >> n_vozlisc; // Preberemo �tevilo vozli��

    // Branje vozli��
    for (int i = 0; i < n_vozlisc; i++) {
        std::string s;
        std::getline(file, s); // Preberemo vrstico
        std::replace(s.begin(), s.end(), ';', ' '); // Zamenjaj ';' z ' '
        std::replace(s.begin(), s.end(), ',', ' '); // Zamenjaj ',' z ' '
        std::istringstream iss(s);
        int node_id;
        double x;
        double y;
        iss >> node_id >> x >> y; // Preberemo ID vozli��a in njegove koordinate

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
    int n_celice = std::stoi(cells_string[1]); // Preberemo �tevilo celic

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
        iss >> cell_id >> node1_id >> node2_id >> node3_id >> node4_id; // Preberemo ID celice in njene vozli��a
        vector<int> celica;

        // Dodamo vozli��a celice v vektor
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

    int n_pogoji = 0; // �tevilo robnih pogojev
    std::string del2;
    std::string del3;

    iss3 >> del2; // Presko�i prvi del
    iss3 >> del3; // Presko�i drugi del
    iss3 >> n_pogoji; // Preberemo �tevilo robnih pogojev

    // Branje robnih pogojev
    for (int i = 0; i < n_pogoji; i++) {
        std::string s;
        string tip_pogoja;
        std::getline(file, s); // Preberemo vrstico
        std::istringstream iss4(s);
        std::string npd1;
        std::string npd2;
        iss4 >> npd1; // Presko�i prvi del
        iss4 >> npd2; // Presko�i drugi del
        iss4 >> tip_pogoja; // Preberemo tip robnega pogoja

        // Branje datoteke, �e imamo opravka z znano temperaturo na robu
        if (tip_pogoja == "temperatura") {
            tipi_robnih_pogojev.push_back(0); // Dodamo tip robnega pogoja
            std::string rp_temp;
            std::getline(file, rp_temp); // Preberemo temperaturo
            std::istringstream issrp1(rp_temp);
            double temperatura;
            std::string del4;

            issrp1 >> del4; // Presko�i del
            issrp1 >> temperatura; // Preberemo temperaturo

            vrednosti_robnih_pogojev.push_back(temperatura); // Dodamo temperaturo v vektor
            vrednosti_prestopa_toplote.push_back(-1); // Dodamo placeholder za toplotni prestop
        }
        // Branje datoteke, �e imamo opravka z toplotnim tokom na robu
        else if (tip_pogoja == "toplotni") {
            tipi_robnih_pogojev.push_back(1); // Dodamo tip robnega pogoja
            std::string rp_q;
            std::getline(file, rp_q); // Preberemo toplotni tok
            std::istringstream issrp2(rp_q);
            double toplotni_tok = 0;

            std::string del4;
            issrp2 >> del4; // Presko�i del
            issrp2 >> toplotni_tok; // Preberemo toplotni tok

            vrednosti_robnih_pogojev.push_back(toplotni_tok); // Dodamo toplotni tok v vektor
        }

        int st_vozlisc_v_rp = 0; // �tevilo vozli�� v robnem pogoju
        std::string st_vozlisc_rp;
        std::getline(file, st_vozlisc_rp); // Preberemo �tevilo vozli�� v robnem pogoju
        std::istringstream strp(st_vozlisc_rp);
        strp >> st_vozlisc_v_rp; // Preberemo �tevilo vozli��

        vector<int> vozlisca_v_robnem_pogoju; // Vektor za shranjevanje vozli�� robnega pogoja

        for (int i = 0; i < st_vozlisc_v_rp; i++) {
            std::string id_voz_rp;
            std::getline(file, id_voz_rp); // Preberemo ID vozli��a
            std::istringstream idvozrp(id_voz_rp);
            int id_vozlisca = 0;
            idvozrp >> id_vozlisca; // Preberemo ID vozli��a

            vozlisca_v_robnem_pogoju.push_back(id_vozlisca); // Dodamo vozli��e v vektor
        }

        vozlisca_robnih_pogojev.push_back(vozlisca_v_robnem_pogoju); // Dodamo vozli��a robnega pogoja v vektor

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

    std::cout << "Vozli��a robnih pogojev: " << std::endl;
    print2DVector(vozlisca_robnih_pogojev);

    // Konec branja datoteke

    // Sledi blok kode, ki preveri sosednjost vozli�� in to zapi�e v obliki vektorjev
    double deltaX = 1.0; // Razdalja med vozli��i v X smeri
    double deltaY = 1.0; // Razdalja med vozli��i v Y smeri
    double k = 24.0; // Koeficient toplotne prevodnosti

    vector<vector<int>> sosednja_vozlisca; // Vektor za shranjevanje sosednjih vozli��

    for (int node_id = 0; node_id < n_vozlisc; node_id++) {
        vector<int> node_i_neighbours = { -1,-1,-1,-1 }; // Inicializacija sosednjih vozli��

        for (int nd = 0; nd < n_celice; nd++) {
            vector<int> trenutna_celica = celice[nd]; // Pridobimo trenutne celice

            int vozlisce1 = trenutna_celica[0];
            int vozlisce2 = trenutna_celica[1];
            int vozlisce3 = trenutna_celica[2];
            int vozlisce4 = trenutna_celica[3];

            // Preverimo, ali je vozli��e del trenutne celice
            if (node_id == vozlisce1 || node_id == vozlisce2 || node_id == vozlisce3 || node_id == vozlisce4) {
                for (int vozl = 0; vozl < 4; vozl++) {
                    int sosednje_vozlisce = trenutna_celica[vozl]; // Pridobimo sosednje vozli��e

                    int pozicija = 0; // Inicializacija pozicije

                    // Preverimo, ali je sosednje vozli��e razli�ne od trenutnega vozli��a
                    if (sosednje_vozlisce != node_id) {
                        double x_obravnavano_vozl = X[node_id]; // X koordinata trenutnega vozli��a
                        double y_obravnavano_vozl = Y[node_id]; // Y koordinata trenutnega vozli��a

                        double x_sosed = X[sosednje_vozlisce]; // X koordinata sosednjega vozli��a
                        double y_sosed = Y[sosednje_vozlisce]; // Y koordinata sosednjega vozli��a

                        // Dolo�imo pozicijo sosednjega vozli��a glede na trenutnega
                        if ((x_obravnavano_vozl - x_sosed) < 1e-9 && (x_obravnavano_vozl - x_sosed) > -(1e-9)) {
                            if ((y_obravnavano_vozl - y_sosed) > 0.0) {
                                pozicija = 1; // Sosednje vozli��e je zgoraj
                            }
                            else {
                                pozicija = 3; // Sosednje vozli��e je spodaj
                            }
                        }
                        else if ((y_obravnavano_vozl - y_sosed) < 1e-9 && (y_obravnavano_vozl - y_sosed) > -(1e-9)) {
                            if ((x_obravnavano_vozl - x_sosed) > 0.0) {
                                pozicija = 0; // Sosednje vozli��e je levo
                            }
                            else {
                                pozicija = 2; // Sosednje vozli��e je desno
                            }
                        }
                        else {
                            pozicija = -1; // Sosednje vozli��e ni v pravilni poziciji
                        }

                        // Dodamo sosednje vozli��e v ustrezno pozicijo
                        if (pozicija != -1) {
                            node_i_neighbours[pozicija] = sosednje_vozlisce;
                        }
                    }
                }
            }
        }
        sosednja_vozlisca.push_back(node_i_neighbours); // Dodamo sosednja vozli��a v vektor
    }

    std::cout << "Vektor sosednjih vozli��: " << std::endl;
    print2DVector(sosednja_vozlisca);

    int n = n_vozlisc; // �tevilo vozli��

    // Ustvarimo matriko koeficientov A in vektor neznank b
    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0)); // Matriko A inicializiramo z ni�lami
    std::vector<double> b(n, 0.0); // Vektor b inicializiramo z ni�lami

    // Izpolnimo matriko A in vektor b
    for (int node_id = 0; node_id < n_vozlisc; node_id++) {
        std::vector<int> sosedi = sosednja_vozlisca[node_id]; // Pridobimo sosednja vozli��a
        int levi_sosed = sosedi[0];
        int spodnji_sosed = sosedi[1];
        int desni_sosed = sosedi[2];
        int zgornji_sosed = sosedi[3];

        // Preverimo, ali je vozli��e v notranjosti
        if (levi_sosed != -1 && spodnji_sosed != -1 && desni_sosed != -1 && zgornji_sosed != -1) {
            A[node_id][levi_sosed] = 1.0; // Dolo�imo koeficiente za sosednje vozli��a
            A[node_id][spodnji_sosed] = 1.0;
            A[node_id][desni_sosed] = 1.0;
            A[node_id][zgornji_sosed] = 1.0;
            A[node_id][node_id] = -4.0; // Dolo�imo koeficient za trenutno vozli��e
        }
        else { // Vozli��a na robu
            int tip_robnega_pogoja = 0; // Inicializiramo tip robnega pogoja
            double vrednost = 0; // Inicializiramo vrednost

            // Preverimo robne pogoje
            for (int robni_pogoj_id = 0; robni_pogoj_id < 5; robni_pogoj_id++) {
                std::vector<int> vozlisca_v_trenutnem_rp = vozlisca_robnih_pogojev[robni_pogoj_id];

                for (int id_vozlisce_rp = 0; id_vozlisce_rp < vozlisca_v_trenutnem_rp.size(); id_vozlisce_rp++) {
                    int vozlisce_v_trenutnem_rp = vozlisca_v_trenutnem_rp[id_vozlisce_rp];

                    if (node_id == vozlisce_v_trenutnem_rp) {
                        tip_robnega_pogoja = tipi_robnih_pogojev[robni_pogoj_id]; // Dolo�imo tip robnega pogoja
                        vrednost = vrednosti_robnih_pogojev[robni_pogoj_id]; // Dolo�imo vrednost robnega pogoja
                    }
                }
            }

            // Obdelava Dirichletovega robnega pogoja
            if (tip_robnega_pogoja == 0) {
                A[node_id][node_id] = 1.0; // Dolo�imo koeficient za robni pogoj
                b[node_id] = vrednost; // Dolo�imo vrednost vektorja b
            }
            // Obdelava Neumannovega robnega pogoja
            else if (tip_robnega_pogoja == 1) {
                int stevilo_sosedov = 0; // �tevilo sosedov

                for (int st = 0; st < 4; st++) {
                    if (sosedi[st] != -1) {
                        stevilo_sosedov++; // Pre�tejemo sosednje vozli��a
                    }
                }

                // Dolo�imo koeficiente za robne pogoje
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
        fileIDA << A[i][0]; // Zapi�i prvi element vrstice

        for (int j = 1; j < A[i].size(); j++) {
            fileIDA << ", " << A[i][j]; // Zapi�i preostale elemente vrstice
        }

        fileIDA << std::endl; // Prehod na novo vrstico
    }
    fileIDA.close(); // Zapri datoteko

    std::ofstream fileIDb;
    fileIDb.open("b.csv"); // Odpri datoteko za shranjevanje vektorja b

    for (int i = 0; i < b.size(); i++) {
        fileIDb << b[i]; // Zapi�i element vektorja b
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

    int st_iter = 1000; // Najve�je �tevilo iteracij

    // Nastavitev �tevila niti
    omp_set_num_threads(16); // Nastavi �tevilo niti za delovanje
    omp_set_dynamic(1); // Omogo�i dinami�no dodeljevanje niti

    // Za�etek merjenja �asa sistema ena�b
    auto start_time = std::chrono::high_resolution_clock::now();

    // Gauss-Seidel metoda s SOR
    for (int iter = 0; iter < st_iter && !converged; iter++) {
        converged = true;

    //Successive Over-Relaxation (SOR) je iterativna metoda za re�evanje linearnih sistemov, 
    //ki izbolj�a hitrost konvergence Gauss-Seidelove metode z uporabo relaksacijskega faktorja 
    //omega, ki prilagaja posodobitve vrednosti, da se pospe�i konvergenca in zmanj�a �tevilo potrebnih iteracij

        // Rde�i vozli (parni indeksi)
#pragma omp parallel for schedule(dynamic) shared(A, b, T, n) reduction(&&: converged)
        for (int jj = 0; jj < n; jj++) {
            if (jj % 2 == 0) { // Rde�i vozli
                double d = b[jj];
                for (int ii = 0; ii < n; ii++) {
                    if (jj != ii) {
                        d -= A[jj][ii] * T[ii]; // Izra�unaj novo vrednost
                    }
                }
                double T_new = d / A[jj][jj]; // Izra�unaj novo temperaturo
                double T_updated = omega * T_new + (1 - omega) * T[jj]; // Uporabi SOR
                if (abs(T_updated - T_old[jj]) > max_error) {
                    converged = false;
                }
                T[jj] = T_updated;
            }
        }

        // �rni vozli (neparni indeksi)
#pragma omp parallel for schedule(dynamic) shared(A, b, T, n) reduction(&&: converged)
        for (int jj = 0; jj < n; jj++) {
            if (jj % 2 != 0) { // �rni vozli
                double d = b[jj];
                for (int ii = 0; ii < n; ii++) {
                    if (jj != ii) {
                        d -= A[jj][ii] * T[ii]; // Izra�unaj novo vrednost
                    }
                }
                double T_new = d / A[jj][jj]; // Izra�unaj novo temperaturo
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

    // Konec merjenja �asa
    auto end_time = std::chrono::high_resolution_clock::now();
    
    // Izpis trajanja
    std::chrono::duration<double> time_duration = end_time - start_time;
    std::cout << "�as Gauss-Seidl metode: " << time_duration.count() << " sekund" << std::endl;

    // Zapis rezultatov v .vtk datoteki
    std::ofstream fileID;
    fileID.open("rezultat.vtk"); // Odpri datoteko za shranjevanje rezultatov

    fileID << "# vtk DataFile Version 3.0\n"; // Zapis verzije VTK datoteke
    fileID << "Mesh_1\n"; // Ime mre�e
    fileID << "ASCII\n"; // Oblika datoteke
    fileID << "DATASET UNSTRUCTURED_GRID\n"; // Tip podatkov
    fileID << "POINTS " << n_vozlisc << " float\n"; // Zapis �tevila vozli��
    for (int koordinata_id = 0; koordinata_id < n_vozlisc; koordinata_id++) {
        fileID << X[koordinata_id] << " " << Y[koordinata_id] << " 0\n"; // Zapis koordinat vozli��
    }
    fileID << "\n";
    fileID << "CELLS " << n_celice << " " << n_celice * 5 << "\n"; // Zapis �tevila celic
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
    fileID << "POINT_DATA " << n_vozlisc << "\n"; // Zapis podatkov o vozli��ih
    fileID << "SCALARS Temperature float 1\n"; // Zapis temperature kot skalarne vrednosti
    fileID << "LOOKUP_TABLE default\n"; // Dolo�itev privzete tabele
    for (int koordinata_id = 0; koordinata_id < n_vozlisc; koordinata_id++) {
        fileID << T[koordinata_id] << "\n"; // Zapis temperature za vsako vozli��e
    }

    fileID.close(); // Zapri datoteko

    return 0; // Konec programa
}