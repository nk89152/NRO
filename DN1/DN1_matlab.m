% NALOGA 1

% Nastavimo prikaz večjega števila decimalk
format long;

% Odpremo datoteko v načinu za branje podatkov
fid = fopen('naloga1_1.txt', 'r');

% Prva vrstica ne vsebuje vrednosti, zato jo preskočimo
fgets(fid);

% Preberemo informacijo o številu vrstic in številu podatkov v vsaki vrstici
dimenzije = fscanf(fid, 'stevilo preostalih vrstic: %d; stevilo podatkov v vrstici: %d', [1 2]);
stVrstic = dimenzije(1);  % Število preostalih vrstic
stPodatkov = dimenzije(2); % Število podatkov v posamezni vrstici

% Preberemo podatke v obliki matrike 
podatki = fscanf(fid, '%f', [stPodatkov, stVrstic]);
podatki = podatki'; % Matriko preberemo v obliki vektorja vrstic za lažje delo

% Datoteko zapremo
fclose(fid);

% Pretvorimo matriko 'podatki' v enodimenzionalni vektor 't' za nadaljnjo obdelavo
t = podatki(:);


% NALOGA 2

% Odpremo drugo datoteko za branje podatkov
fid = fopen('naloga1_2.txt', 'r');

% Preberemo prvo vrstico, ki določa skupno število vrednosti v drugi datoteki
stVrednosti = fscanf(fid, 'stevilo_podatkov_P: %d', 1);

% Ustvarimo prazen vektor 'P', kamor bomo shranjevali prebrane vrednosti
P = zeros(stVrednosti, 1);

% Vrednosti shranimo v vektor P z uporabo zanke for
for i = 1:stVrednosti
    P(i) = fscanf(fid, '%f', 1);  % Preberemo naslednjo vrednost kot decimalno število
end

% Datoteko zapremo
fclose(fid);

% Izrišemo graf odvisnosti P od t
plot(t, P, 'b-', 'LineWidth', 1); % Izrišemo graf z modro črto in širino črte 2
xlabel('t [s]'); % Oznaka x-osi predstavlja čas v sekundah
ylabel('P [W]'); % Oznaka y-osi predstavlja moč v vatih
title('Graf funkcije P(t)'); % Naslov grafa prikazuje odvisnost moči od časa
grid on; % Vključimo mrežo za boljšo preglednost grafa


% NALOGA 3

% Določimo najmanjšo in največjo vrednost vektorja 't' za izračun intervalov
tmin = min(t);
tmax = max(t);

% Izračunamo število intervalov, ki jih potrebujemo za trapezno integracijo
n = length(t) - 1;

% Izračunamo širino posameznega intervala
dx = (tmax - tmin) / n;

% Nastavimo začetno vrednost za seštevek trapezov
integral_vrednost = 0;

% Z zanko for seštejemo površino trapezov, s čimer dobimo približek integrala
for i = 1:n
    integral_vrednost = integral_vrednost + (dx / 2) * (P(i) + P(i + 1));
end

% Prikažemo vrednost ročno izračunanega integrala
disp(['Ročno izračunan integral:',num2str(integral_vrednost, '%.10f')]);

% Izračunamo in prikažemo vrednost izračunanega integrala s funkcijo trapz
trapz_integral = trapz(t, P);
disp(['Izračunan integral s funkcijo trapz: ', num2str(trapz_integral, '%.10f')]);
