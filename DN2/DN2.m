clc;
clear;

% METODA: "scatteredInterpolant"

tic;
% Odpremo datoteko 'vozlisca_temperature_dn2.txt' za branje
fileID = fopen('vozlisca_temperature_dn2.txt', 'r');

% Preskočimo prve štiri vrstice z informacijami in opisom
for vrstica = 1:4
    fgetl(fileID);
end

% Uporabimo textscan za branje podatkov (x, y, temperatura), kjer vejica služi kot ločilo
vsebina = textscan(fileID, '%f %f %f', 'Delimiter', ',', 'CollectOutput', true);
fclose(fileID);

% Podatke shranimo kot matriko za enostaven dostop
podatki = vsebina{1};

% Razdelimo podatke na posamezne spremenljivke: x, y, temperatura
x_koordinate = podatki(:, 1); % Prvi stolpec: x-koordinate
y_koordinate = podatki(:, 2); % Drugi stolpec: y-koordinate
temperature = podatki(:, 3); % Tretji stolpec: temperature

% Kreiramo interpolacijsko funkcijo
interp_funkcija = scatteredInterpolant(x_koordinate, y_koordinate, temperature, 'linear', 'none');

% Definiramo točko, kjer želimo izračunati temperaturo
x_tocka1 = 0.403;
y_tocka1 = 0.503;

% Izračunamo temperaturo v izbrani točki
T_tocka1 = interp_funkcija(x_tocka1, y_tocka1);

% Prikažemo rezultat
fprintf('Temperatura v točki (%.3f, %.3f) je približno %.3f.\n', x_tocka1, y_tocka1, T_tocka1);
toc;

% METODA: "griddedInterpolant"

tic;
% Rekonstruiramo mrežo iz podatkov
x_mreza = unique(x_koordinate); % Mreža v x-smeri
y_mreza = unique(y_koordinate); % Mreža v y-smeri

% Pretvorimo podatke o temperaturi v obliko mreže
T_mreza = reshape(temperature, numel(x_mreza), numel(y_mreza));

% Kreiramo interpolacijsko funkcijo
interp_mreza = griddedInterpolant({x_mreza, y_mreza}, T_mreza, 'linear', 'none');

% Definiramo drugo točko za izračun temperature
x_tocka2 = 0.403;
y_tocka2 = 0.503;

% Izračunamo temperaturo v tej točki
T_tocka2 = interp_mreza(x_tocka2, y_tocka2);

% Prikažemo rezultat
fprintf('Temperatura v točki (%.3f, %.3f) je približno %.3f.\n', x_tocka2, y_tocka2, T_tocka2);
toc;

% LASTNA METODA: Najbližji sosed

tic;
% Definiramo točko, kjer želimo določiti temperaturo
x_tocka3 = 0.403;
y_tocka3 = 0.503;

% Izračunamo razdalje med vsemi vozlišči in točko
razdalje = sqrt((x_koordinate - x_tocka3).^2 + (y_koordinate - y_tocka3).^2);

% Poiščemo indeks vozlišča z najmanjšo razdaljo
[~, indeks_min] = min(razdalje);

% Temperatura najbližjega vozlišča
T_najblizje = temperature(indeks_min);

% Prikažemo rezultat
fprintf('Temperatura v točki (%.3f, %.3f) po metodi najbližjega soseda je %.3f.\n', x_tocka3, y_tocka3, T_najblizje);
toc;

% Najvišja temperatura in njene koordinate

% Poiščemo najvišjo temperaturo in njen indeks
[naj_temp, indeks_max] = max(temperature);

% Koordinate točke z najvišjo temperaturo
x_naj = x_koordinate(indeks_max);
y_naj = y_koordinate(indeks_max);

% Izpišemo rezultat
fprintf('Najvišja temperatura je %.3f °C in se nahaja pri koordinatah (x = %.3f, y = %.3f).\n', naj_temp, x_naj, y_naj);
