% Naloži matriko A in vektor b iz datotek
A = readmatrix('A.csv'); % Naloži matriko A iz A.csv
b = readmatrix('b.csv'); % Naloži vektor b iz b.csv

% Pretvorba v sparse format za boljšo učinkovitost
A = sparse(A);
b = sparse(b);

% Validacija dimenzij vhodnih podatkov
if size(A, 1) ~= size(b, 1)
    error('Dimenzije matrike A in vektorja b niso skladne!');
end

% Parametri
tol = 1e-8;
max_iter = 1000;

% Začetek merjenja časa
start_time = tic;

% Reševanje sistema
[T, flag, relres, iter, resvec] = gmres(A, b, [], tol, max_iter);

% Merjenje časa
time_duration = toc(start_time);

% ===== IZPIS REZULTATOV =====
fprintf('\n=== REZULTATI REŠEVANJA SISTEMA ===\n');
fprintf('Velikost sistema: %d x %d\n', size(A,1), size(A,2));

% Izpis konvergence
if flag == 0
    fprintf('\nSTATUS: Uspešna konvergenca\n');
else
    fprintf('\nSTATUS: Metoda ni konvergirala (flag = %d)\n', flag);
end

% Izpis časa in iteracij
fprintf('\nČasovna statistika:\n');
fprintf('- Čas izvajanja: %.6f sekund\n', time_duration);
fprintf('- Število iteracij: %d\n', iter);
fprintf('- Relativni ostanek: %e\n', relres);

% Preveri točnost rešitve
residual = norm(A*T - b);
fprintf('\nTOČNOST REŠITVE:\n');
fprintf('- Norma ostanka ||Ax-b||: %e\n', residual);

% Pretvori rešitev v stolpični vektor (če ni že stolpični)
T = full(T(:)); % Pretvorimo sparse obliko v polno, če je potrebno

% Shranjevanje rešitve v CSV datoteko
csv_filename = 'resitev_T.csv'; % Ime datoteke
writematrix(T, csv_filename); % Shrani rešitev v CSV datoteko brez dodatnih znakov

fprintf('\nRešitev je bila uspešno shranjena v datoteko: %s\n', csv_filename);

