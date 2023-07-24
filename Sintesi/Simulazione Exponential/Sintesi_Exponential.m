clc
clear
close

%% Quadrilatero Articolato Generatore di Funzione
%  sintesi del terzo ordine con funzione di terzo grado

f = 40;             % lunghezza telaio in cm
x_m = 380;          % punto di precisione
size = 500;         % dimensione array

delta_theta = 40;   % escursione massima angolo input
delta_phi = 20;     % escursione massina angolo output
theta = linspace(0, delta_theta, size);
phi = linspace(0, delta_phi, size);

figure(1);
plot(theta, phi);
title("Angolo di output in funzione dell'angolo di input");
xlabel("θ [°]");
ylabel("ϕ [°]");
legend("ϕ = f(θ)")
grid on;
set(gca, "xlim", [0 40]);
set(gca, "ylim", [0 20]);

x = [0, 300, 500];
y = [0, 20, 150];
funct = fit(x',y','exp1');

%p = [funct.p1, funct.p2, funct.p3, funct.p4];
p = [funct.a, funct.b];

% grafico funzione da generare tramite sintesi
figure(2);
plot(funct, x, y);
title("Funzione da generare tramite sintesi")
xlabel("L [mm]");
ylabel("δ [mm]");
legend("δ_i = f(L_i)","δ = f(L)");
grid on;
set(gca, "xlim", [0 500]);
set(gca, "ylim", [0 150]);

x_data = linspace(0, 500, size);
y_data = zeros(size,1);
for i=1:size
    y_data(i) = p(1)*exp(p(2)*x_data(i));
end

% calcolo fattori di scala
K_theta = delta_theta / x_data(size);
K_phi = delta_phi / (y_data(size) - y_data(1));

% calcolo R, R', R" per sintesi del terzo ordineS
R = K_phi / K_theta * (p(1)*p(2)*exp(p(2)*x_data(x_m)));
R_1 = K_phi / (K_theta)^2 * (p(1)*p(2)*p(2)*exp(p(2)*x_data(x_m)));
R_2 = K_phi / (K_theta)^3 * (p(1)*p(2)*p(2)*p(2)*exp(p(2)*x_data(x_m)));

% lunghezza segmento e
e = -R*f/(R-1);

% diametro della circonferenza di Carter-Hall
d_c = 3*f*(R^2 * (1-R)^2 + R_1^2)/((1-R)*(R*(1-R)^3 + 2*R^2*(1-R)^2 + 3*R_1^2 + R_2*(1-R)));

% angoli per la procedura di sintesi del terzo ordine usando la
% circonferenza di Carter-Hall
psi = atan2(R*(1-R), R_1);
lambda = 10*pi/180;

% coordinate del centro di istantanea rotazione del moto relativo
P24x = d_c*(cos(lambda+psi))^2;
P24y = d_c/2*(sin(2*(lambda+psi)));

%% Grafici Sintesi
% dati ottenuti tramite simulazione su Solidworks
% problemi:
% -scala diversa
% -segni diversi

load("Workspace_Exponential");

%% Lambda = 5
% rescale dei dati
SpostamentoAngolareForcellone1 = rescale(SpostamentoAngolare1deg, 0, max(SpostamentoAngolare1deg) - min(SpostamentoAngolare1deg));
SpostamentoAngolareCedente1 = rescale(SpostamentoAngolare2deg, 0, max(SpostamentoAngolare2deg) - min(SpostamentoAngolare2deg));
CompressioneMolla1 = rescale(Distanza7mm, 0, max(Distanza7mm) - min(Distanza7mm));
SpostamentoForcellone1 = rescale(Distanza5mm, 0, max(Distanza5mm) - min(Distanza5mm));
AngoloTrasmissione1 = SpostamentoAngolare3deg;

% compressione molla negativa -> compressione molla positiva
CompressioneMolla1 = - CompressioneMolla1 + max(CompressioneMolla1);

%% Lambda = 7.5
% rescale dei dati
SpostamentoAngolareForcellone2 = rescale(SpostamentoAngolare1deg1, 0, max(SpostamentoAngolare1deg1) - min(SpostamentoAngolare1deg1));
SpostamentoAngolareCedente2 = rescale(SpostamentoAngolare2deg1, 0, max(SpostamentoAngolare2deg1) - min(SpostamentoAngolare2deg1));
CompressioneMolla2 = rescale(Distanza7mm1, 0, max(Distanza7mm1) - min(Distanza7mm1));
SpostamentoForcellone2 = rescale(Distanza5mm1, 0, max(Distanza5mm1) - min(Distanza5mm1));
AngoloTrasmissione2 = SpostamentoAngolare3deg1;

% compressione molla negativa -> compressione molla positiva
CompressioneMolla2 = - CompressioneMolla2 + max(CompressioneMolla2);

%% Lambda = 10
% rescale dei dati
SpostamentoAngolareForcellone3 = rescale(SpostamentoAngolare1deg2, 0, max(SpostamentoAngolare1deg2) - min(SpostamentoAngolare1deg2));
SpostamentoAngolareCedente3 = rescale(SpostamentoAngolare2deg2, 0, max(SpostamentoAngolare2deg2) - min(SpostamentoAngolare2deg2));
CompressioneMolla3 = rescale(Distanza7mm2, 0, max(Distanza7mm2) - min(Distanza7mm2));
SpostamentoForcellone3 = rescale(Distanza5mm2, 0, max(Distanza5mm2) - min(Distanza5mm2));
AngoloTrasmissione3 = SpostamentoAngolare3deg2;

% compressione molla negativa -> compressione molla positiva
CompressioneMolla3 = - CompressioneMolla3 + max(CompressioneMolla3);

%%

figure(3);
plot(theta, phi, SpostamentoAngolareForcellone1, SpostamentoAngolareCedente1, SpostamentoAngolareForcellone2, SpostamentoAngolareCedente2, SpostamentoAngolareForcellone3, SpostamentoAngolareCedente3);
title("Angolo di output in funzione dell'angolo di input");
xlabel("θ [°]");
ylabel("ϕ [°]");
legend("ϕ = f(θ)", "λ = 5", "λ = 7.5", "λ = 10")
grid on;
set(gca, "xlim", [0 40]);
set(gca, "ylim", [0 20]);

figure(4);
plot(x_data, y_data, SpostamentoForcellone1, CompressioneMolla1, SpostamentoForcellone2, CompressioneMolla2, SpostamentoForcellone3, CompressioneMolla3);
title("Funzione generata tramite sintesi");
xlabel("L [mm]");
ylabel("δ [mm]");
legend("δ = f(L)", "λ = 5", "λ = 7.5", "λ = 10");
grid on;
set(gca, "xlim", [0 500]);
set(gca, "ylim", [0 200]);

VM1 = zeros(length(SpostamentoForcellone1),1);
VM2 = zeros(length(SpostamentoForcellone2),1);
VM3 = zeros(length(SpostamentoForcellone3),1);
for i=1:length(SpostamentoForcellone1)
    VM1(i) = CompressioneMolla1(i)/SpostamentoForcellone1(i);
    VM2(i) = CompressioneMolla2(i)/SpostamentoForcellone2(i);
    VM3(i) = CompressioneMolla3(i)/SpostamentoForcellone3(i);
end

figure(5);
plot(SpostamentoAngolareForcellone1, VM1, SpostamentoAngolareForcellone2, VM2, SpostamentoAngolareForcellone3, VM3);
title("Vantaggio meccanico variabile");
xlabel("θ [°]");
ylabel("VM [ ]");
legend("λ = 5", "λ = 7.5", "λ = 10");
grid on;

figure(6);
plot(SpostamentoAngolareForcellone1, AngoloTrasmissione1, SpostamentoAngolareForcellone2, AngoloTrasmissione2, SpostamentoAngolareForcellone3, AngoloTrasmissione3);
title("Angolo di trasmissione variabile");
xlabel("θ [°]");
ylabel("μ - 90 [°]");
legend("λ = 5", "λ = 7.5", "λ = 10");
grid on;