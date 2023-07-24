clc
clear
close

%% Quadrilatero Articolato Generatore di Funzione
%  sintesi del terzo ordine con funzione di terzo grado

f = 40;             % lunghezza telaio in cm
x_m = 1;         % punto di lavoro

delta_theta = 40;   % escursione massima angolo input
delta_phi = 20;     % escursione massina angolo output
theta = linspace(0, delta_theta, 100);
phi = linspace(0, delta_phi, 100);
[funct_angles, error_angles] = fit(theta', phi', fittype('poly1'));

figure(1);
plot(funct_angles, theta, phi);
title("Angolo di output in funzione dell'angolo di input");
xlabel("θ [°]");
ylabel("ϕ [°]");
legend("Data","ϕ = f(θ)")
grid on;
% set(gca, "xlim", [0 35]);
% set(gca, "ylim", [0 20]);

% punti da generare tramite polinomio di terzo grado
x = [0, 100, 200, 400];
y = [0, 5, 25, 125];
[funct, error] = fit(x', y', fittype('poly3'));

disp(funct)

% funzione da generare tramite sintesi
figure(2);
plot(funct, x, y);
title("Funzione da generare tramite sintesi")
xlabel("L [cm]");
ylabel("δ [cm]");
legend("Data","δ = f(L)");
grid on;
% set(gca, "xlim", [0 38]);
% set(gca, "ylim", [0 10]);

% calcolo fattori di scala
K_theta = delta_theta / (max(x) - min(x));
K_phi = delta_phi / (max(y) - min(y));

% calcolo R, R', R" per sintesi del terzo ordine
x_data = linspace(0,max(x),1000);
R = K_phi / K_theta * (funct.p1*3*x_data(x_m)^2 + funct.p2*2*x_data(x_m) + funct.p3);
R_1 = K_phi / (K_theta)^2 * (funct.p1*3*2*x_data(x_m) + funct.p2*2);
R_2 = K_phi / (K_theta)^3 * (funct.p1*3*2);

% a = funct.a;
% b = funct.b;
% R = K_phi / K_theta * (a*b*exp(b*x_data(x_m)));
% R_1 = K_phi / (K_theta)^2 * (a*b*b*exp(b*x_data(x_m)));
% R_2 = K_phi / (K_theta)^3 * (a*b*b*b*exp(b*x_data(x_m)));

% lunghezza segmento e
e = -R*f/(R-1);

% diametro della circonferenza di Carter-Hall
d_c = 3*f*(R^2 * (1-R)^2 + R_1^2)/((1-R)*(R*(1-R)^3 + 2*R^2*(1-R)^2 + 3*R_1^2 + R_2*(1-R)));

% angoli per la procedura di sintesi del terzo ordine usando la
% circonferenza di Carter-Hall
psi = atan2(R*(1-R), R_1);
lambda = 30*pi/180;

% coordinate del centro di istantanea rotazione del moto relativo
P24x = d_c*(cos(lambda+psi))^2;
P24y = d_c/2*(sin(2*(lambda+psi)));

%% Grafici Sintesi
% dati ottenuti tramite simulazione su Solidworks
% problemi:
% -scala diversa
% -segni diversi
% -unità di misura diverse

% rescale dei dati
SpostamentoAngolareForcellone = rescale(SpostamentoAngolare1deg, 0, max(SpostamentoAngolare1deg) - min(SpostamentoAngolare1deg));
SpostamentoAngolareCedente = rescale(SpostamentoAngolare2deg, 0, max(SpostamentoAngolare2deg) - min(SpostamentoAngolare2deg));
CompressioneMolla = rescale(Distanza6mm, 0, max(Distanza6mm) - min(Distanza6mm));
SpostamentoForcellone = rescale(Distanza5mm, 0, max(Distanza5mm) - min(Distanza5mm));

% compressione molla negativa -> compressione molla positiva
CompressioneMolla = - CompressioneMolla + max(CompressioneMolla);

% conversione da [mm] a [cm]
CompressioneMolla = CompressioneMolla/10;
SpostamentoForcellone = SpostamentoForcellone/10;

figure(3);
plot(funct_angles, SpostamentoAngolareForcellone, SpostamentoAngolareCedente);
title("Angolo di output in funzione dell'angolo di input");
xlabel("θ [°]");
ylabel("ϕ [°]");
legend("Fitted", "Generated");
grid on;
%set(gca, "xlim", [0 35]);
%set(gca, "ylim", [0 20]);

VM = zeros(length(SpostamentoForcellone),1);
for i=1:length(SpostamentoForcellone)
    VM(i) = SpostamentoForcellone(i)/CompressioneMolla(i);
end

figure(4);
plot(funct, SpostamentoForcellone, CompressioneMolla);
title("Funzione generata tramite sintesi");
xlabel("L [cm]");
ylabel("δ [cm]");
legend("Fitted", "Generated");
grid on;
%set(gca, "xlim", [0 38]);
%set(gca, "ylim", [0 10]);

figure(5);
plot(VM);
title("Vantaggio meccanico variabile nel tempo");
xlabel("L/δ [cm]");
ylabel("VM");
legend("VM = L/δ");
grid on;

%% Punto di lavoro

for x_m = 1:1:500
    % calcolo R, R', R" per sintesi del terzo ordine
    R = K_phi / K_theta * (funct.a*funct.b*exp(funct.b*x_data(x_m)));
    R_1 = K_phi / (K_theta)^2 * (funct.a*funct.b*funct.b*exp(funct.b*x_data(x_m)));
    R_2 = K_phi / (K_theta)^3 * (funct.a*funct.b*funct.b*funct.b*exp(funct.b*x_data(x_m)));
    
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
    
    if (e > 0 && d_c > 0 && P24y > 0)
        disp(x_m);
        disp(e);
        disp(d_c);
        disp(P24x);
        disp(P24y);
    end
end
%% Valori di R, R', R" 

R_vector = zeros(length(x_data),1);

for x_m = 1:1:length(x_data)
    R_vector(x_m) = K_phi / K_theta * (a*3*x_data(x_m)^2 + b*2*x_data(x_m) + c);
    R_1vector(x_m) = K_phi / (K_theta)^2 * (a*3*2*x_data(x_m) + b*2);
    R_2vector(x_m) = K_phi / (K_theta)^3 * (a*3*2);
end

figure(5)
plot(R_vector);

figure(6);
plot(R_1vector);

figure(7);
plot(R_2vector);