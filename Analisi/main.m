clc
clear
close all

%% Data

% L1 = 1.1;
% L2 = 4.7;
% L3 = 2.6;
% L4 = 6.1;

L1 = 25.63;
L2 = 24.24;
L3 = 41.30;
L4 = 40;

theta_1min = 100*pi/180;
theta_1max = 150*pi/180;
theta_delta = 1*pi/180;
theta_i = [theta_1min; 1; -1; pi];

theta_1 = [theta_1min:theta_delta:theta_1max]';
size = length(theta_1);

epsilon = 1e-3;
max_iterations = 1e3;

[f, j] = FourBarLinkage(L1,L2,L3,L4);

%% Analysis of Angular Positions

theta_2 = zeros(size, 1);
theta_3 = zeros(size, 1);

i = 1;

while i <= size

    theta_i(1) = theta_1(i);
    [theta_est, iterations] = MethodNewtonRaphson(theta_i, epsilon, max_iterations, f, j);
    
    theta_2(i) = theta_est(2);
    theta_3(i) = theta_est(3);

    i = i + 1;

end

%% Analysis of Angular Velocities

syms w2 w3

w_1 = ones(size, 1);
w_2 = zeros(size, 1);
w_3 = zeros(size, 1);

i = 1;

while i <= size

    jacob = j([theta_1(i); theta_2(i); theta_3(i); pi]);
    b = [L1*sin(theta_1(i)); -L2*cos(theta_1(i))];
    
    w = solve(jacob*[w2; w3] == b*w_1(i));
    w_2(i) = double(w.w2);
    w_3(i) = double(w.w3);

    i = i + 1;

end

%% Analysis of Angular Accelerations

syms alpha2 alpha3

alpha_1 = zeros(size, 1);
alpha_2 = zeros(size, 1);
alpha_3 = zeros(size, 1);

i = 1;

while i <= size

    jacob = j([theta_1(i); theta_2(i); theta_3(i); pi]);
    b = [L1*alpha_1(i)*sin(theta_1(i))+L1*cos(theta_1(i))*w_1(i)^2 + L2*cos(theta_2(i))*w_2(i)^2 + L3*cos(theta_3(i))*w_3(i)^2;
         -L1*alpha_1(i)*cos(theta_1(i))+L1*sin(theta_1(i))*w_1(i)^2 + L2*sin(theta_2(i))*w_2(i)^2 + L3*sin(theta_3(i))*w_3(i)^2];
    
    alpha = solve(jacob*[alpha2; alpha3] == b);
    alpha_2(i) = double(alpha.alpha2);
    alpha_3(i) = double(alpha.alpha3);

    i = i + 1;

end

%% Analysis of Angular Jerks

syms jerk2 jerk3

jerk_1 = zeros(size, 1);
jerk_2 = zeros(size, 1);
jerk_3 = zeros(size, 1);

i = 1;
while i <= size

    jacob = j([theta_1(i); theta_2(i); theta_3(i); pi]);
    b = [+L1*jerk_1(i)*sin(theta_1(i)) + 3*L1*alpha_1(i)*w_1(i)*cos(theta_1(i)) - L1*w_1(i)^3*sin(theta_1(i)) + 3*L2*alpha_2(i)*w_2(i)*cos(theta_2(i)) - L2*w_2(i)^3*sin(theta_2(i)) + 3*L3*alpha_3(i)*w_3(i)*cos(theta_3(i)) - L3*w_3(i)^3*sin(theta_3(i));
         -L1*jerk_1(i)*cos(theta_1(i)) + 3*L1*alpha_1(i)*w_1(i)*sin(theta_1(i)) + L1*w_1(i)^3*cos(theta_1(i)) + 3*L2*alpha_2(i)*w_2(i)*sin(theta_2(i)) + L2*w_2(i)^3*cos(theta_2(i)) + 3*L3*alpha_3(i)*w_3(i)*sin(theta_3(i)) + L3*w_3(i)^3*cos(theta_3(i))];
    
    jerk = solve(jacob*[jerk2; jerk3] == b);
    jerk_2(i) = double(jerk.jerk2);
    jerk_3(i) = double(jerk.jerk3);

    i = i + 1;
end

%% Plots

% Plot of Angular Positions
figure(1);
subplot(1,2,1);
plot(theta_1, theta_2);
title("Angolo θ_2 in funzione dell'angolo θ_1");
xlabel("θ_1 [rad]");
ylabel("θ_2 [rad]");
legend("θ_2");
grid on;

subplot(1,2,2);
plot(theta_1, theta_3);
title("Angolo θ_3 in funzione dell'angolo θ_1");
xlabel("θ_1 [rad]");
ylabel("θ_3 [rad]");
legend("θ_3");
grid on;

% Plot of Angular Velocities
figure(2);
subplot(1,2,1);
plot(0:1:size-1, w_2);
title("Velocità Angolare ω_2 in funzione della Velocità Angolare ω_1 = 1 [rad/s]");
xlabel("ω_1 [rad/s]");
ylabel("ω_2 [rad/s]");
legend("ω_2");
grid on;

subplot(1,2,2);
plot(0:1:size-1, w_3);
title("Velocità Angolare ω_3 in funzione della Velocità Angolare ω_1 = 1 [rad/s]");
xlabel("ω_1 [rad/s]");
ylabel("ω_3 [rad/s]");
legend("ω_3");
grid on;

% Plot of Angular Accelerations
figure(3);
subplot(1,2,1);
plot(0:1:size-1, alpha_2);
title("Accelerazione Angolare α_2 in funzione dell'Accelerazione Angolare α_1 = 1 [rad/s^2]");
xlabel("α_1 [rad/s^2]");
ylabel("α_2 [rad/s^2]");
legend("α_2");
grid on;

subplot(1,2,2);
plot(0:1:size-1, alpha_3);
title("Accelerazione Angolare α_3 in funzione dell'Accelerazione Angolare α_1 = 0 [rad/s^2]");
xlabel("α_1 [rad/s^2]");
ylabel("α_3 [rad/s^2]");
legend("α_3");
grid on;

% Plot of Angular Jerks
figure(4);
subplot(1,2,1);
plot(0:1:size-1, jerk_2);
title("Jerk Angolare J_2 in funzione del Jerk Angolare J_1 = 0 [rad/s^3]");
xlabel("J_1 [rad/s^3]");
ylabel("J_2 [rad/s^3]");
legend("J_2");
grid on;

subplot(1,2,2);
plot(0:1:size-1, jerk_3);
title("Jerk Angolare J_3 in funzione del Jerk Angolare J_1 = 0 [rad/s^3]");
xlabel("J_1 [rad/s^3]");
ylabel("J_3 [rad/s^3]");
legend("J_3");
grid on;