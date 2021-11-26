%% Bearbeitungsbogen
clear all
close all
clc

%% Parameter
l = 0.2;                                                           % Laenge l
g = 9.81;                                                          % Erdbeschleunigung
omega = 0.5*sqrt(g/l);                                             % Rotationsgeschwindigkeit
AnzahlSchritte = 10^4;                                             % Anzahl der Zeitschritte
r = 2*pi-0                                                         % Range von 0 - 2*pi
h = r/AnzahlSchritte;                                              % Zeitschrittweite

x_0 = [pi/6;0];                                                    % Anfangsbedingungen in Zustandsvektor
t = 0:h:r;                                                         % Zeitvektor

x_Eu_expl(:,1) = x_0;                                              % Anfangsbedingungen fuer Euler explizit
x_Eu_impl(:,1) = x_0;                                              % Anfangsbedingungen fuer Euler implizit
x_RuKu(:,1) = x_0;                                                 % Anfangsbedingungen fuer Runge-Kutta-Verfahren

%% Aufgabe 1: Systemmatrix
SystMatr = [0, 1; omega^2-(g/l), 0];                               % Systemmatrix

%% Aufgabe 2: exakte Loesung
x_exakt = [(pi/12) * (exp(sqrt(omega^2-g/l)*t) + exp(-sqrt(omega^2-g/l)*t)) ; 
           (sqrt(omega^2-g/l)* (pi/12) * exp(sqrt(omega^2-g/l)*t)) + (-sqrt(omega^2-g/l)* (pi/12) * exp(-sqrt(omega^2-g/l)*t))];

for n = 'ausfuellen'
    
    %% Aufgabe 3: Euler explizit
    x_Eu_expl(:,n+1) = 'ausfuellen';
    
    %% Aufgabe 4: Euler implizit
    x_Eu_impl(:,n+1) = 'ausfuellen';
    
    %% Aufgabe 5: Runge-Kutta-Verfahren
    k1(:,n) = 'ausfuellen';
    k2(:,n) = 'ausfuellen';
    k3(:,n) = 'ausfuellen';
    k4(:,n) = 'ausfuellen';
    x_RuKu(:,n+1) = 'ausfuellen';
    
end

%% Aufgabe 6
l_Eu_expl = 'ausfuellen';
l_Eu_impl = 'ausfuellen';
l_RuKu = 'ausfuellen';

e_Eu_expl = 'ausfuellen';
e_Eu_impl = 'ausfuellen';
e_RuKu = 'ausfuellen';

%% Aufgabe 7: Plots
'ausfuellen';

%% Aufgabe 8: Plots
'ausfuellen';