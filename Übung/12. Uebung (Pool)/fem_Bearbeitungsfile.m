clear all;
close all;
clc

%% Parameter
L = 'ausf�llen';
F = 'ausf�llen';
E = 'ausf�llen';
AnzahlElemente = 'ausf�llen';
LElemente = (L/AnzahlElemente);       % L�nge eines Elemntes
Knoten = (0:L/AnzahlElemente:L);      % Knotenvektor
AElement = (0.001-0.0009*Knoten);     % Querschnittsverlauf

%% Aufgabe 2
%'ausf�llen'

%% Aufgabe 3
N = 'ausf�llen';

%% Aufgabe 4
B = 'ausf�llen';
C = 'ausf�llen';

%% Aufgabe 5
for b='ausf�llen'                     % Schleife f�r die Anzahl an gew�hlten Elementen
   K{b} = 'ausf�llen';
end

%% Aufgabe 6
KGesamt = zeros(length(Knoten),length(Knoten));
for b='ausf�llen'
    KGesamt(b:b+1,b:b+1) = KGesamt(b:b+1,b:b+1)+K{b};
end

%% Aufgabe 7
fb = 'ausf�llen';
ua = 'ausf�llen';

Kaa = 'ausf�llen';
Kab = 'ausf�llen';
Kba = 'ausf�llen';
Kbb = 'ausf�llen';

ub = 'ausf�llen';
fa = 'ausf�llen';

%% Aufgabe 8
%'ausf�llen'