clear all;
close all;
clc

%% Parameter
L = 'ausfüllen';
F = 'ausfüllen';
E = 'ausfüllen';
AnzahlElemente = 'ausfüllen';
LElemente = (L/AnzahlElemente);       % Länge eines Elemntes
Knoten = (0:L/AnzahlElemente:L);      % Knotenvektor
AElement = (0.001-0.0009*Knoten);     % Querschnittsverlauf

%% Aufgabe 2
%'ausfüllen'

%% Aufgabe 3
N = 'ausfüllen';

%% Aufgabe 4
B = 'ausfüllen';
C = 'ausfüllen';

%% Aufgabe 5
for b='ausfüllen'                     % Schleife für die Anzahl an gewählten Elementen
   K{b} = 'ausfüllen';
end

%% Aufgabe 6
KGesamt = zeros(length(Knoten),length(Knoten));
for b='ausfüllen'
    KGesamt(b:b+1,b:b+1) = KGesamt(b:b+1,b:b+1)+K{b};
end

%% Aufgabe 7
fb = 'ausfüllen';
ua = 'ausfüllen';

Kaa = 'ausfüllen';
Kab = 'ausfüllen';
Kba = 'ausfüllen';
Kbb = 'ausfüllen';

ub = 'ausfüllen';
fa = 'ausfüllen';

%% Aufgabe 8
%'ausfüllen'