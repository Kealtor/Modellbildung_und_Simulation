clear all;
close all;
clc

%% Parameter
L = 1; %[m]
F = 20000; %{N}
E = 70*(10^9);  %[Pa]
LElemente = (L/AnzahlElemente);       % Länge eines Elemntes
AnzahlElemente = 5;
Knoten = (0:L/AnzahlElemente:L);      % Knotenvektor
AElement = 0.001-0.0009*Knoten;     % Querschnittsverlauf

%% Aufgabe 2
syms u(x);
ode = diff(u) == F/(E*(0.001-0.0009*x));  %du = strain dx; du = F/(EA)dx
cond = u(0) == 0;
uSol(x) = dsolve(ode,cond);
%% Aufgabe 3
N = [1-(Knoten'/LElemente) Knoten'/LElemente];

%% Aufgabe 4
%B = diff(N);
B = diff(N)./diff(([Knoten;Knoten])');
C = E;

%% Aufgabe 5

for b=1:AnzahlElemente                    % Schleife für die Anzahl an gewählten Elementen
   A = (AElement(b)+AElement(b+1))/2;  
   K{b} = A*LElemente*(B(b,:))'*C*B(b,:);
   
end

%% Aufgabe 6
KGesamt = zeros(length(Knoten),length(Knoten));
for b=1:AnzahlElemente  
    KGesamt(b:b+1,b:b+1) = KGesamt(b:b+1,b:b+1)+K{b};
end

%% Aufgabe 7
fb = [zeros(AnzahlElemente-1,1);F];
ua = 0;   

Kaa = KGesamt(1,1);
Kab = KGesamt(1,2:end); 
Kba = KGesamt(2:end,1);
Kbb = KGesamt(2:end,2:end);

ub = (Kbb^(-1))*(fb-Kba*ua);
fa = Kaa * ua + Kab * ub;  

u_disp_vec = [ua;ub];
force_vec = [fa;fb];

%% Aufgabe 8
%uSol_L = uSol(1)

%displacement plot
subplot(2,1,1)
x=linspace(0,L,10);
plot(x,uSol(x),'-r','linewidth',1)
hold on
plot(Knoten,u_disp_vec,'-+b','linewidth',1)
title('Displacement Comparison')
legend('Exact Soln','Numerical Soln')
xlabel('Position [m]')
ylabel('Displacement [m]')

%stress diagram 
stress = E*diff(uSol);
stress_num=E*diff(u_disp_vec)./diff((Knoten)');

subplot(2,1,2)
x=linspace(0,L,10);
plot(x,stress(x),'-r','linewidth',1)
hold on
stairs(Knoten,[stress_num;stress_num(end)],'-+b','linewidth',1)  %because the size isnt the same also since it is a stairs the last value is the same
title('Stress Comparison')
legend('Exact Soln','Numerical Soln')
xlabel('Position [m]')
ylabel ('Stress [Pa]')

