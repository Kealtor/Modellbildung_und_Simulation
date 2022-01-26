%% Bearbeitungsbogen
clear all
close all
clc

%% Parameter
l = 0.2;                                                           % Laenge l
g = 9.81;                                                          % Erdbeschleunigung
omega = 2*sqrt(g/l);                                             % Rotationsgeschwindigkeit
AnzahlSchritte = 10^4;                                             % Anzahl der Zeitschritte
r = 2*pi-0;                                                        % Range von 0 - 2*pi
h = r/AnzahlSchritte;                                              % Zeitschrittweite

x_0 = [pi/6;0];                                                    % Anfangsbedingungen in Zustandsvektor
t = 0:h:r;                                                         % Zeitvektor

x_Eu_expl(:,1) = x_0;                                              % Anfangsbedingungen fuer Euler explizit
x_Eu_impl(:,1) = x_0;                                              % Anfangsbedingungen fuer Euler implizit
x_RuKu(:,1) = x_0;                                                 % Anfangsbedingungen fuer Runge-Kutta-Verfahren

%% Aufgabe 1: Systemmatrix
SystMatr = [0, 1; omega^2-(g/l), 0];                               % Systemmatrix
I = eye(size(SystMatr,1));                                         % Einheitsmatrix

gamma = det(SystMatr);
alpha = 0.5*trace(SystMatr);

lambda1 = sqrt(omega^2-(g/l)); 
lambda2 = -sqrt(omega^2-(g/l));

disp(['Lambda1: ',num2str(lambda1)]);
disp(['Lambda2: ',num2str(lambda2)]);

%disp(gamma);
%disp(alpha);
disp(['h*Lambda1: ',num2str(h*lambda1)]);
disp(['h*Lambda2: ',num2str(h*lambda2)]);


%% Aufgabe 2: exakte Loesung
x_exakt = [(pi/12) * (exp(sqrt(omega^2-g/l)*t) + exp(-sqrt(omega^2-g/l)*t)) ; 
           (sqrt(omega^2-g/l)* (pi/12) * exp(sqrt(omega^2-g/l)*t)) + (-sqrt(omega^2-g/l)* (pi/12) * exp(-sqrt(omega^2-g/l)*t))];

for n = 1:AnzahlSchritte
    
    %% Aufgabe 3: Euler explizit
    x_Eu_expl(:,n+1) = I*x_Eu_expl(:,n) + (h.*SystMatr)*x_Eu_expl(:,n);
    
    %% Aufgabe 4: Euler implizit
    x_Eu_impl(:,n+1) = inv(I-h*SystMatr) * x_Eu_impl(:,n);
    
    %% Aufgabe 5: Runge-Kutta-Verfahren
    k1(:,n) = h.*(SystMatr*x_RuKu(:,n));
    k2(:,n) = h.*(SystMatr* (x_RuKu(:,n)+0.5.*k1(:,n)));
    k3(:,n) = h.*(SystMatr* (x_RuKu(:,n)+0.5.*k2(:,n)));
    k4(:,n) = h.*(SystMatr* (x_RuKu(:,n)+k3(:,n)));
    x_RuKu(:,n+1) =  x_RuKu(:,n)+(1/6).*(k1(:,n)+2*k2(:,n)+2*k3(:,n)+k4(:,n));
    
end

%Local
l_Eu_expl = x_exakt - x_Eu_expl;
l_Eu_impl = x_exakt - x_Eu_impl;
l_RuKu = x_exakt - x_RuKu;

%Global

for i=1:AnzahlSchritte
    e_Eu_expl(1,i) = max(abs(x_Eu_expl(1,1:i)-x_exakt(1,1:i)));
    e_Eu_expl(2,i) = max(abs(x_Eu_expl(2,1:i)-x_exakt(2,1:i)));
    
    e_Eu_impl(1,i) = max(abs(x_Eu_impl(1,1:i)-x_exakt(1,1:i)));
    e_Eu_impl(2,i) = max(abs(x_Eu_impl(2,1:i)-x_exakt(2,1:i)));
    
    e_RuKu(1,i) = max(abs(x_RuKu(1,1:i)-x_exakt(1,1:i)));
    e_RuKu(2,i) = max(abs(x_RuKu(2,1:i)-x_exakt(2,1:i)));

end

%% Aufgabe 7: Plots
figure(1);  
subplot(2,1,1);
plot(t,x_exakt(1, :),'Color','g','LineWidth',3, 'LineStyle', ':');     
hold on;                         
plot(t,x_Eu_expl(1,:),'Color','r','LineWidth',1);
plot(t,x_Eu_impl(1,:),'Color','b','LineWidth',1); 
plot(t,x_RuKu(1,:),'Color','k','LineWidth',1);     
legend('exakt','Expl', 'Impl', 'RuKu');
title('Exact phi and Simulation(Integration methods) Plot')
xlabel('Zeit');
hold off

subplot(2,1,2); 
plot(t,x_exakt(2, :),'Color','g','LineWidth',3, 'LineStyle', ':');     
hold on;                         
plot(t,x_Eu_expl(2,:),'Color','r','LineWidth',1);
plot(t,x_Eu_impl(2,:),'Color','b','LineWidth',1); 
plot(t,x_RuKu(2,:),'Color','k','LineWidth',1);     
legend('exakt','Expl', 'Impl', 'RuKu');
title('Exact phi dot and Simulation(Integration methods) Plot');
xlabel('Zeit');
hold off
%% Aufgabe 8: Plots
figure(2); 
subplot(2,1,1);
plot(t,l_Eu_expl(2,:),'b-','LineWidth',1);
hold on
plot(t,l_Eu_impl(2,:),'g-','LineWidth',1);
plot(t,l_RuKu(2,:),'r-','LineWidth',1);
legend('Expl', 'Impl', 'Ruku');
title('Local Error');
hold off

subplot(2,1,2);
plot(t(1:end-1),e_Eu_expl(2,:),'b-','LineWidth',4);
hold on
plot(t(1:end-1),e_Eu_impl(2,:),'g-','LineWidth',2);
plot(t(1:end-1),e_RuKu(2,:),'r-','LineWidth',1);
legend('Expl', 'Impl', 'Ruku');
title('Global Error');
hold off