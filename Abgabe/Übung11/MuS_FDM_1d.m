function MuS_FDM_1d
%
% Modellbildung und Simulation
% Übung zum Kapitel: Zeitkontinuierliche Modelle mit verteilten Parametern
%
% Dr.-Ing. Balazs Pritz, pritz@kit.edu
%
% 1D Testfall: Konvektions-Diffusionsgleichung mit Dirichlet-Randbedingungen
%
% Diskretisierungsschema: Finite Differenzen Methode
% 
% Quelle: Joel H. Ferziger, Milovan Peric:
% Numerische Strömungsmechanik, Springer, 2008

clear all

% Zeitschritt
dt=0.0013;
% Konvektionsgeschwindigkeit
u0=2.0;
% Anzahl der Knotenpunte
nx=11;
% Materialeigenschaften
% Dichte
% Luft
rho=1.0;
% Wasser
%rho=998;

% Diffusionskoeffizient für phi
gamma=1;

% Maximale Anzahl von Zeitschritten
% (Simulationszeit=timestep*maxnt)
% (Es ist auch möglich, dass die Simulation nicht nur durch maxnt beendet
% wird, sondern die Änderung der Lösung von zwei Zeitschritten berechnet
% wird und anhand dieser ein Abbruchkriterium definiert wird.)
maxnt=200;

% Alternativ können Sie die Simulationszeit definieren und maxnt ausrechnen
% lassen.



% Räumliche Auflösung
% Problemgebiet: x=[0,1]
%
% Es wird hier ein äquidistantes Netz verwendet.
% Wenn Sie lokale Verfeinerung benutzen möchten, dann müssen Sie
% für x(i) einen Vektor definieren, und in den Gleichungen dx durch
% x(i)-x(i-1), x(i+1)-x(i) oder x(i+1)-x(i-1) ersetzen.
%


nxm1=nx-1;

% Zellengröße
dx=1/(nxm1);



% Initialisierung
phi_uds=zeros(nx,1);
phi_cds=zeros(nx,1);

% Randbedingungen für phi
phi_uds(1)=0;
phi_uds(nx)=1;

phi_cds(1)=0;
phi_cds(nx)=1;


% Peclet-Zahl
% Für die Änderung der Pe-Zahl können Sie natürlich beliebig rho, u oder gamma ändern
peclet=rho*u0/gamma;


% Die exakte Lösung (wird mit besserer Auflösung gerechnet)
xel=zeros(100,1);
el=zeros(100,1);
for i=1:100
    xel(i)=0.01*i-0.01;
    el(i)=(exp(xel(i)*peclet) - 1.)/(exp(peclet) - 1.);
end

% Zur Beurteilung der Stabilität wird DCFL gerechnet.
% (Wenn Sie eine lokale Verfeinerung verwenden, müssen Sie nach dem größten Wert suchen.)
dcfl=2*gamma*dt/(dx*dx*rho) + u0*dt/dx;


% Die Position der Stützstellen wird berechnet
x=zeros(nx,1);
for i=1:nx
    x(i)=dx*i-dx;
end


% Initialisierung von Vektoren
rhs_uds=zeros(nx,1);
rhs_cds=zeros(nx,1);


% Zeitschleife
for nt=1:maxnt
    
    for i=2:nxm1
        % Konvektiver Term mit Upwind Differenzen Schema
        konv_uds=rho*u0*((phi_uds(i)-phi_uds(i-1))/dx);
        
        % Konvektiver Term mit zentralem Differenzen Schema
        konv_cds=rho*u0*((phi_cds(i+1)-phi_cds(i-1))/(2*dx));
        
        % Diffusiver Term mit CDS
        diff_uds=gamma* (phi_uds(i+1)-2*phi_uds(i)+phi_uds(i-1))/dx^2 ;
        diff_cds=gamma* (phi_cds(i+1)-2*phi_cds(i)+phi_cds(i-1))/dx^2 ;

        rhs_uds(i)=diff_uds-konv_uds;
        rhs_cds(i)=diff_cds-konv_cds;

    end
    
    % Kalkuliere neue Werte für phi in der Zeit mit einem
    % expliziten Euler Schema
    for i=2:nxm1
        phi_uds(i)= phi_uds(i)+rhs_uds(i)*dt; 
        phi_cds(i)= phi_cds(i)+rhs_cds(i)*dt;
    end

    % Simulationszeit
    zeit=nt*dt;

    % Ergebnisse werden dargestellt
    % Nach jeder n-ten Iteration soll das Ergebnis geplottet werden
    % (werden die Ergebnisse nicht in jeder Iteration dargestellt, läuft
    % die Simulation schneller)
    n=5;
    if rem(nt,n)==0
        figure(1)
        plot(xel,el,'-b',x,phi_uds,'-ro',x,phi_cds,'-kx');
        legend('exakt','uds','cds','Location','Eastoutside');
        axis([0 1 -0.2 1])
        text(0.1,0.9,['Peclet-Zahl= ',num2str(peclet)])
        text(0.1,0.85,['DCFL-Zahl= ',num2str(dcfl)])
        text(0.1,0.8,['Zeit= ',num2str(zeit)])
        %drawnow %in neuerer Version von Matlab soll diese Funktion zusätzlich benutzt werden
    end
   
end
%Ende Zeitschleife

end
% Ende Funktion

