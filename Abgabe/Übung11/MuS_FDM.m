import MuS_FDM_1.*

nx=200; % Anzahl der Knotenpunte
u0=4.0; % Konvektionsgeschwindigkeit
rho=1.0; %% Dichte Luft
gamma=0.1;% Diffusionskoeffizient für phi
%rho=998; %% Dichte Wasser 
start = 0.00001 ;
stop  = 0.0001;
step  = 0.00001;
dx = 1/(nx-1)
%dcfl=2*gamma*dt/(dx*dx*rho) + u0*dt/dx;
results = zeros(1,3);
dt_max= 1/ (2*gamma/(rho*dx*dx) + u0/dx)
%dt=0.001; % Zeitschritt
i=1;
for dt=start:step:stop
    
    [phi_uds,phi_cds,x,el,xel,peclet,dcfl,zeit] = MuS_FDM_1(dt,nx,rho,u0,gamma);
    nxm1 = nx-1;
    num_idx = 1 + 10*nxm1;
    
    el_filtered = el(1:10:num_idx);
    e_uds = (el_filtered - phi_uds) ./el_filtered;
    e_cds = (el_filtered - phi_cds) ./el_filtered;
    
    e_uds = fillmissing(e_uds,'constant',0);
    e_cds = fillmissing(e_cds,'constant',0);
    
    MAPE_uds = sum(abs(e_uds))/nx*100;
    MAPE_cds = sum(abs(e_cds))/nx*100;
%     disp(dt)
%     disp([MAPE_uds,MAPE_cds,dt])
    results(i,:) = [MAPE_uds,MAPE_cds,dt];
    i=i+1;
    figure(1)
    plot(xel,el,'-b',x,phi_uds,'-ro',x,phi_cds,'-kx');
    legend('exakt','uds','cds','Location','Eastoutside');
    axis([0 1 -0.2 1])
    text(0.1,0.9,['Peclet-Zahl= ',num2str(peclet)])
    text(0.1,0.85,['DCFL-Zahl= ',num2str(dcfl)])
    text(0.1,0.8,['Zeit= ',num2str(zeit)])
end

a =min(results,[],1);
[ridx_uds cidx_uds]= find(results==a(1));
[ridx_cds cidx_cds] = find(results==a(2));
opt_dt_uds = results(ridx_uds,3)
opt_dt_cds = results(ridx_cds,3)