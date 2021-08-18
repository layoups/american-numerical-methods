clear all; 

sigma = 0.8;
r     = 0.25; 
E     = 10.0; 
D0    = 0.2; 
k     = (r-D0)/(0.5*sigma^2); 
T     = 1;                 

SRight = 30.0;                 
SLeft  = 1e-9;                  
xLeft  = log(SLeft/E); 
xRight = log(SRight/E); 

Nx = 1000;
dx = (xRight-xLeft)/Nx; 

a  = 1.0; 

dtau = a*dx^2;
tau_Max = (0.5*sigma^2)*T; 
M       = ceil(tau_Max/dtau);

[u,xgrid] = crank_fd_PSOR(@tran_payoff_call, @u_m_inf_call, @u_p_inf_call, r-D0, sigma, xLeft, xRight, Nx, tau_Max, M );

S   = E*exp( xgrid ); 
t   = 0;                        
tau = 0.5*(sigma^2)*(T-t);     

Spow = (S.^(0.5*(1-k))); 
Smat = repmat( Spow(:).', [M+1, 1] ); 
V  = (E^(0.5*(1+k))) * Smat * exp( -(1/4)*((k+1)^2)*tau ).*u; 

fh=gcf; 
figure(fh); as=plot( S, V(end,:), '-or' ); grid on; hold on; xlabel( 'S' ); ylabel('C'); 

[C,P] = blsprice(S, E, r, T, sigma, D0);  

figure(fh); bss=plot( S, C, '-k', 'LineWidth', 2 ); 
figure(fh); es = plot( S, max(S - E, 0), '-b', 'LineWidth', 2 ); 
 
legend( [ es,as,bss ], {'Payoff', 'American Call', 'Black-Scholes Analytic Solution'}, 'location', 'northwest' ); 
axis( [0,30,0,22] ); 

Sgrid = [ 1e-9:5:30 ]; 
Vgrid_eu = interp1( S, C(end,:), Sgrid, 'linear', 'extrap' );
Vgrid_am = interp1( S, V(end,:), Sgrid, 'linear', 'extrap' ); 
Vgrid_payoff = max(Sgrid - E, 0);

fprintf('row 1 = S\n');
fprintf('row 2 = American Call\n');
fprintf('row 3 = European Call\n'); 
[Sgrid; Vgrid_am; Vgrid_eu]