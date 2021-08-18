sigma = 0.4;
r     = 0.1; 
E     = 10.0; 
k     = r/(0.5*sigma^2); 
T     = 4/12;                

SRight = 30.0;                
SLeft  = 1e-9;                 
xLeft  = log(SLeft/E); 
xRight = log(SRight/E); 

Nx = 1000;
dx = (xRight-xLeft)/Nx; 

a  = 0.5; 

dt = a*dx^2;
tau_Max = (0.5*sigma^2)*T; 
M       = ceil(tau_Max/dt);

S_values = [ SLeft:SRight ];
[C, P] = blsprice(S_values, E, r, T, sigma); 

Sgrid = [ 1.e-9:2:18 ]; 
Vgrid_eu = interp1( S_values, P(end,:), Sgrid, 'linear', 'extrap' );

[u,xgrid] = crank_fd_PSOR(@tran_payoff_put, @u_m_inf_put, @u_p_inf_put, r, sigma, xLeft, xRight, Nx, tau_Max, M );

S   = E*exp( xgrid ); 
t   = 0;                       
tau = 0.5*(sigma^2)*(T-t);     

Spow = (S.^(0.5*(1-k))); 
Smat = repmat( Spow(:).', [M+1, 1] ); 
V    = (E^(0.5*(1+k))) * Smat * exp( -(1/4)*((k+1)^2)*tau ).*u;  

Vgrid_am = interp1( S, V(end,:), Sgrid, 'linear', 'extrap' ); 

fprintf('row 1 = S\n');
fprintf('row 2 = European Put\n');
fprintf('row 3 = American Put\n'); 
[ Sgrid; Vgrid_eu; Vgrid_am ]