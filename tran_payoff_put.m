function g = tran_payoff_put(xgrid,tau,k)
xgrid = xgrid(:); 
nx    = length(xgrid); 
g = exp( 0.25*((k+1)^2)*tau ) * max( [ exp( 0.5*(k-1)*xgrid ) - exp( 0.5*(k+1)*xgrid ), zeros(nx,1) ], [], 2 ); 

