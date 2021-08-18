function [u,loops] = PSOR_solver( u, b, g, a, omega, n_eps )
loops = 0; 
error = 1e6; 
nelts = length(u); 
u(1)   = max( u(1), g(1) ); 
u(end) = max( u(end), g(end) ); 

while( error>n_eps )
  error = 0; 
 
  for n=2:(nelts-1)
    y     = ( b(n) + a * ( u(n-1) + u(n+1) ) )/(1+2*a); 
    y     = max( u(n) + omega * ( y - u(n) ), g(n) );
    error = error + (u(n)-y)^2;
    u(n)  = y;                                 
  end
    
  loops = loops+1; 

end
