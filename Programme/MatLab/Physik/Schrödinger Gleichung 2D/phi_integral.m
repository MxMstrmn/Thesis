k = 1;  
con = -9.0486e+03; 

I = [0 15 100];
N = [200 20]; 

[u, g_u]       = integrate(I,N,4); 
[phi, g_phi]    = integrate([-pi pi],1000,6); 

[PHI,U] = meshgrid(phi,u);

%conv_factor = @(k,k1)      (4*k.^4)./(k1.^2+k.^2).^2 ; 
%fcn         = @(k,k1,phi) con *k1 .*(k1.^2 + k.^2 -2*k.*k1.*cos(phi)).^(-0.5).*conv_factor(k,k1) ; 

fcn = @(u,phi) 4*u .*(u.^2+1-2*u.*cos(phi)).^(-0.5) .* (1+u.^2).^(-2);

(fcn(U,PHI)*g_phi)'*g_u; 
