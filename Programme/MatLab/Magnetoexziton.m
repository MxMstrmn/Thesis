const       = constants ; 
h           = const.hbar ; 
n           = 300 ;

lambda      = 2;
coufac      = (2*lambda/pi)^0.5 ;                              % Coulomb Vorfaktor



[n2,n1]     = meshgrid(0:n) ; 

VC_ij       = @(n)   -coufac *exp(gamma_prefactor(0:n)) .*F32(0:n) ;
Hmx_ii      = @(n1,n2)   lambda*(2*n1+1)  .*(n1==n2);
Inh_ii      = @(n1,n2,phi) (-phi -1i*0.2)  *(n1==n2);       % phi = hw-Eg / EB

b           = ones(n+1,1) ;
X           = [] ;
phi         = linspace(-10,45,1000);

H           = Hmx_ii(n1,n2) + VC_ij(n) ;
for i=1:length(phi)
    A       = H + Inh_ii(n1,n2,phi(i)) ;
    x_n     = linsolve(A,b) ;
    x       = sum(x_n)*lambda/pi;
    X       = [X x] ;
end

