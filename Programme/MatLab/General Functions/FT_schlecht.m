function [K,y] = FT_schlecht(I,N,method)
% I & N are values for the output function
Z=1;
a0=0.0529;
switch method 
    case {1}
        psi = @(r) sqrt(4*Z/a0^3)*exp(-Z.*r/a0)*sqrt(1/4 /pi);
    case {2}
        psi = @(r) sqrt(Z/8 /a0^3)*(2-Z.*r/a0).*exp(-Z.*r/2/a0)*sqrt(1/4 /pi);
    case {3}    
        psi = @(r) 1 /81 /sqrt(3*pi) *a0^(-3/2) *(27-18.*r/a0 + 2.*r.^2/a0^2).*exp(-r/3/a0);
end
% Data for Fourierintegral
I_F = [0 20 100 500]; 
N_F = [70 20 20]; 

y=[];
for i=1:N
    K = linspace(I(1)+eps,I(2),N) ; 
    k = K(i) ;
    % Bei dem Vorfaktor habe ich keine Ahnung
    vorfaktor = 4*pi/(2*pi)^(3/2); 
    weight = @(r) vorfaktor*r/k .* sin(k*r) ; 
    integrand = @(r) weight(r).*psi(r) ; 
    [ri, g]  = integrate(I_F,N_F,4,integrand);
    y(i) = integrand(ri)'*g ;
   
end
if y(1)<0
    y= -y; 
end