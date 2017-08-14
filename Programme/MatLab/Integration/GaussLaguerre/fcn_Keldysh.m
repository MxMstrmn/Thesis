function out = fcn_Keldysh(x,B)
c       = constants();

% Gen. Laguerre Polynomial
L           = @(n,a,x)  lf_function(length(x),n,a,x(:));
magnlen     = @(B)      sqrt(c.hbar/(B*c.e)) ; 

% Parameter 
alpha   = 1; 
beta    = 2.1; %alpha*2*pi*0.66; 
B       = B * c.unitB;

% Keldysh function for Gauss-Laguerre Integration
f   = @(n,a,x)    (L(n,a,x)).^2 ./(alpha+beta/magnlen(B)*sqrt(2*x(:))) ;


n   = getGlobaln ;
n1  = getGlobaln1 ; 
out = f(n,n1-n,x);

% lf_function gibt die Werte für alle Ordnungen der Laguerre Pol. aus,
% beginnend mit n=0 
% % EDIT: lf_function gibt nur noch höchste Ordnung aus 
% % out = out(:,n+1);

end

