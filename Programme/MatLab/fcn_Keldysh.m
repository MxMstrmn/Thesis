function out = fcn_Keldysh(x)

% Gen. Laguerre Polynomial
L   = @(n,a,x)    lf_function(length(x),n,a,x(:));

% Parameter
A   = 1; 
B   = A*2*pi*0.66; 

% Keldysh function for Gauss-Laguerre Integration
f   = @(n,a,x)    (L(n,a,x)).^2 ./(A+B*sqrt(2*x(:))) ;


n   = getGlobaln ;
n1  = getGlobaln1 ; 
out = f(n+1,n1-n,x);

% lf_function gibt die Werte für alle Ordnungen der Laguerre Pol. aus,
% beginnend mit n=0 
% % EDIT: lf_function gibt nur noch höchste Ordnung aus 
% % out = out(:,n+1);

end

