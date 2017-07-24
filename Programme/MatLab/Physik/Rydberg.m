function energy = Rydberg(n,m,eps,dim)
c      = constants; 
mu     = m*c.me;

switch dim 
    case {3}
        energy = -mu*c.e^4/(32*pi^2*c.eps0^2*eps^2*c.hbar^2).*1./n.^2; 
    case {2}
        energy = -mu*c.e^4/(32*pi^2*c.eps0^2*eps^2*c.hbar^2).*1./(n-0.5).^2;
end