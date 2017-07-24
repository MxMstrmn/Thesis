%==============================================================
%          Untersuchung von 1/sqrt(k²+k'²-2kk'cos(y))
%   Frage: Wo können bei der num. Beh. Probleme auftreten?
%==============================================================

% Das ist der Ausdruck der in der Wurzel im Nenner steht, fcn muss also für
% alle k,k',phi >0 sein. Diese Ausgangsfunktion wird im Schritt zu fcn2 für
% die einfachere Handhabung umgeschrieben. 
fcn     = @(k,k1,phi) k.^2 + k1.^2 -2*k.*k1.*cos(phi);

% Für das Verhältnis durch k1^2 teilen, k1 = 0 ist ausgeschlossen, das ist
% analytisch zu betrachten (k=0 oder k'=0, nicht beide, ist kein Problem, 
% da fcn >0)
fcn2    = @(u,phi) u.^2 + 1 -2*u.*cos(phi); 


phi = linspace(-pi,pi,100);
u   = linspace(0,1,100); 

[PHI,U] = meshgrid(phi,u); 
z = fcn2(U,PHI); 

surface(PHI,U,z,'LineStyle','none','FaceLighting','phong')%'EdgeColor','none',
xlabel('\phi')
ylabel('k/k1')
zlabel('Funktionswert')
xlim([-pi,pi])
ylim([0,1])

%==============================================================
%                           Ergebnis
%==============================================================
% Jeweils die 2D Plots ansehen, z(phi) & z(u), dann sofort ersichtlich,
% dass Problem bei k=k' v phi=0, genauso ist ein Problem bei k=k'=0
