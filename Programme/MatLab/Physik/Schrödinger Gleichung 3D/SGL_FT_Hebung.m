%=====================================================%
%      FT - SGL mit Hebung der Singularität           %
%=====================================================%

function [EW, states, k] = SGL_FT_Hebung(I,N)
% Einführung der Konstanten mit Hilfe der constants Funktion/ Struktur
c       = constants; 
m1      = 2*c.me;
m2      = m1;
coulomb = -c.e^2/(4*pi^2*c.eps0) ; 
mu      = m1*m2/(m2+m1);


% Die verschiedenen Potentialfunktionen und die kintetische Energie.
% v beschreibst das "normale" Potential, welches man auch ohne die Hebung
% der Sigularität betrachtet. v1 ist das selbe Potential multipliziert mit dem
% Konvergenzerzeugenden Faktor 4k⁴/(k²+k'²)². vcons beschreibst das analytisch 
% ausgerechnete Intergral, welches den Anteil der Singularität beschreibt (coumlomb*pi*k).
% E ist der Anteil der kintetischen Energie.

v       = @(k,k1)   (coulomb*k1./k .* log(abs((k+k1)./(k-k1)))).^(1-eq(k,k1))-eq(k,k1) ;  
v1      = @(k,k1)    v(k,k1).* (4*k.^4)./(k1.^2+k.^2).^2 ;      
E       = @(k,k1,m) (c.hbar*k).^2/(2*m).*eq(k,k1) ; 
vcons   = @(k,k1)    coulomb*pi*k.*eq(k,k1) ;

% B E M E R K U N G 
%=============================
% Bei repmat sehr aufmerksam sein, denn hier spielt es eine Rolle wie die
% Gewichte multipliziert werden.
%=============================
[k,g] = integrate(I,N,4); 
[K,K1] = meshgrid(k);
weight = repmat(g',size(g));
dim = size(K);


% B E M E R K U N G 
%=============================
% Auch hier ist es sehr wichtig zu wissen, welche Werte welchen
% Matrixelementen zugeordnet werden. Die ersten Werte von K gehören zu k(1)
% und damit die Reihen der Matritzen k sind und die Spalten k' muss man die 
% Matritzen transponieren. Nur dann ist auch das Skalarprodukt mit den
% Gewichtungen g in H_V_ii korrekt. 
%=============================
H_kin       = reshape( E(K(:),K1(:),mu)   ,dim); 
H_V_ii      = reshape( vcons(K(:),K1(:))  ,dim) - diag(reshape(v1(K(:),K1(:)),dim)'*g);
H_V_ij      = reshape( v(K(:),K1(:))      ,dim)'.*weight;

% Zusammenstellen des Hamiltonians 
H = H_kin + H_V_ii +  H_V_ij;


% Bestimmung der Eigenfunktionen und Eigenwerte
[states, values] = eig(H,'vector'); 

[EW, idx]   = sort(values);
states      = states(:,idx);
disp(EW(1))
norm = sqrt(4*pi*(states.^2)'*(k.^2.*g));
for i=1:dim(1); states(:,i) = states(:,i)*1/norm(i)*sign(states(1,i)); end

end

