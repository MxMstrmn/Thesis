%=====================================================%
%           FT - Schrödinger-equations                %
%=====================================================%
function [EW, states, k] = SGL_FT(I,N)
% constants
hbar = 0.6582119514;
e    = 1.602176565e5;
me   = 2*5.686e-3;  
mh   = me;
eps0 = 1.41844e6;
coulomb = - e^2/(4*pi^2*eps0) ; 

% possible values for interpolation points
% I = [0,10,100,200,1000];
% N = [200,400,300,100];


% potential & energy functions
V = @(k,k1) (coulomb * k1./k .* log(abs((k+k1)./(k-k1)))).^(1-eq(k,k1))-eq(k,k1); 
E = @(k,k1,m) (hbar*k).^2./(2*m).*eq(k,k1) ;

% comparison of different integrators
results = struct();

for i=5:5
% k values and weigthings g
[k, g]  = integrate(I,N,i); 
[K, K1] = meshgrid(k);
weight  = repmat(g',length(k),1);

%=====================================================================%
% N O T I Z 
%
% Es ist sehr wichtig wie man die Gewichte an das Potential
% dranmultipliziert, denn trotz identischer Eigenwerte verändern sich die
% Wellenfunktionen! Das lässt sich dann aus dem Verlauf aber auch erahnen!
%
%=====================================================================%

% Hamiltonian
dim = size(K);
H = reshape(  V(K(:),K1(:)).*weight(:) + E(K(:),K1(:),me) + E(K(:),K1(:),mh),  dim);


% Trapez und Kepler gehen nicht, da Randwerte in den x_i doppelt vorkommen

[states, energy] = eig(H,'vector');
%============================================================%
% Just in case one wants to compare the different integrators
%
% switch i
%     case {1} %Rechteck
%         results.H_Rechteck = energy;
%     case {2} %Trapez
%         results.H_Trapez = energy;
%     case {3} %Kepler
%         results.H_Kepler = energy;
%     case {4} %GTS
%         results.H_GTS = energy;
%     case {5} %Gauss
%         results.H_Gauss = energy;
%
%============================================================%

end

% try
%     result = sort([results.H_Rechteck results.H_Trapez results.H_Kepler results.H_GTS results.H_Gauss]);
% catch
%     %disp('Nicht alle Integratoren benutzt')
% end

[EW, idx]   = sort(energy);
states      = states(:,idx);
disp(EW(1))

end


