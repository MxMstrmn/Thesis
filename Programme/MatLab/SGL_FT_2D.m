% Pysikalische Konstanten
c = constants; 
factor = -0.5*c.e^2 /c.eps0 ;  
mu = c.me; 

% Dimension der k Matrix und Stützstellen für das phi-Intergral
I = [0 2000 10000]; 
N_k = [200 10c0]; 
N_phi = [50];

[k,g_k] = integrate(I,N_k,4);
[phi, g_phi] = integrate([-pi pi],N_phi,6);

% Erzeugen des 3D Gitters (k',k,phi)
[K1,K,PHI] = meshgrid(k,k,phi);

% Potentialfunktion und kintetische Energie
V_ij        = @(k,k1,phi) ((factor*k1./sqrt(k1.^2+k.^2-2*k.*k1.*cos(phi))).^(1-eq(k,k1))-eq(k,k1)) ;
T_ii        = @(k,k1)       0.5*c.hbar^2/mu*k.^2       .*eq(k,k1) ;

% Terme die durch Addition der "1" hinzukamen
conv_factor = @(k,k1)      (4*k.^4)./(k1.^2+k.^2).^2 ; 
v1_ij        = @(k,k1,phi)  -V_ij(k,k1,phi).*conv_factor(k,k1) ;
v1_ii        = @(k,k1)       factor*12.00152*k .*eq(k,k1); 

% Anpassung der Gewichte auf die 3Dimensionen
dim_k           = sum(N_k);
dim_phi         = sum(N_phi);
weight_k1       = repmat(g_k',[dim_k,1,1]);
weight_phi      = permute(repmat(g_phi',[dim_k,1,dim_k]),[1,3,2]);

% Zusammenfügen der einzelnen Bestandteile für: 
% das Potential 
h_nurPot =      sum(V_ij(K,K1,PHI).*weight_phi,3).*weight_k1 ;
h_mitKF  = diag(sum(v1_ij(K,K1,PHI).*weight_phi,3)*g_k) ; 
h_Int    =          v1_ii(K(:,:,1),K1(:,:,1)); 

H_V      = h_nurPot + h_mitKF + h_Int;

% die kinetische Energie 
H_T      = T_ii(K(:,:,1),K1(:,:,1));

% Der Hamiltonian 
H = H_V + H_T;

%%
[states, eig_val]   = eig(H,'vector'); 
[eig_val, idx]      = sort(eig_val);
states              = states(:,idx);

disp(eig_val(1))

