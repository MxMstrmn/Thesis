% Pysikalische Konstanten
c = constants; 
factor = -0.5*c.e^2 /c.eps0 /pi.^2 * 0.25 ;  
mu = c.me; 

% Dimension der k Matrix und Stützstellen für das phi-Intergral
I = [0 1 10 100 1000 10000]; 
N_k = [2 2 2 2 2];
% Der Energiewert konvergiert schon mit [20 20 20 10 10], aber für die
% Plots sind etwas mehr Stützstellen schäner anzusehen
N_phi = [50];

[k,g_k] = integrate(I,N_k,4);
[phi, g_phi] = integrate([-pi pi],N_phi,6);

% Erzeugen des 3D Gitters (k',k,phi)
% Es muss ein solches Gitter erstellt werden, da die Faltung des Potentials
% mit der Wellenfunktion 
[K1,K,PHI] = meshgrid(k,k,phi);

% Potentialfunktion und kintetische Energie
V_ij        = @(k,k1,phi) ((factor*k1./sqrt(k1.^2+k.^2-2*k.*k1.*cos(phi))).^(1-eq(k,k1))-eq(k,k1)) ; % 0 entlang der Diagonalen
T_ii        = @(k,k1)       0.5*c.hbar^2/mu*k.^2       .*eq(k,k1) ;

% Terme die durch Addition der "1" hinzukamen
conv_factor  = @(k,k1)      (4*k.^4)./(k1.^2+k.^2).^2 ; 
v1_ij        = @(k,k1,phi)  -V_ij(k,k1,phi).*conv_factor(k,k1) ;
v1_ii        = @(k,k1)       factor*12.00152*k .*eq(k,k1); 

% Anpassung der Gewichte auf die 3Dimensionen
dim_k           = sum(N_k);
dim_phi         = sum(N_phi);
weight_k1       = repmat(g_k',[dim_k,1,1]);
weight_phi      = permute(repmat(g_phi',[dim_k,1,dim_k]),[1,3,2]);

% Zusammenfügen der einzelnen Bestandteile für: 
% das Potential 
h_nurPot =      sum(V_ij(K,K1,PHI).*weight_phi,3).*weight_k1 ;  % Dranmultpl. der Gewichte, phi Int ausgeführt, Hier ist Faltung also Psi noch im Integral
h_mitKF  = diag(sum(v1_ij(K,K1,PHI).*weight_phi,3)*g_k) ;       % Hier keine Faltung deswegen k Int ausgeführt
h_Int    =          v1_ii(K(:,:,1),K1(:,:,1));                  % Wegen Phi sind K und K1 3Dim, benötigt man aber nicht für die Integralmatrix, daher (:,:,1)

H_V      = h_nurPot + h_mitKF + h_Int;

% die kinetische Energie 
H_T      = T_ii(K(:,:,1),K1(:,:,1));

% Der Hamiltonian 
H = H_V + H_T;

%%
% Bestimmung der Eigenwerte (eig_val) samt Normierung der Wellenfunktionen
% (states). Beides beginnend mit dem Grundzustand (sort).  
[states, eig_val]   = eig(H,'vector'); 
[eig_val, idx]      = sort(eig_val);
states              = states(:,idx);
norm = sqrt(2*pi*(states.^2)'*(k.*g_k));
for i=1:dim_k; states(:,i) = states(:,i)*1/norm(i)*sign(states(1,i)); end

% Anzeigen der Grundzustandsenergie 
disp(eig_val(1))
%%
% Fouriertransformation der analytischen 2D Wellenfunktionen 
fcn1 = @(r) 4/c.a0.*exp(-2/c.a0*r)*1/(4*pi).^(1/3) ; 
fcn2 = @(r) 4/c.a0/3 /sqrt(3) .*(1-4*r/3 /c.a0) .*exp(-2/3 /c.a0*r)*1/(4*pi).^(1/3) ; 

[x1,y1] = FT2([0 2 10],[200 200],[0 150],[60],fcn1);
[x2,y2] = FT2([0 2 10],[200 200],[0 150],[60],fcn2);

%%
% Darstellen der Ergebnisse
en_num = eig_val(eig_val<0);
en_ana = Rydberg(1:15,1,1,2)';


figure 
rows = 3;
subplot(rows,1,1); num1 = plot(k,states(:,1),'b.',x1,y1,'r-',k,0*k,'k-'); xlim([0 150]); ylim([0 0.025]); xlabel('$k$ in nm$^{-1}$','interpreter','latex'); ylabel('$\Psi_1 (r)$','interpreter','latex')
subplot(rows,1,2); num2 = plot(k,states(:,2),'b.',x2,y2,'r-',k,0*k,'k-'); xlim([0 40]);  ylim([-0.05 0.125]); xlabel('$k$ in nm$^{-1}$','interpreter','latex'); ylabel('$\Psi_2 (r)$','interpreter','latex')
subplot(rows,1,3); num3 = plot(k,states(:,3),'b.',k,0*k,'k-');            xlim([0 20]);  ylim([-0.1 0.25]); xlabel('$k$ in nm$^{-1}$','interpreter','latex'); ylabel('$\Psi_3 (r)$','interpreter','latex')
legend([num1(1) num1(2)],sprintf('Eigenfunktion zum EW=%4.1f eV', en_num(1)/1000),'Analytisches Ergebnis');
legend([num2(1) num2(2)],sprintf('Eigenfunktion zum EW=%4.1f eV', en_num(2)/1000),'Analytisches Ergebnis');
legend([num3(1)        ],sprintf('Eigenfunktion zum EW=%4.1f eV', en_num(3)/1000))

%subplot(4,1,4); plot(k,states(:,1:3),'r-',k,0*k,'k-');          xlim([0 8]);

figure 
num4 = plot([0 1],[en_num en_num]*1e-3,'b-',[1 2],[en_ana en_ana]*1e-3,'r-');
set(gca,'xticklabel',[],'xtick',[]); 
ylabel('Energieeigenwerte in eV','interpreter','latex')
legend([num4(1) num4(end)],'numerisch','analytisch','location','best')

