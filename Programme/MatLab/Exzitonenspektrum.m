Input   = struct( 'I_k',   [0 1 10 100 1000 10000], 'N_k',   [100 100 100 100 100], ...
                  'I_phi', [-pi pi],     'N_phi',  50 ); 

%[X , hw]  = Exzitonenanregung(Input,'Spektrum',  'Coulomb',-700:50) ; 
%[EW, ~ ]  = Exzitonenanregung(Input,'Eigenwerte','Coulomb',-700:50) ; 

plot(hw, X/max(X), 'r', ...
    'linewidth',        2.5); hold on
stem(EW(EW<0), ones(1, length(EW(EW<0))), 'r--', ...
    'markeredgecolor', 'none', ...
    'linewidth',        1.3)

xlim([min(hw) max(hw)]);    xlabel('$E-E_G$ in meV', 'interpreter','latex');
ylim([0 1]);                ylabel('Im$(\chi)$'    , 'interpreter','latex')
set(gca, 'fontsize', 13)
%set(gcf, 'units',   'centimeters', 'position', [0 0 13 5])
 
