
figure

EXMX      = 'MX' ; 

switch EXMX
    case 'EX'
        %==================================================================================
        Input_ex    = struct(   'I_k',   [0 1 10 100 1000 10000], 'N_k',   [10 10 10 10 10], ...
                                'I_phi', [-pi pi],                'N_phi',  50 ) ; 
        %==================================================================================
        [X , hw]    = Exzitonenanregung(Input_ex,'Spektrum',  'Coulomb',-700:50) ; 
        [EW, ~ ]    = Exzitonenanregung(Input_ex,'Eigenwerte','Coulomb',-700:50) ; 

        plot_i      = plot(hw, imag(X)/max(imag(X)), 'r',                'linewidth', 2.5) ; hold on
        stem_i      = stem(EW(EW<0), ones(1, length(EW(EW<0))), 'r--',   'linewidth', 1.3, ... 
                            'markeredgecolor', 'none') ;
        
        xlim([min(hw) max(hw)]);    xlabel('$E-E_G$ in meV', 'interpreter','latex');
        ylim([0 1]);                ylabel('Im$(\chi)$'    , 'interpreter','latex');
        grid on;                    set(gca,'GridLineStyle', '--', 'fontsize',13);
    case 'MX'
        %==================================
        lambda      = [0.25 0.5 1 2 4 8]  ; 
        phi         = linspace(-8,20,500) ;
        %==================================
        
        for i=1:6   
        ax          = subplot(3,2,i);  
        Input_mx    = struct( 'n',300, 'lambda',lambda(i), 'phi',phi ) ;                 
        
        [X , hw]    = Magnetoexziton(Input_mx,'Spektrum',  'Coulomb') ;
        [EW, ~ ]    = Magnetoexziton(Input_mx,'Eigenwerte','Coulomb') ;
        
        plot_i      = plot(ax,hw, imag(X)/max(imag(X)), 'r',    'linewidth', 2.5) ; hold on
        stem_i      = stem(ax,EW, ones(1, length(EW)), 'r--',   'linewidth', 1.3, ...
                            'markeredgecolor', 'none') ;
        
        xlim([min(hw) max(hw)]); ylim([0 1]); 
        
        if i<5
            set(gca,'XTickLabel',[]); 
        else; xlabel(ax, '$(\hbar \rm{\omega}- E_G)/E_B$', 'interpreter','latex'); end
        
        if mod(i,2)==0
            set(gca,'YTickLabel',[]);  
        else; ylabel(ax, 'Im($\rm{\chi}$)', 'interpreter','latex'); end
        
        grid(ax,'on'); set(gca,'GridLineStyle', '--', 'fontsize',13);
        legend(plot_i,['\lambda = ' num2str(lambda(i))],'interpreter','latex');
        end
end

 
