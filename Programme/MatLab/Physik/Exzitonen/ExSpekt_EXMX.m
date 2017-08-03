% clear all
% clc 
figure
 
EXMX      = 'MX' ; 
UNIT      = 'NONE';


switch EXMX
    case 'EX'
        %==================================================================================
        Input_ex    = struct(   'I_k',   [0 1 10 100 1000 10000], 'N_k',   [10 10 10 10 10], ...
                                'I_phi', [-pi pi],                'N_phi',  50 ) ; 
        %==================================================================================
        [X , hw]    = Exziton_EX(Input_ex,'Spektrum',  'Coulomb',-700:50) ; 
        [EW, ~ ]    = Exziton_EX(Input_ex,'Eigenwerte','Coulomb',-700:50) ; 

        plot_i      = plot(hw, imag(X)/max(imag(X)), 'r',                'linewidth', 2.5) ; hold on
        stem_i      = stem(EW(EW<0), ones(1, length(EW(EW<0))), 'r--',   'linewidth', 1.3, ... 
                            'markeredgecolor', 'none') ;
        
        xlim([min(hw) max(hw)]);    xlabel('$E-E_G$ in meV', 'interpreter','latex');
        ylim([0 1]);                ylabel('Im$(\chi)$'    , 'interpreter','latex');
        grid on;                    set(gca,'GridLineStyle', '--', 'fontsize',13);
    case 'MX'
        %==================================
        n           = 300; 
        lambda      = [1 2 4 8 16 25]  ;  
        phi         = linspace(-4000,100,700) ;
        Potential   = 'Coulomb' ; 
        Method      = 'Ana' ; 
        EB3D        = 300; 
        %==================================
        
        for i=1:3
        ax          = subplot(3,2,i);   
        Input_mx    = {n,lambda(i),phi,Potential,Method};
        
        
        
        [X , hw]    = Exziton_MX(Input_mx{:},'Spektrum') ;
        [EW, ~ ]    = Exziton_MX(Input_mx{:},'Eigenwerte') ;
        
        if strcmp(UNIT,'Si')
            hw = hw*EB3D;
            EW = EW*EB3D; 
        end
        
        
        plot_i      = plot(ax,hw, imag(X)/max(imag(X)), 'b',    'linewidth', 2.0) ; hold on
        stem_i      = stem(ax,EW, ones(1, length(EW)), 'b--',   'linewidth', 1.1, ...
                            'markeredgecolor', 'none') ;
        
        xlim([min(hw) max(hw)]); ylim([0 1.1]); 
        
        if i<5
            %set(gca,'XTickLabel',[]);
        elseif strcmp(UNIT,'SI')
            xlabel(ax, '$(\hbar \rm{\omega}- E_G)$ in meV', 'interpreter','latex')
        else
            xlabel(ax, '$(\hbar \rm{\omega}- E_G)/E_B$', 'interpreter','latex')
        end
        
        if mod(i,2)==0
            %set(gca,'YTickLabel',[]);
        else
            ylabel(ax, 'Im($\rm{\chi}$)', 'interpreter','latex')
        end
        
        grid(ax,'on'); set(gca,'GridLineStyle', '--', 'fontsize',13);
        legend(plot_i,['\lambda = ' num2str(lambda(i))],'interpreter','latex');
        end
end
set(gcf, 'Position', [100, 100, 750, 900])
