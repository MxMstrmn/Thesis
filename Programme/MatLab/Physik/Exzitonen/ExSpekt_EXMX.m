% clear all
% clc 
figure
%hold on;
 
EXMX      = 'MX' ; 
UNIT      = 'SI';


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
        n           = 20; 
        LambdaOrB   = 150;%1.5*[1 2 4 6 10 16]  ;  % Depends on your input, either magnetic field B or dimensionless quantity lambda 
        phi         = linspace(-700,100,600) ;
        Potential   = 'Keldysh' ; 
        Method      = 'Num' ; 
        %==================================
        
        for i=1:min(length(LambdaOrB),6)
        ax          = subplot(3,2,i);   
        Input_mx    = {n,LambdaOrB(i),phi,Potential,Method};
        
        
        
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
        
        if strcmp(UNIT, 'SI')
            legend(plot_i,['B = ' num2str(LambdaOrB(i))],'interpreter','latex');
        else
            legend(plot_i,['\lambda = ' num2str(LambdaOrB(i))],'interpreter','latex');
        end
        
        
        grid(ax,'on'); set(gca,'GridLineStyle', '--', 'fontsize',13);
        
        end
end
set(gcf, 'Position', [100, 100, 750, 900])
