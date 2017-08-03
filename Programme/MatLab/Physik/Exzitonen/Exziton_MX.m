function [out1, out2] = Exziton_MX(n,lambda,phi,Potential,Method,Output) 

disp(['Berechnung: ' Output ' des ' Potential '-Potentials'])
dim                 = n+1 ; 

% lambda           = a0^2*B*c.e/c.hbar*unitB;

c               = constants(); 
unitB           = 1/c.e*1e-3 ;   
B               = 150*4*lambda * unitB ; 
disp(B/unitB)
% MoS2
mu              = 0.46*0.41/0.87*c.me ;
a0              = 1.1016 ; 
% % GaAs
% mu              = 0.457*0.067/0.524*c.me ;  
% a0              = 12.3 ;

l               = sqrt(c.hbar/c.e/B) ; 
lambda          = a0^2/l^2 
EB              = c.hbar.^2 /2 / mu / a0^2 ;
% dimensionless
CONS            = EB*(2*lambda/pi)^0.5 ; 
% SI 
 
wc              =  c.e*B/mu ;

switch Potential
    case 'Coulomb'
            %CONS        =  c.e^2 /4 /pi /c.eps0 /c.eps /sqrt(2) *sqrt(c.e*B/c.hbar);
        switch Method
            case 'Ana'
                %CONS    = CONS/sqrt(pi) ;
                V_ij    = - CONS *exp(gammaPrefactor(0:n)) .*F32(0:n) ;
            case 'Num'     
                try
                    V_ij        = - CONS *csvread(['VC_ij_' num2str(n) '.dat']);
                    
                catch
                    disp('Berechne Coulomb-Matrix erst mit gaussLaguerre.m')
                    Output = 'issMatrix' ; 
                end
                % Vergleich
                % imagesc(exp(gammaPrefactor(0:20)) .*F32(0:20)./V_ij_Num(0:20))
        end 
 

    case 'Keldysh'
        CONS            =  c.e^2 /4 /pi /c.eps0 /sqrt(2) *sqrt(c.e*B/c.hbar);
        try 
            V_ij        =  - CONS *csvread(['VK_ij_' num2str(n) '.dat']);
            
        catch
            disp('Berechne Keldysh-Matrix erst mit gaussLaguerre.m')
            Output = 'issMatrix' ;
        end
end

Hmx_ii      = @(n)      diag     (  0.5*c.hbar*wc*(2*(0:n)+1) ) ;
Inh_ii      = @(n,phi)  eye(dim)*( -phi    -1i*20     ) ;      



switch Output
    case 'Spektrum'
        H           = Hmx_ii(n) + V_ij ;
        b           = ones (     dim    ,1) ;
        X           = zeros(length(phi) ,1) ;
        
        for i=1:length(phi)
            A       = H + Inh_ii(n,phi(i)) ;
            x_n     = linsolve(A,b) ;
            X(i)    = sum(x_n)*lambda/pi ;
        end
        out1        = X ;
        out2        = phi ;
    case 'Eigenwerte'
        H              = Hmx_ii(n) + V_ij ;
        [states, EW]   = eig(H,'vector');
        [EW, idx]      = sort(EW);
        states         = states(:,idx);
        EW             = EW(EW<5000) ;
        states         = states(:,EW<5000) ;
        [nn,g_nn]      = integrate([0 dim],dim,4) ; 
        norm           = sqrt(2*pi*(states.^2)'*(nn.*g_nn));
        for i=1:length(EW); states(:,i) = states(:,i)*1/norm(i)*sign(states(1,i)); end
         
        % Anzeigen der Grundzustandsenergie
%         disp(EW(1))
        out1            = EW ;
        out2            = states ;
    case 'issMatrix'
        disp('Programm  gestoppt')        
        out1            = [] ;
        out2            = [] ;
end



