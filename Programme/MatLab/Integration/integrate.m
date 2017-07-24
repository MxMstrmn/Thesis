function [x,g] = integrate(I, N, method, fcn)
    % I     :   Intervalle des bestimmten Integrals
    % n     :   Anzahl der "Intervalle" pro Intervall (? Funktionsaufrufen)
    % fcn   :   Integrand, anonyme Funktion
    % method:   Gibt Integrationsmethode an 
    
    if nargin < 4
        fcn = @(x) 0*x;
    end
    method_names = {'Rechteck', 'Trapez', 'Kepler', 'GTS Methode', 'Gauß-Chebyshev'};
    y = zeros(length(N),1);     % function value
    x = []  ;      % x values
    g = []  ;      % weighting
    
    for i=1:length(N)
        a = I(i);
        b = I(i+1); 
        n = N(i);
        switch method
            case {1} %Rechteck
                step = (b-a)/n ;
                x_i = a+step/2 : step : b-step/2 ;
                g_i = step *ones(n,1) ;
                             
            case {2} %Trapez
                step = (b-a)/(n-1) ;
                x_i = a : step : b ;
                g_i = step *[1/2; ones(n-2,1); 1/2] ;
                           
            case {3} %Kepler
                % N muss ungerade sein und >2
                step = 2*(b-a)/(n-1) ;
                x_i = a : step/2 : b ;
                g_i = 1/6*step *[1 ; 2*(mod(1:n-2,2)+1)' ; 1] ;
                           
            case {4} %GTS Methode
                phi_step = 2*pi/(n+1) ;
                phi_i = phi_step*(1:n) ;
                x_i = (b-a)/(2*pi) *(phi_i-sin(phi_i)) +a ;
                g_i = (b-a)/(2*pi) *phi_step *(1 - cos(phi_i)') ;
                    
            case {6} %GTS Methode
                phi_step = 2*pi/(n+1) ;
                phi_i = phi_step*(1:n) ;                
                a1 = (3*a-b)/2;
                c  = (a+b)/2 ; 
                b1 = (3*b-a)/2;
                x_1 = (c-a1)/(2*pi) *(phi_i-sin(phi_i)) +a1 ;
                g_1 = (c-a1)/(2*pi) *phi_step *(1 - cos(phi_i)') ;      
                x_2 = (b1-c)/(2*pi) *(phi_i-sin(phi_i)) +c ;
                g_2 = (b1-c)/(2*pi) *phi_step *(1 - cos(phi_i)') ;  
                x_i = [x_1  x_2];
                g_i = [g_1  g_2]; 
                g_i = g_i(x_i>a & x_i<b);
                x_i = x_i(x_i>a & x_i<b);
            
            otherwise %Gauss-Chebyshev
                phi_step = pi/n;
                phi_i = phi_step*(1:n)-phi_step/2 ;
                x_i = (b-a)/2 *(cos(phi_i)+1) + a ;
                g_i = (b-a)/2 *phi_step *sin(phi_i)' ;

                
        end
        % y(i) = fcn(x_i)*g_i;
        % merge weighting vectors
        x = [x ; x_i(:)];
        g = [g ; g_i(:)];
    end
    out = sum(y(:));
    % KLEINERE AUSWERTUNGEN 
    %sprintf('Die %s-Methode ergibt: I=%f', method_names{method},out)
    %figure 
    %plot(x,g,'o');
end
