function V_ij = gaussLaguerre(dim, Potential)

% % Beispiel
% 
%   dim         = 20 ;
%   Potential   = 'Keldysh' ; 
% 
% % 

formatSpec  = '%f';


switch Potential 
    case {'Coulomb'}
        
        fcn     = 'fcn_Coulomb';
        f       = @(x) fcn_Coulomb(x) ;
        
%         Für den fall das Keldysh und Coulomb mit unterschiedlich vielen 
%         Stützstellen konvergien 
        order   = @(n) n+1 ;
    
    case {'Keldysh'}
        
        fcn     = 'fcn_Keldysh';
        f       = @(x) fcn_Keldysh(x) ; 
        order   = @(n) n+1 ;
end

V_ij       = zeros(dim+1) ;
for jj=0:dim
    for ii = 0:jj
        n           = ii;
        n1          = jj;
        setGlobaln(n)  ;
        setGlobaln1(n1) ;
        beta        = n1-n-0.5 ;
         
        disp(['Ordnung ist ' num2str(order(n)) ' bei n=' num2str(n)])
        gen_laguerre_rule(order(n),beta,0,1,fcn)
        
        H           = fscanf(fopen([fcn '_w.txt'],'r'),formatSpec);
        ak          = fscanf(fopen([fcn '_x.txt'],'r'),formatSpec);
        Int         = fscanf(fopen([fcn '_r.txt'],'r'),formatSpec);
        
        const       = sqrt(pi)*factorial(n)/factorial(n1);
        
        GAUSSLAGUERRE       = const*(H.'*f(ak));
        fprintf('\n GAUSSLAGUERRE = %0.4g \n',GAUSSLAGUERRE)
        
        V_ij(ii+1,jj+1)    = GAUSSLAGUERRE ; 
        V_ij(jj+1,ii+1)    = V_ij(ii+1,jj+1) ; 
        
        fclose all;
    end
end

end


