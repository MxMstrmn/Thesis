function V = V_ij_Num(n)
try 
    V   = csvread(['VC_ij_' num2str(max(n+1)) '.dat']);
catch
    disp('Berechne Matrix mit gaussLaguerre.m')
% % %     
%     Mein erster Versuch, funktioniert nicht gut, gehe zu gaussLaguerre.m
%     V = zeros(length(n));
%     for jj = n
%         for ii = 0:jj
%             L    = @(x,n,a) nchoosek(n+a,n).*hypgeomNum(-n,a+1,x) ;
%             
%             f    = @(x,ii,jj) 1./sqrt(2*x).*exp(-x).*(x).^(jj-ii).*(L(x,ii,jj-ii)).^2;
%             
%             %         Eigener Integrator
%             %        [k,g_k]=integrate([0 1000],400,4);
%             %        V(ii+1,jj+1) = factorial(ii)/factorial(jj)*dot(f(k,ii,jj),g_k) ;
%             
%             %         Quad Integrator
%             V(ii+1,jj+1) = sqrt(2*pi)*factorial(ii)/factorial(jj)*quadgk(@(x)f(x,ii,jj),0,5000);
%             
%             V(jj+1,ii+1) = V(ii+1,jj+1) ;
%         end
%         %     disp(jj)
%     end
%     csvwrite(['MatrixFiles/V_ij_' num2str(max(n)) '.dat'], V)  
% % % 

end


