function V = V_ij_Num(n)
V = zeros(length(n));
for jj = n
    for ii = 0:jj
       L    = @(x,n,a) nchoosek(n+a,n).*hypgeomNum(-n,a+1,x) ; 
       
       f    = @(x,ii,jj) 1./sqrt(2*x).*exp(-x).*(x).^(jj-ii).*(L(x,ii,jj-ii)).^2;
       
%        [k,g_k]=integrate([0 1000],400,4); 
%        V(ii+1,jj+1) = factorial(ii)/factorial(jj)*dot(f(k,ii,jj),g_k) ; 
       V(ii+1,jj+1) = sqrt(2*pi)*factorial(ii)/factorial(jj)*quadgk(@(x)f(x,ii,jj),0,5000);
       
       V(jj+1,ii+1) = V(ii+1,jj+1) ;
    end
%     disp(jj)
end

csvwrite(['MatrixFiles/V_ij_' num2str(max(n)) '.dat'], V)
