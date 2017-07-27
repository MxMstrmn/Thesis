function A = gammaPrefactor(n)

A = zeros(max(n)+1) ; 
try
    A = csvread(['MatrixFiles/Gamma_' num2str(max(n)) '.dat']);
%     disp('Vorfaktor wurde schon berechnet')
catch
    disp('Vorfaktor muss noch berechnet werden')
    for ii=n
        for jj=0:ii
            A(jj+1,ii+1) = gammaLn(ii-jj+0.5) + gammaLn(jj+0.5) - gammaLn(ii-jj+1) - gammaLn(jj+1);
            A(ii+1,jj+1) = A(jj+1,ii+1);
        end
    end
    if isempty(A(isnan(A))) == 0; disp('Dimension zu gross'); end 
    csvwrite(['MatrixFiles/Gamma_' num2str(max(n)) '.dat'], A)
    disp('Vorfaktor wurde berechnet')
end
end



    




                   