function A = gamma_prefactor(n)

A = zeros(max(n)+1) ; 
try
    A = csvread(['MatrixFiles/Gamma_' num2str(max(n)) '.dat']);
%     disp('Vorfaktor wurde schon berechnet')
catch
    disp('Vorfaktor muss noch berechnet werden')
    for ii=n
        for jj=0:ii
            A(jj+1,ii+1) = gamma_ln(ii-jj+0.5) + gamma_ln(jj+0.5) - gamma_ln(ii-jj+1) - gamma_ln(jj+1);
            A(ii+1,jj+1) = A(jj+1,ii+1);
        end
    end
    if isempty(A(isnan(A))) == 0; disp('Dimension zu groﬂ'); end 
    csvwrite(['MatrixFiles/Gamma_' num2str(max(n)) '.dat'], A)
    disp('Vorfaktor wurde berechnet')
end
end



    




                   