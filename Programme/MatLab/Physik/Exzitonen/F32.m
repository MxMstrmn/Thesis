function A = F32(n)

A = zeros(length(n));
try
    A = csvread(['MatrixFiles/F32_' num2str(max(n)) '.dat']);
%     disp('F32 wurde schon berechnet')
catch
    disp('F32 muss noch berechnet werden')
    for ii=n
        for jj=0:ii
            % hypgeomNum() zeigt im Bereich n=0:300 eine max. Abweichung von
            % 0.05% gegen�ber hypergeom()
            A(jj+1,ii+1) = hypgeomNum([-jj ii-jj+0.5 0.5],[ii-jj+1 -jj+0.5],1);
            A(ii+1,jj+1) = A(jj+1,ii+1);
        end
%         if mod(ii,10)==0
%         disp(ii)
%         end
    end
%csvwrite(['MatrixFiles/F32_' num2str(max(n)) '.dat'], A)
disp('F32 wurde berechnet')
end


    

