function out = gauss(funk,a,b,N) 
    phi_i = pi/(2*N):pi/N:(2*N-1)*pi/(2*N); 
    x_i = (b-a)/2 *(cos(phi_i)+1) +a ; 
    g_i = (b-a)/2 * sin(phi_i) *pi/N ; 
    
    out = funk(x_i)*g_i';
    
    figure
    plot(x_i,funk(x_i),'bx-')
    legend('Integration mit Gauß')
end