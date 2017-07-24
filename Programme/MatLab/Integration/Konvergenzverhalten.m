fcn  = @(x) x.^2 ;
fcn2 = @(x) sin(x) ; 
y = zeros(10,5) ;
for i=1:10
    for j=1:5
        y(i,j)=integrate(fcn2,0,pi,2*i+1,j);
    end
end
N=3:2:21;
plot(N,y(:,1),N,y(:,2),N,y(:,3),N,y(:,4),N,y(:,5)) ; 
legend('Rechteck','Trapez','Kepler','Gauß-Cheb.','GTS','Location','southeast')
