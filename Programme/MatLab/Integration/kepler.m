function out= kepler(funk,a,b,N)

% For Schleife
step =  (b-a)/N ;
out  = 0; 

for i=1:N
    x_left  = a + step*i -step ; 
    x_right = a + step* i ; 
    x_mid   = a + step*i -step/2 ;
    out  = out + step*(funk(x_left)+4*funk(x_mid)+funk(x_right))/6 ;
end 

% Vektoriell

x_i = a:step/2:b;
g_i = 1/6 * [1 2*(mod(1:2*N-1,2)+1) 1];

out = (b-a)/N*dot(funk(x_i),g_i);

  figure
  plot(x_i,funk(x_i),'bx-')
  legend('Integration mit Kepler')
