
function out=rechteck(funk,a,b,N)

% For Schleife
step =  (b-a)/N ;
out  = 0; 

for i=1:N
    x_mid   = a + step*i -step/2 ;
    x_left  = a + step*i -step ; 
    x_right = a + step* i ; 
    
    out  = out + (b-a)/N*funk(x_mid) ;
end 

% Vektoriell

x_i = a+step/2:step:b-step/2 ;
y_i = funk(x_i) ;
out = (b-a)/N*sum(y_i);

 