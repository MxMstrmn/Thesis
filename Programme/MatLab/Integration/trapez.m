function out= trapez(funk,a,b,N)

% For Schleife
step =  (b-a)/N ;
out  = 0; 

for i=1:N
    x_left  = a + step*i -step ; 
    x_right = a + step* i ; 
    
    out  = out + step*(funk(x_left)+funk(x_right))/2 ;
end 

% Vektoriell
x_i    = a:step:b ;
g_i    = [1/2  ones(1,length(x_i)-2) 1/2];
out = (b-a)/N*dot(funk(x_i),g_i);
