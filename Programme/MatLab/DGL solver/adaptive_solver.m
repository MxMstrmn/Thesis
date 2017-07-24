function [t,y] = adaptive_solver(RHS, y0,t0,tend, delta)
% Dormand Prince
c = [0; 1/5; 3/10; 4/5; 8/9; 1; 1]; 

A = [0 0 0 0 0 0 0; 
    1/5 0 0 0 0 0 0; 
    3/40 9/40 0 0 0 0 0;
    44/45 -56/15 32/9 0 0 0 0; 
    19372/6561 -25360/2187 64448/6561 -212/729 0 0 0;
    9017/3168 -355/33 46732/5247 49/176 -5103/18656 0 0;
    35/384 0 500/1113 125/192 -2187/6784 11/84 0];

b1 = [35/384 0 500/1113 125/192 -2187/6784 11/84 0];
b2 = [5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40];

y(:,1) = y0;
t(1) = t0;

maxdelta = 100*delta;
mindelta = 0.01*delta; 
increase = 1.3;
decrease = 0.5;
precision = 1e-7; 

maxIter = 100;
Iter = 0;
i = 2; 

while t(i-1)<tend
    k = zeros(length(y0),length(c));
    for j=1:length(c)
        k(:,j) = RHS(t(i-1)+delta*c(j), y(:,i-1)+delta*(A(j,:)*k')');
    end
    
    err = max((b1-b2)*k');
    if err < precision
        y(:,i) = y(:,i-1) + delta*(b1*k')' ;
        t(i) = t(i-1)+delta;
        
        delta = increase*delta;
        if delta > maxdelta
            delta = maxdelta;
        end
        i = i+1;
        Iter = 0; 
    else
        delta = decrease*delta;
        if delta < mindelta
            delta = mindelta;
        end
    end
    Iter = Iter+1;
    if Iter > maxIter
        disp('Delta zu groß gewählt')
        break
    end
end
 

y(:,i:end) = [];
t(i:end) = [];  
    
end
