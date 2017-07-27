function out = gamma_ln(x)
a = 0 ; 
threshold = 150;
if max(max(x))>threshold 
    while max(max(x))>threshold
        a = a+log(x-1) ;
        x = x-1 ;
        %disp(x)
    end
    out = a+log(gamma(x));
else
    out = log(gamma(x)); 
end
