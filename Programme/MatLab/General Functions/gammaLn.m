function out = gammaLn(x)
a = 0 ; 
threshold = 150;
if max(max(x))>threshold 
    while max(max(x))>threshold
        a = a+log(x-1) ;
        x = x-1 ;
    end
    out = a+log(gamma(x));
else
    out = log(gamma(x)); 
end
