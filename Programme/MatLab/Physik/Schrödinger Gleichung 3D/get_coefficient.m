function out = get_coefficient(a,b,c,d)

c1 = log(b/d)/(c/a -1 -log(a/c)) ;
c2 = c1 + log(b/(a^c1)) ; 
c3 = c1/a ; 

out = [c1 c2 c3];
end