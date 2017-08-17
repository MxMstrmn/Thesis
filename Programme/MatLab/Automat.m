B = [1.5 8 15 20 30 100] ; 

dim = 70;

for Bi=B
    gaussLaguerre(dim,'Keldysh',Bi)
end

gaussLaguerre(dim,'Coulomb',Bi)