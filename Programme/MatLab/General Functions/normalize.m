function norm =  normalize(x,y, weight,dim)
y     = y(:);
weight  = weight(:); 

% in 2D 
switch dim
    case {2}
        norm    = 2*pi*(x.*y.*conj(y))'*weight;
end

% Wenn die Funktion erweitert werden soll, dann einfach einen switch case
% einbauen 