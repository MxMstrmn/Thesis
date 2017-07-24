function vec = vecspace(VAL, n)
    dim = length(VAL(:,1));
    I   = length(n);
    vec = zeros(dim,sum(n)-(I-1));
    k=1;
    for i=1:I
        for j=1:dim
            vec(j,k:k+n(i)-1) = linspace(VAL(j,i),VAL(j,i+1),n(i));
        end
        k = k+n(i)-1;
    end
end