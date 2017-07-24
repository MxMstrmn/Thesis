function [k,FT_fcn] = FT(I_r, N_r, I_k, N_k, fcn)

[r,g_r] = integrate(I_r,N_r,4);
[k] = integrate(I_k,N_k,4);

[K,R] = meshgrid(k,r);

weight = @(k,r) 4*pi/(2*pi)^(3/2)*r ./k .*sin(k.*r);
integrand = @(k,r) fcn(r).*weight(k,r);

FT_fcn = (integrand(K,R))'*g_r;
if FT_fcn(1)<0
    FT_fcn = -FT_fcn; 
end