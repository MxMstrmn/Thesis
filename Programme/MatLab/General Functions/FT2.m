function [k,FT_k,g_k] = FT2(I_r, N_r, I_k, N_k, fcn)

% Die Funktion wird in den k raum transformiert, weshalb die Gewichte nicht
% benötigt werden. Da sich die phi-Integration im Gegensatz zum 3D Fall
% nicht ausführen lässt, müssen hier 2 Integrale ausgeführt werden, also
% wird ein 3D Gitter gebraucht. 
[k, g_k]             = integrate(I_k,N_k,4) ;
[r, g_r]        = integrate(I_r,N_r,4) ;
[phi, g_phi]    = integrate([0 2*pi],400,4) ;

[R,K,PHI] = meshgrid(r,k,phi);

% Integrand bestimmen 
% In weight sind alle Teile des Integranden, die fest sind, also Normerhaltender
% Vorfaktor, die Jacobideterminante und die e-Funktion in Polarkoordinaten
weight = @(k,r,phi) 1/(2*pi)*r .*exp(1j*k.*r.*sin(phi));
integrand = @(k,r,phi) fcn(r).*weight(k,r,phi);

% Fouriertransformierte bestimmen
% Zuerst noch HIlfsparameter bestimmen, dann die phi Summation durchführen,
% dafür 3D-Gewichtsmatrix weight_phi, und dann noch die r Integration
% mitsamt der Gewichte. 
dim_k = length(k);
dim_r = length(r);
weight_phi = permute(repmat(g_phi',[dim_k,1,dim_r]),[1,3,2]);
FT_k = sum((integrand(K,R,PHI).*weight_phi),3)*g_r;

if FT_k(1)<0
    FT_k = -FT_k; 
end

FT_k = FT_k/sqrt(normalize(k,FT_k,g_k,2));

% figure
% plot(k,FT_k)
