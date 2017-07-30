function F = hypgeomNum(ak,bk,z)
% Reicht für meine Anwendung und ergibt das richtige für diesen Fall, hat
% aber im Komplexen Problem und ergibt Nan, wenn hypergeom() eine komplexe
% Zahl als Ergbnis ausgibt 

F       = 1; 
k       = 1;
temp    = 1; 

while abs(temp) > 1e-20
   
   temp = temp*prod(ak)./prod(bk)/k ;
%    disp(temp)
   F    = F + temp*z.^k;
   ak   = ak+1;
   bk   = bk+1; 
   k    = k+1; 
   
end
