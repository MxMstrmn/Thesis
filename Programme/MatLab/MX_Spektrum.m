mu          = (0.41*const.e*0.46)/(0.41+0.46) ;             % m*e = 0.46me & m*h = 0.41mh
EB          = h^2/(2*mu*const.a0^2) ;
B           = 3*1e3*const.e ;
lambda      = const.a0^2*const.e*B/h ;
w_c         = (const.e *B) /mu ;                            % Zyklotronfrequenz des Exzitons