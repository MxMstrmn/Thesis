function S = constants

% Physikalische Konstanten

% Einheiten in [ nm, ps, pA, meV ]
CONST(1).c0         = 299792.458;                                   % nm / ps
CONST(1).e          = 1.602176565e5;                                % pA * ps
CONST(1).kb         = 0.08617343;                                   % meV / K
CONST(1).me         = 5.686e-3        ;                             % meV * ps^2 / nm^2
CONST(1).hbar       = 0.6582119514;                                 % meV * ps
CONST(1).eps0       = 1.41844e6;                                    % pA^2 * ps^2 / meV / nm
CONST(1).eps        = 4.508;                                        % eps_r in Einheit von eps0
CONST(1).factor     = CONST(1).e^2 / 2 / CONST(1).eps0;
CONST(1).Ry         = CONST(1).me * CONST(1).e^4 / 8 / CONST(1).eps0^2 ...
                      / (CONST(1).hbar*2*pi)^2;                     % 1 Rydberg in meV
CONST(1).a0         = 4 * pi * CONST(1).eps0 * CONST(1).hbar^2 ...
                      / CONST(1).me / CONST(1).e^2;                 % Bohrradius in nm
CONST(1).Gamma      = 10;                                           % Linienbreite in meV
CONST(1).E_G        = 0;                                            % Bandlücke des HL    

S = CONST;  
end

% Eventuell nuetzlich fuer MX gerechnet mit Einheiten
%
% const       = constants ;
% h           = const.hbar ;                                      
% B           = 1; 
% mu          = (0.46*const.me*0.41) /(0.46+0.41);                % m*e=0.46me, m*h=0.41me  
% w_c         = ( const.e*B)/mu ;                                 % Zykloronfrequnez des Exzitons
% coulomb     = (-const.e^2) /(4*pi*const.eps*const.eps0) ;       % Coulomb Vorfaktor
