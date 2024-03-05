function resnorm2 = RB_compute_resnorm2(mu, Xrb, Respart_ll, Respart_la, Respart_aa)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RB_compute_resnorm2 :
% calculer norme du residu RB au carre 
%        
% INPUT * mu (taille 2): valeur parametre
%       * Xrb (taille N): solution reduite pour cette valeur du parametre
%       * Respart_ll: cellarray Ql x Ql (de scalaires) des interactions L-L
%       * Respart_la: cellarray Ql x Qa des vecteurs (taille N x 1) des interations L-A
%       * Respart_aa: cellarray Qa x Qa des matrices (taille N x N) des interations A-A
%
% OUTPUT - resnorm2: valeur de la norme du residu au carre
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% valeur des coefficients
% ------------------------
theta_a = PARAMETRIC_thetaA(mu);
theta_l = PARAMETRIC_thetaL(mu);

Qa = length(theta_a);
Ql = length(theta_l);

% Contribution des interaction RHS-RHS
% ------------------------------------
resnorm2 = 0.0;


% Contribution des interaction RHS-LHS
% ------------------------------------


% Contribution des interaction LHS-LHS
% ------------------------------------

error('RB_compute_resnorm2() not yet implemented');


end

