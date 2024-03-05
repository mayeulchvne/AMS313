function resnorm2 = RB_compute_resnorm2_stabilized(mu, Xrb,Rmat, f_parallel, f_perp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RB_compute_resnorm2_stabilized :
% calculer norme du residu stabilise au carre 
%    
% INPUT * mu: valeur parametre
%       * Xrb: solution reduite pour cette valeur du parametre
%       * Rmat: matrice N*Qa x N*Qa 
%       * f_parallel:  matrice N*Qa x Ql
%       * f_perp:  matrice Ql x Ql 
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
N  = length(Xrb);

% Calcul du vecteur c(mu), t.q. X^{-1}A(mu)UUrb(mu)= Qmat*Rmat*c(mu)
% ------------------------------------
c = zeros(N*Qa, 1);
for n=1:N
    c( ((n-1)*Qa +1):(n*Qa) ) = Xrb(n)*theta_a;
end

% Contribution provenant de Range(Q)
% ------------------------------------
vec= Rmat*c - f_parallel*theta_l;
resnorm2= vec'*vec;

% Contribution provenant de Range(Q) orthogonal
% ------------------------------------
resnorm2= resnorm2 + theta_l'*f_perp*theta_l;


end

