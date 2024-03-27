function [ maxresnorm, mu_star ] = RB_argmax_resnorm_stabilized(Xi, Arb_decomp, Lrb_decomp, Rmat, f_parallel, f_perp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RB_argmax_resnorm_stabilized :
% trouver le parametre qui maximise la norme du residu stabilise
%          
% INPUT * Xi: espace des parametres (taille n_train x 2)
%       * Arb_decomp: cellarray des matrices de l'operateur reduit (Qa matrices de taille N x N)
%       * Lrb_decomp: cellarray des vecteurs du RHS reduit (Ql vecteurs taille N x 1)
%       * Rmat: matrice N*Qa x N*Qa pour le residu stabilise
%       * f_parallel:  matrice N*Qa x Ql pour le residu stabilise
%       * f_perp:  matrice Ql x Ql pour le residu stabilise
%
% OUTPUT - maxresnorm: maximum (sur Xi) de la norme du residu 
%        - mu_star (taille 2):  parametre maximisant
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% nombre de points espace des parametres
nb_points = size(Xi, 1);

mu_star = Xi(1,:);
maxresnorm = 0;

NumberOfNegRes = 0;

for i=1:nb_points
        mu = Xi(i,:);
        %calcul solution base réduite
        Xrb = RB_solve(mu, Arb_decomp, Lrb_decomp);
        %résidu calcul efficace
        res2 = RB_compute_resnorm2_stabilized(mu, Xrb,Rmat, f_parallel, f_perp);
        
        if (res2 > maxresnorm)
            maxresnorm = sqrt(res2);
            mu_star = mu;
        elseif (res2<0)
            NumberOfNegRes = NumberOfNegRes+1;
        end
end

end

