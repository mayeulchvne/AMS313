function [ maxresnorm, mu_star, NumberOfNegRes ] = RB_argmax_resnorm(Xi, Arb_decomp, Lrb_decomp, ...
                                       Respart_ll, Respart_la, Respart_aa)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RB_argmax_resnorm :
% trouver le parametre mu_star qui maximise le residu 
%          
% INPUT * Xi: espace des parametres (taille n_train x 2)
%       * Arb_decomp: cellarray des matrices de l'operateur reduit (Qa matrices de taille N x N)
%       * Lrb_decomp: cellarray des vecteurs du RHS reduit (Ql vecteurs taille N x 1)
%       * Respart_ll: cellarray Ql x Ql (de scalaires) des interactions L-L
%       * Respart_la:  cellarray Ql x Qa des vecteurs (taille N x 1) des interations L-A
%       * Respart_aa:  cellarray Qa x Qa des matrices (taille N x N) des interations A-A
%
% OUTPUT - maxresnorm: maximum (sur Xi) de la norme du residu 
%        - mu_star (taille 2): valeur du parametre maximisant
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nombre de points espace des parametres
n_train= size(Xi, 1);

mu_star = Xi(1,:);
maxresnorm = 0;

NumberOfNegRes = 0;

for i=1:n_train
        mu = Xi(i,:);
        %calcul solution base réduite
        Xrb = RB_solve(mu, Arb_decomp, Lrb_decomp);
        %résidu calcul efficace
        res2 = RB_compute_resnorm2(mu, Xrb, Respart_ll, Respart_la, Respart_aa);
        
        if (res2 > maxresnorm)
            maxresnorm = sqrt(res2);
            mu_star = mu;
        elseif (res2<0)
            NumberOfNegRes = NumberOfNegRes+1;
        end
end

end

