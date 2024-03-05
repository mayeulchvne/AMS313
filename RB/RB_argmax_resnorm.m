function [ maxresnorm, mu_star ] = RB_argmax_resnorm(Xi, Arb_decomp, Lrb_decomp, ...
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

error('RB_argmax_resnorm() not yet implemented');

end

