function [ Xrb ] = RB_solve(mu, Arb_decomp, Lrb_decomp)

% nombre de termes decomp affines 
Qa = length(Arb_decomp);
Ql = length(Lrb_decomp);

% taille de l'espace d'approximation Haute fidelite
NbDof = size(Arb_decomp{1},1);

Arb = sparse(NbDof,NbDof);
theta = PARAMETRIC_thetaA(mu);

for q=1:Qa
    Arb = Arb + theta(q)*Arb_decomp{q};
end

Lrb = sparse(NbDof,1);
theta = PARAMETRIC_thetaA(mu);

for q=1:Ql
    Lrb = Lrb + theta(q)*Lrb_decomp{q};
end

Xrb = Arb \ Lrb;


end