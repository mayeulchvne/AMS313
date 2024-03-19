function [ Arb_decomp, Lrb_decomp ] = RB_reduced_decomp(AA_decomp, LL_decomp, PP)


Qa = length(AA_decomp);
Ql = length(LL_decomp);

Arb_decomp = cell(Qa);
Lrb_decomp = cell(Ql);

for i=1:Qa
    Arb_decomp{i} = PP'*AA_decomp{i}*PP;
end

for i=1:Ql
    Lrb_decomp{i} = PP'*LL_decomp{i};
end

end