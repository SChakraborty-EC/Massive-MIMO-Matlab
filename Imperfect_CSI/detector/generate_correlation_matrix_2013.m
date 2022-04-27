function [Rtx, Rrx] = generate_correlation_matrix_2013(par)
Rtx = eye(par.MT, par.MT);
Rrx = eye(par.MR, par.MR);
if par.correlation == 1
for p = 1:par.MT
    for q  = 1:par.MT
        if p<q
            Rtx(p,q) = ( par.zeta*exp(1i*par.theta))^(q-p);
            Rtx(q,p) = conj(Rtx(p,q));
        end
    end
end
for p = 1:par.MR
    for q  = 1:par.MR
        if p<q
            Rrx(p,q) = ( par.zeta*exp(1i*par.theta))^(q-p);
            Rrx(q,p) = conj(Rrx(p,q));
        end
    end
end
end