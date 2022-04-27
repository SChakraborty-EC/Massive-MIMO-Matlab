function x = MMSE_init(par,H,y,N0)
switch (par.MMSE_init) % select algorithms
    case 'ZF', % zero-forcing detection
        xhat = ZF(par,H,y);        
    case 'bMMSE', % biased MMSE detector
        xhat = real_bMMSE(par,H,y,N0);        
    case 'ML', % ML detection using sphere decoding
        xhat = ML(par,H,y);
    case 'SIMO' % SIMO lower bound
        xhat = SIMO(par,H,y,N0,s);               
    case 'Gauss-Seidel' % Gauss-Seidel detector
        xhat = real_Gauss_Seidel(par,H,y,N0);        
    case 'ADMIN' % ADMM-based Infinity Norm detector
        xhat = real_ADMIN(par,H,y,N0);        
    case 'OCDBOX' % co-ordinate descent (optimized) detector
        xhat = real_OCDBOX(par,H,y);        
    case 'Neumann' % coordinate descent
        xhat = real_Neumann(par,H,y,N0);        
    case 'Conjugate-Gradient' % conjugate gradient detector
        xhat = real_CG(par,H,y,N0);        
    case 'richardson' % conjugate gradient detector
        xhat = real_richardson(par,H,y,N0);        
    case 'SOR' % conjugate gradient detector
        xhat = real_SOR(par,H,y,N0);        
    case 'Newton_iteration' % conjugate gradient detector
        xhat = real_Newton_iteration(par,H,y,N0);        
    case 'jacobi' % conjugate gradient detector
        xhat = real_jacobi(par,H,y,N0);  
    case 'chebyshev' % conjugate gradient detector
        xhat = real_chebyshev(par,H,y,N0); 
    case 'Newton_richardson' % conjugate gradient detector
        xhat = real_Newton_richardson(par,H,y,N0);       
    case 'Newton_chebyshev' % conjugate gradient detector
        xhat = real_Newton_chebyshev(par,H,y,N0);     
    case 'Newton_GS' % conjugate gradient detector
        xhat = real_Newton_GS(par,H,y,N0);  
    case 'MMSEapprox', % biased MMSE detector
        xhat = real_MMSE_approx(par,H,y, N0);        
    otherwise,
        error('par.detector type not defined.')
end


[~,idxhat] = min(abs(xhat*ones(1,length(par.RealSymbols))-ones(2*par.MT,1)*par.RealSymbols).^2,[],2);
x = par.RealSymbols(idxhat)';

% [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% x = par.symbols(idxhat)';
end


