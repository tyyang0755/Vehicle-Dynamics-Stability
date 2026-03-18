function [F,G] = predmat(A,B,C,Np)
% Predmat build the prediction matrix F and G offline 
% Y_N = Fx_0 + GU_N
    nx = size(A,1);      % number of states
    nu = size(B,2);      % number of inputs
    ny = size(C,1);      % number of outputs
    
    F = zeros(ny*Np, nx);
    G = zeros(ny*Np, nu*Np);
    
    for i = 1:Np
        
        % F matrix block
        F((i-1)*ny+1:i*ny,:) = C*A^i;
        
        % G matrix blocks
        for j = 1:i
            
            G((i-1)*ny+1:i*ny , (j-1)*nu+1:j*nu) = C*A^(i-j)*B;
            
        end
        
    end


end