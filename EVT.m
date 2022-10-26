function [VAR] = EVT(u,beta,shape,N,Nu,pVaR)
    VAR = -u+beta/shape*((pVaR*N/Nu)^(-shape)-1);
end