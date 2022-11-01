function [WHS_VaR2] = WHSVaR2(Returns,pVaR,eda)
       N = length(Returns);
       tau = sort(1:N,'descend');
    
       n_tau = ((eda.^(tau-1)*(1-eda))/(1-eda^N))'; 
       d1 = [Returns,n_tau];
       y1 = sortrows(d1(:,:),1);
       c1 = cumsum(y1(:,2),1); 
       y2 = [y1,c1];
       
       index = y2(:,3)<=pVaR;
       index2 =y2(index,:);
       if  isempty(index2)
        index3 = [0,0,0];
       else
        index3 = index2(end,:);
       end
     c_we = index3(:,3);

       WHS_VaR2 = -quantile(Returns,c_we);    
end