function price = GBM_Price(K,x,s,r,Lstar)
    
    if x <= Lstar
        price = max(K-x,0);
    else
        price = (K-Lstar)*(x/Lstar)^(-2*r/(s*s));
    end
end