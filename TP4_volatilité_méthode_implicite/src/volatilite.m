clear

function volatilite=volatilite(S,E,r,T,C0)
 
  d1=@(sigma) (log(S./E) + (r+sigma.^2).*T)./(sigma.*sqrt(T));
  d2=@(sigma) d1(sigma) - sigma .* sqrt(T);

  Ee = E.*exp(-r .* T);
  N1=@(sigma) 0.5*erfc(-d1(sigma)./sqrt(2));
  N2=@(sigma) 0.5*erfc(-d2(sigma)./sqrt(2));
  c=@(sigma) S.* N1(sigma) - Ee.*N2(sigma);
  F=@(sigma) C0-c(sigma);

  a=0; a1=a;
  b=10; b1=b;
  d=(a+b)/2;
  i=0;
    
    
  while i<300
    if F(a)*F(d)>0
      a=d; 
      else b=d;
      end
  d=(a+b)/2;
  i=i+1;  
  end

  plot([a1 b1], F([a1 b1]))
end





volatilite=d;