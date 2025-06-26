function F=yexL(x3)
global f1total;
global K;
global L;
global c;
global CL;
global f2;
for i=1:K
  f(i)=(L(i)*c^2*CL(i)+L(i)*c*x3(i))/(x3(i)*f2(i)+c*CL(i)*(x3(i)+f2(i)))-x3(K+1);
end
i=1:K;
f(K+1)=sum(x3(i))-f1total;
F=f;
end
