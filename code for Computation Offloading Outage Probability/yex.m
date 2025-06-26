function F=yex(x3)
global newk1;
global f1total;
global K;
global alpha;
global L;
global B;
global thre;
global c;
global dd0;
global v;
for i=1:K
  m(i)=2^(x3(i)*alpha(i)*L(i)/(B*(thre*x3(i)-alpha(i)*L(i)*c)))-1;
  f(i)=1-marcumq(newk1(i),sqrt(2*m(i)/(dd0(i)*v(i))))-x3(K+1);
end
i=1:K;
f(K+1)=sum(x3(i))-f1total;
F=f;
end
