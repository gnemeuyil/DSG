function theta = vectang(a,b)
[m1,n1]=size(a);
[m2,n2]=size(b);
pi = 3.1415926;

if (m1==m2) && (n1==1) &&(n2==1)
    s = sum(a.*b)/norm(a)/norm(b);
    theta = acos(s)/pi*180;
elseif (n1==n2) && (m1==1) &&(m2==1)
    s = sum(a.*b,2)/norm(a)/norm(b);
    theta = acos(s)/pi*180;
else
    error('Input should be vectors with equal size');
end