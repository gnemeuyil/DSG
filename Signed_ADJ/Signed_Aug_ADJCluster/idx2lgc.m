function lgc = idx2lgc(Idx)
[m,n]=size(Idx);
a= unique(Idx);
[ma na]=size(a);

lgc = zeros(m,ma);
for i = 1:ma
    lgc(:,i) = logical(Idx == a(i));
end;