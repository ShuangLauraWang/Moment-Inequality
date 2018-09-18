function A=PossibleOutcome(n,m)
A=ones(n^m,n);
a=(1:n)';

for i=1:m
A(:,i)=repmat(kron(a,ones(n^(i-1),1)),n^m/n^i,1);
end

end