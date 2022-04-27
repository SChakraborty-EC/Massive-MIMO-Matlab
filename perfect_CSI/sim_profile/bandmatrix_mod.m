function G=bandmatrix_mod(Data,k)
n = length(Data);
G = zeros(n,n);
for i=1:n
    for j=1:n
        if abs(j-i) <=k
            G(i,j)=Data(i,j);
        elseif ((j==i+ k+1) && (mod(i,2) ==1)) ||((i==j+ k+1) && (mod(j,2) ==1))
            G(i,j)=Data(i,j);
        end
    end
end