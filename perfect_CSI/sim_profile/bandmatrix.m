function G=bandmatrix(Data,k)
n = length(Data);
G = zeros(n,n);
for i=1:n
    for j=1:n
        if abs(j-i) <=k
            G(i,j)=Data(i,j);
        end
    end
end