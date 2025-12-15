a=zeros(10,10);
b=zeros(10,10);
for i=1:10
    for j=1:10
        a(i,j)=Heaviside(i-j);
        b(i,j)=intHeaviside(i-j);
    end
end
