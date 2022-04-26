rng('shuffle')

clear conna
clear x
clear y

nrow=6;
ncol=6;

for j=1:nrow
	for k=1:ncol
		x(j,k)=2*(j-1)+rand(1);
                y(j,k)=2*(k-1)+rand(1);
end
end

voronoi(x,y)

ax=reshape(x,nrow*ncol,1);
ay=reshape(y,nrow*ncol,1);

X=[ax ay];

dt=delaunayTriangulation(X);
[V,R]=voronoiDiagram(dt);


fv=zeros(length(V)-1,1);
fx=fv;
fy=0;

for j=1:ncol-2

for k=1:nrow-2

   n=nrow*j+k+1;
   ln=length(R{n});

for i=1:ln-1
	a=R{n}(i);b=R{n}(i+1);
        conna(a,b)=1;
        conna(b,a)=1;
        fv(a)=1;fv(b)=1;
end

a=R{n}(ln);b=R{n}(1);
conna(a,b)=1;
conna(b,a)=1;
fv(a)=1;fv(b)=1;

end
end

for k=2:length(fv)
	fx(k)=fv(k)*(fy+1);
        fy=fy+fv(k);
end

ca=zeros(fy);

for j=1:length(fx)
	for k=1:length(fx)
	if (fx(j) > 0) && (fx(k) > 0)
			  ca(fx(j),fx(k))=conna(j,k);
end 
end 
end

Va=zeros(fy,2);

for j=1:length(fx)
	if (fx(j) > 0)
	  Va(fx(j),1)=V(j,1);
Va(fx(j),2)=V(j,2);

end
end

hold off

vorinfo

network
