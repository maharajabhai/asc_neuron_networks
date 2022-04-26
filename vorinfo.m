nodes=length(Va);
tot=nodes+2;

diam0=8;

vin=nodes+1;
vout=nodes+2;

conn=zeros(tot,tot);
a=zeros(nodes,nodes);
b=zeros(nodes,1);
alpha=zeros(tot,tot);
dist=zeros(tot,tot);
poss=zeros(tot,2);
nodepos=zeros(tot,3);
q=zeros(tot,tot);
press=zeros(tot,1);

[ymin,Imin]=min(Va(:,2));
[ymax,Imax]=max(Va(:,2));

	    xmin=Va(Imin,1);
	    xmax=Va(Imax,1);

posio=[xmax 2*ncol; xmin 0];

poss=10*[Va;posio];

for j=1:length(ca)
	for k=1:length(ca)
		conn(j,k)=ca(j,k);
end
end

conn(Imax,nodes+1)=1;
conn(Imin,nodes+2)=1;
conn(nodes+1,Imax)=1;
conn(nodes+2,Imin)=1;

pio=[50;35];

visc=3;
aflow=pi/(128*visc);

for j=1:tot
	for k=1:tot
		ax=poss(j,1);ay=poss(j,2);
                bx=poss(k,1);by=poss(k,2);

                if conn(j,k) > 0
                dist(j,k)= ( (ax-bx)^2+(ay-by)^2)^(1/2);

                alpha(j,k)=aflow*(diam0^4)/dist(j,k);
                alpha(k,j)=alpha(j,k);

                end
        end
end

nvess=0;
for j=1:tot
  for k=1:tot
	  if conn(j,k)>0
           nvess=nvess+1;
          end
  end
end

nvess=nvess/2;

for j=1:nodes
	for k=1:nodes
            a(j,k)=-alpha(j,k);
        end

        for k=1:tot
	    a(j,j)=a(j,j)+alpha(j,k);
        end
end

for j=1:nodes
	for k=1:2
		b(j)=b(j)+alpha(j,k+nodes)*pio(k);
        end 
end

%npress=inv(a)*b;
npress=a\b;
  
for j=1:nodes
	press(j)=npress(j);
    end

for j=1:2
	press(nodes+j)=pio(j);
end

for j=1:tot
	for k=1:tot
		if (conn(j,k) > 0) 
	            q(j,k)=(press(j)-press(k))*alpha(j,k);
                end
        end
end

n=1;

vess=zeros(nvess,3);

for j = 1:tot
  for k = 1:tot
    if q(j,k)>0
      vess(n,1)=j;
      vess(n,2)=k;
      vess(n,3)=q(j,k);
      n=n+1;
    end
  end
end

for j=1:tot
        nodepos(j,1)=poss(j,1);
        nodepos(j,2)=0;
        nodepos(j,3)=poss(j,2);
end

%csvwrite('nodepos.csv',nodepos)
%csvwrite('vess.csv',vess)








