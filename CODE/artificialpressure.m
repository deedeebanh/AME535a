function rms=artificialpressure(nx,ny)

dx=1/(nx-1);
dy=1/(ny-1);

n=nx*ny;

count=1;
N=zeros(ny,nx);
for i=1:nx
    for j=1:ny
        N(j,i)=count;
        count=count+1;
    end
end

x=0:dx:1;
y=0:dy:1;
A=matrix(ny,nx,dx,dy);
for j=1:ny
    for k=1:nx  
       pex(j,k)=cos(pi*x(k))*cos(pi*y(j)); 
       Dex(j,k)=-2*pi^2*cos(pi*x(k))*cos(pi*y(j));
    end;
end;
pex=reshape(pex,[],1);
for j=1:ny
    for k=1:nx
        l=j+(k-1)*ny;
        rhs(l)=Dex(j,k);        
        rhs(l)=rhs(l)-A(l,n);               
    end
end
rhs(l)=1;
A(:,l)=0;
A(l,:)=0;
A(l,l)=rhs(l);

jpvt(1:n)=0;
[AJ,jpvt]=fact(n,A,jpvt);
[rd]=solve(n,AJ,jpvt,rhs);
zz=0;
for i=1:nx
    for j=1:ny
        zz=zz+1;
        p(j,i)=rd(zz);
        moddex(j,i)=rhs(zz);
    end
end
sum=0;
for j=1:n
    sum = sum+ (pex(j)-p(j))^2;    
    err(j)=pex(j)-p(j);
end
rms=sqrt(sum/n);
% x=1:n;
% plot(x,err)
