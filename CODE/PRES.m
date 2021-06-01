function [p u v]=PRES(u,v,p,nx,ny,dx,dy,dt,uin,iter)

% count=1;
% N=zeros(ny,nx);
% for i=1:nx
%     for j=1:ny
%         N(j,i)=count;
%         count=count+1;
%     end
% end

n=nx*ny;
cx=1/(2*dx);
cy=1/(2*dy);
D=zeros(ny,nx);
for j=2:ny-1 %rows
    for k=2:nx-1 %cols
        D(j,k)=1/dt*((u(j,k)-u(j,k-1))/dx+(v(j,k)-v(j-1,k))/dy);
    end
end
%%
rhs(1:n)=0;
A=matrix(ny,nx,dx,dy);
for j=1:ny
    for k=1:nx
        l=j+(k-1)*ny;
        rhs(l)=D(j,k);        
        rhs(l)=rhs(l)-A(l,n);               
    end
end
rhs(l)=1;
A(:,l)=0;
A(l,:)=0;
A(l,l)=rhs(l);
jpvt(1:n)=0;
ww=1;
if(mod(iter,ww)==0 || iter==1)
    [AJ,jpvt]=fact(n,A,jpvt);
end
[rd]=solve(n,AJ,jpvt,rhs);
%%
zz=0;
for i=1:nx
    for j=1:ny
        zz=zz+1;
        p(j,i)=rd(zz);
    end
end
%%
for j=2:ny-1
    for k=2:nx
        if (k==nx)
            p(:,nx+1)=p(:,nx-1);
        end
        u(j,k)=u(j,k)-dt*cx*(p(j,k+1)-p(j,k-1));
    end;
end;
for j=2:ny-1
    for k=2:nx-1
        v(j,k)=v(j,k)-dt*cy*(p(j+1,k)-p(j-1,k));
    end
end
u(2:ny-1,1)=uin;

end
