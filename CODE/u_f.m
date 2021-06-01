function [u_f] = u_f(u,v,nx,ny,dx,dy,re)
% Calculate the Nonlinear step
% for Dudt=Fu 
u_f=u;
cx=1/2/dx;
cy=1/2/dy;
for j=2:ny-1 %rows
    for k=2:nx %cols        
        if (k==nx)
            u(j,k+1)=u(j,k-1);
        end
        a=u(j,k)*(u(j,k+1)-u(j,k-1))*cx;
        b=v(j,k)*(u(j+1,k)-u(j-1,k))*cy;
        c=3/re;
        u_f(j,k)=-a-b+c;
    end
end
