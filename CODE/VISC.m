function u=VISC(u,nx,ny,dx,dy,re,dt,uin,ch)

r=zeros(ny,nx);
b=zeros(5,nx);

%% In x-direction
ccxa=dt/(2*re*dx*dx);
ccya=dt/(2*re*dy*dy);
for j=2:nx %cols
    for k=2:ny-1 %rows
        if (ch==1)
            u(k,nx+1)=u(k,nx-1);
        else
            u(k,nx+1)=0;
        end            
        r(k,j)=u(k,j)+ccxa*(u(k,j-1)-2*u(k,j)+u(k,j+1));
    end
end

for k=2:ny %rows
    for j=2:nx %cols
        jm=j-1;
        b(2,jm)=-ccxa;
        b(3,jm)=1+(2*ccxa);
        b(4,jm)=-ccxa;
        rrt(jm)=r(k,j);
    end
    if(ch==1)
        u(k,nx+1)=u(k,nx-1);
        b(2,jm)=b(2,jm)+b(4,jm);
    else        
        rrt(jm)=rrt(jm)-b(4,jm)*u(k,nx);
    end;  
    rrt(1) = rrt(1) - b(2,1)*u(k,1);
    b(4,jm)=0;
    b(2,1)=0;
    % solve tridiagonal system of eqns
    b=banfac(b,nx-1);
    ddt=bansol(rrt,b,nx-1);
    for j=2:nx
        jm=j-1;
        u(k,j)=ddt(jm);
    end;
end

%% tridiagonal systems in the y-direction
for j=2:nx+1 %cols
    for k=2:ny-1 %rows
        if (ch==1)
            u(k,nx+1)=u(k,nx-1);
        else
            u(k,nx+1)=0;
        end
        r(k,j)=u(k,j)+ccya*(u(k-1,j)-2*u(k,j)+u(k+1,j));
    end
end
for j=2:nx+1 %cols
    for k=2:ny-1 %rows
        km=k-1;
        b(2,km)=-ccya;
        b(3,km)=1+(2*ccya);
        b(4,km)=-ccya;
        rrt(km)=r(k,j);
    end;
    rrt(km)=rrt(km)-b(4,km)*u(ny,j);    
    % ensure RHS gets contribution from BC
    rrt(1) = rrt(1) - b(2,1)*u(1,j);
    b(4,km)=0;
    b(2,1)=0;    
    % solve tridiagonal system of eqns
    b=banfac(b,ny-2);
    ddt=bansol(rrt,b,ny-2);
    % increment t
    for k=2:ny-1
        km=k-1;
        u(k,j)=ddt(km);
    end;
end; %j

%%
u=u(1:ny,1:nx);
if (ch==1)
    u(2:ny-1,1)=uin;
else
    u(1,:)=0;
    u(:,nx)=0;
    u(:,1)=0;
    u(ny,:)=0;
end