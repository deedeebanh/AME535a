function A=matrix(ny,nx,dx,dy)

ccx=1/(dx*dx);
ccy=1/(dy*dy);
n=nx*ny;

A(1:n,1:n)=0;
ib=ny;
it=ny+1;
for j=1:ny
    for k=1:nx
        l=j+(k-1)*ny;
        A(l,l)=-2*(ccx+ccy);
        
        xr=l+ny;
        if(xr<=n)
            A(l,xr)=ccx;
            if(l<=ny)
                A(l,xr)=2*ccx;
            end
        end
        
        xl=l-ny;
        if(xl>0)
            A(l,xl)=ccx;
            if(l>=n-ny+1)
                A(l,xl)=2*ccx;
            end
        end
        
        yr=l+1;
        if(yr<=n)            
            A(l,yr)=ccy;
            if(yr==it && it<n)
                A(yr,it+1)=2*ccy; 
                A(yr,it-1)=0;
                it=it+ny;  
            end
        end      
        
        yl=l-1;
        if(yl>0)
            A(l,yl)=ccy;
            if(l==ib && ib<n)
                A(l,ib+1)=0;
                A(l,ib-1)=2*ccy;
                ib=ib+ny;
            end                
        end  
    end; %k        
end; %j
A(1,2)=2*ccy;
A(l,l-1)=2*ccy;