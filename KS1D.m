clear; clc; close all;

nx=100; Lx=-1; Rx=1; h=(Rx-Lx)/nx; x=(Lx+0.5*h:h:Rx-0.5*h)'; dt=0.1*h;
max_it=10; ns=max_it/10;

L=-2*eye(nx);
for ix=1:nx-1
    L(ix,ix+1)=1; L(ix+1,ix)=1;
end
L(1,1)=L(1,1)+L(1,2); L(nx,nx)=L(nx,nx)+L(nx,nx-1); L=L/h^2;

u=2*(2.1+cos(pi*x)); v=u;
mass(1)=sum(u,'all')*h;
ene(1)=sum(0.5*v.^2-0.5*v.*(L*v)-u.*v+u.*log(u)-u,'all')*h;

plot(x,u, x,v);
drawnow;

for it=1:max_it
    maxGS=1000; itGS=0; tol=1e-10; ou=u;
    oup=0.5*(ou([2:end,end])+ou); oun=0.5*(ou([1,1:end-1])+ou);
    while true
        itGS=itGS+1;

        B=zeros(nx);
        for ix=1:nx-1
            B(ix,ix)=B(ix,ix)-oup(ix)/u(ix);
            B(ix,ix+1)=oup(ix)/u(ix+1);
        end
        for ix=2:nx
            B(ix,ix)=B(ix,ix)-oun(ix)/u(ix);
            B(ix,ix-1)=oun(ix)/u(ix-1);
        end
        B=B/h^2;

        cof=eye(nx)-dt*B;
        sor=ou+(dt/h^2)*(oup.*(log(u([2:end,end]))-v([2:end,end])) ...
            -(oup+oun).*(log(u)-v)+oun.*(log(u([1,1:end-1]))-v([1,1:end-1])));
        nu=cof\sor;

        err=norm(nu-u,'inf');
        u=nu;
        if (err<tol)||(itGS>maxGS)
            % fprintf('it=%d, err=%1.3e\n',itGS,err)
            break;
        end
    end

    cof=(1+dt)*eye(nx)-dt*L; sor=v+dt*u;
    v=cof\sor;

    mass(it+1)=sum(u,'all')*h;
    ene(it+1)=sum(0.5*v.^2-0.5*v.*(L*v)-u.*v+u.*log(u)-u,'all')*h;

    if mod(it,ns)==0
        clf;
        subplot(1,2,1);
        plot(x,u, x,v);
        grid on;

        subplot(1,2,2);
        plot((0:it)/max_it,mass/mass(1));
        hold on;
        plot((0:it)/max_it,ene/ene(1));
        grid on;
        drawnow;
    end
end
