function UGG=GG(U,delta_x,delta_y,Coord_center,ESUEL)
global Nelement;
%% Preproceeding
global Nface;
global INTFAC;
global ESUP1;
global ESUP2;
global Striangle;
Urect=zeros(2,Nelement);
RHS=zeros(2,Nelement);


%% Proceeding
%------GG
for iface=1:Nface
    iel=INTFAC(1,iface);
    ier=INTFAC(2,iface);
    ip1=INTFAC(3,iface);
    ip2=INTFAC(4,iface);
    n_x=INTFAC(6,iface);
    n_y=INTFAC(7,iface);
    tau_ij=INTFAC(8,iface);
    u_a=0;u_b=0;weight=zeros(2,1);
    %calculate u(a)&u(b)
    for istor=ESUP2(ip1)+1:ESUP2(ip1+1)%loop over elem surrounding ip1
        jelem=ESUP1(istor);
        u_a=u_a+Striangle(1,jelem)*U(1,jelem);
        weight(1,1)=weight(1,1)+Striangle(1,jelem);
    end
    u_a=u_a/weight(1,1);
    
    for istor=ESUP2(ip2)+1:ESUP2(ip2+1)%loop over elem surrounding ip2
        jelem=ESUP1(istor);
        u_b=u_b+Striangle(1,jelem)*U(1,jelem);
        weight(2,1)=weight(2,1)+Striangle(1,jelem);
    end
    u_b=u_b/weight(2,1);

    RHS(1,iel)=RHS(1,iel)+0.5*(u_a+u_b)*n_x*tau_ij;RHS(2,iel)=RHS(2,iel)+0.5*(u_a+u_b)*n_y*tau_ij;
    if ier<=Nelement
    RHS(1,ier)=RHS(1,ier)-0.5*(u_a+u_b)*n_x*tau_ij;RHS(2,ier)=RHS(2,ier)-0.5*(u_a+u_b)*n_y*tau_ij;
    end
end

Urect(1,:)=(RHS(1,:).*delta_x)./Striangle;%Ux*delta_x&Uy*delta_y
Urect(2,:)=(RHS(2,:).*delta_y)./Striangle;
UGG=Urect;


%------LSr
%In P0P1 LSr is equal to LSR

%------HLSr
%In P0P1 HLSr is equal to LSR too

end