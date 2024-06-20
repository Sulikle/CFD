function RHS=CRHSP0(Vcurrent,Ucurrent)
global Flux;
global Nface;
global Nelement;
global gamma;

global INPOEL;
global COORD;
global BCOND;
global ESUEL;
global INTFAC;
global ESUP1;
global ESUP2;
global Striangle;

RHS=zeros(4,Nelement);

    for iface=1:Nface
        iel=INTFAC(1,iface);ier=INTFAC(2,iface);
        %calculate V_n
        ip1=INTFAC(3,iface);%points in face
        ip2=INTFAC(4,iface);
        n_vector=INTFAC(6:7,iface)';
        t_vector=[COORD(1,ip2)-COORD(1,ip1),COORD(2,ip2)-COORD(2,ip1)];
        tau_ij=INTFAC(8,iface);%face length

        v_n_l=dot(Vcurrent(2:3,iel)',n_vector);
        Ma_n_l=v_n_l/Vcurrent(5,iel);
        Vr=JudgeBC(ier,iel,Ma_n_l,v_n_l,n_vector,t_vector,Vcurrent,Ucurrent);
        v_n_r=dot(Vr(2:3,1)',n_vector);
        rho_g=Vr(1,1);u_g=Vr(2,1);v_g=Vr(3,1);p_g=Vr(4,1);a_g=Vr(5,1);
        Ma_n_r=v_n_r/a_g;
                
        Vl=Vcurrent(:,iel);
        Ul=Ucurrent(:,iel);
        Ur=[rho_g;rho_g*u_g;rho_g*v_g;p_g/(gamma-1)+0.5*rho_g*(u_g^2+v_g^2)];
%-flux        
        if Flux==1
           F_ij=vanleerFVS(Ma_n_l,Ma_n_r,v_n_l,v_n_r,Vl,Ul,Vr,Ur,n_vector);
        elseif Flux==2
           F_ij=Roe(v_n_l,v_n_r,Vl,Ul,Vr,Ur,n_vector);
        end
%-rhs
        RHS(:,iel)=RHS(:,iel)-F_ij*tau_ij;
        if ier<=Nelement
            RHS(:,ier)=RHS(:,ier)+F_ij*tau_ij;

        end
    end

end