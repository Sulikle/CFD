function deltat=Deltat(ESUEL,INPOEL,INTFAC,Vcurrent,Ucurrent,Striangle,COORD,BCOND)
%-Auxiliary array
LNOFA=zeros(2,3);
LNOFA(1,1)=2;LNOFA(2,1)=3;
LNOFA(1,2)=3;LNOFA(2,2)=1;
LNOFA(1,3)=1;LNOFA(2,3)=2;
%global paramaters
global Ma_far;
global afa_attack;
global gamma;
global CFL;
global Ndim;
global Nelement;
global Npoint;
global Nnode;
global NBface;
global Nface;
%calculate deltat
ADLTA=zeros(1,Nelement);
for ie=1:Nelement
    for IFAEL=1:Nnode
        if(ESUEL(IFAEL,ie)>ie)
            iel=ie;ier=ESUEL(IFAEL,ie);
        else
            iel=ESUEL(IFAEL,ie);ier=ie;
        end

        ip1=INPOEL(LNOFA(1,IFAEL),ie);%points in face
        ip2=INPOEL(LNOFA(2,IFAEL),ie);
        tau_ij=sqrt((COORD(1,ip1)-COORD(1,ip2))^2+(COORD(2,ip1)-COORD(2,ip2))^2);
        for iface=1:Nface
            if INTFAC(1,iface)==iel&&INTFAC(2,iface)==ier
                n_vector=[INTFAC(6,iface),INTFAC(7,iface)];
            end
        end
        t_vector=[COORD(1,ip2)-COORD(1,ip1),COORD(2,ip2)-COORD(2,ip1)]/norm([COORD(1,ip2)-COORD(1,ip1),COORD(2,ip2)-COORD(2,ip1)]);
        v_n_l=dot(Vcurrent(2:3,iel),n_vector);
        Ma_n_l=v_n_l/Vcurrent(5,iel);

         Vr=JudgeBC(ier,iel,Ma_n_l,v_n_l,n_vector,t_vector,Vcurrent,Ucurrent);
         rho_g=Vr(1,1);u_g=Vr(2,1);v_g=Vr(3,1);p_g=Vr(4,1);a_g=Vr(5,1);
        ADLTA(1,ie)=ADLTA(1,ie)+tau_ij*0.5*(abs(n_vector(1,1)*(Vcurrent(2,iel)+u_g)+n_vector(1,2)*(Vcurrent(3,iel)+v_g))+Vcurrent(5,iel)+a_g);
    end
end
deltat1=CFL*Striangle./ADLTA;
deltat=min(deltat1);
end