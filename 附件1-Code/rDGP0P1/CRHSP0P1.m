function RHS=CRHSP0P1(Ucurrent,Urect,delta_x,delta_y,Coord_center)
global Flux;
global Nface;
global Nelement;
global gamma;
global Mesh;
% Geometry
global INPOEL;
global COORD;
global BCOND;
global ESUEL;
global INTFAC;
global ESUP1;
global ESUP2;
global Striangle;

RHS=zeros(4,Nelement);
tf=[-sqrt(15)/5,0,sqrt(15)/5];Wf=[5/9,8/9,5/9];%-Gauss in face

%calculate RHS_boundary
for iface=1:Nface
        iel=INTFAC(1,iface);ier=INTFAC(2,iface);
        xcl=Coord_center(1,iel);ycl=Coord_center(2,iel);
        if ier<=Nelement
           xcr=Coord_center(1,ier);ycr=Coord_center(2,ier);
        end
        %calculate V_n
        ip1=INTFAC(3,iface);%points in face
        ip2=INTFAC(4,iface);
        x1=COORD(1,ip1);x2=COORD(1,ip2);
        y1=COORD(2,ip1);y2=COORD(2,ip2);
        Coord_face=[x1,x2;y1,y2];


        n_vector=INTFAC(6:7,iface)';
        t_vector=[COORD(1,ip2)-COORD(1,ip1),COORD(2,ip2)-COORD(2,ip1)];
        tau_ij=INTFAC(8,iface);%face length
        xi=0.5*(x1+x2);yi=0.5*(y1+y2);
        if Mesh==1.1||Mesh==1.2||Mesh==1.3||Mesh==1.4
           if ier>Nelement&&BCOND(3,ier-Nelement)==2
                 n_vector=[0-xi,0-yi]/norm([0-xi,0-yi]);
                 t_vector=[-yi,xi]/norm([-yi,xi]);
           end
        end

            B2l=(xi-xcl)/delta_x(iel);B3l=(yi-ycl)/delta_y(iel);
            if ier<=Nelement
               phir=BJlimiter(ier,Ucurrent,Urect,Coord_face,Coord_center,ESUEL);
               B2r=(xi-xcr)/delta_x(ier);B3r=(yi-ycr)/delta_y(ier);
               rho_ir=Ucurrent(1,ier)+phir(1,1)*(Urect(1,1,ier)*B2r+Urect(2,1,ier)*B3r);
               if rho_ir<0
                   rho_ir=1e-9;
               end
               rhou_ir=Ucurrent(2,ier)+phir(2,1)*(Urect(1,2,ier)*B2r+Urect(2,2,ier)*B3r);
               rhov_ir=Ucurrent(3,ier)+phir(3,1)*(Urect(1,3,ier)*B2r+Urect(2,3,ier)*B3r);
               rhoe_ir=Ucurrent(4,ier)+phir(4,1)*(Urect(1,4,ier)*B2r+Urect(2,4,ier)*B3r);
               if rhoe_ir<0
                   rhoe_ir=1e-9;
               end               
               p_ir=(gamma-1)*(rhoe_ir-0.5*rho_ir*((rhou_ir/rho_ir)^2+(rhov_ir/rho_ir)^2));
               if p_ir<0
                   p_ir=1e-9;
               end
               a_ir=sqrt(gamma*p_ir/rho_ir);
               V_ir=[rho_ir;rhou_ir/rho_ir;rhov_ir/rho_ir;p_ir;a_ir];
            else
                V_ir=zeros(5,1);
            end

            phil=BJlimiter(iel,Ucurrent,Urect,Coord_face,Coord_center,ESUEL);
            rho_il=Ucurrent(1,iel)+phil(1,1)*(Urect(1,1,iel)*B2l+Urect(2,1,iel)*B3l);
            rhou_il=Ucurrent(2,iel)+phil(2,1)*(Urect(1,2,iel)*B2l+Urect(2,2,iel)*B3l);
            rhov_il=Ucurrent(3,iel)+phil(3,1)*(Urect(1,3,iel)*B2l+Urect(2,3,iel)*B3l);
            rhoe_il=Ucurrent(4,iel)+phil(4,1)*(Urect(1,4,iel)*B2l+Urect(2,4,iel)*B3l);
            if rho_il<0
                rho_il=1e-9;
            end
            if rhoe_il<0
                   rhoe_il=1e-9;
            end
            p_il=(gamma-1)*(rhoe_il-0.5*rho_il*((rhou_il/rho_il)^2+(rhov_il/rho_il)^2));
            if p_il<0
                p_il=1e-9;
            end
            a_il=sqrt(gamma*p_il/rho_il);

            vl=[rhou_il/rho_il,rhov_il/rho_il];
            v_n_l=dot(vl,n_vector);
            Ma_n_l=v_n_l/a_il;
            V_il=[rho_il;rhou_il/rho_il;rhov_il/rho_il;p_il;a_il];

            U_il=[rho_il;rhou_il;rhov_il;rhoe_il];
            
            Vr=JudgeBC1(ier,iel,Ma_n_l,v_n_l,n_vector,t_vector,V_il,V_ir,U_il);
            
            v_n_r=dot(Vr(2:3,1)',n_vector);

            rho_g=Vr(1,1);u_g=Vr(2,1);v_g=Vr(3,1);p_g=Vr(4,1);a_g=Vr(5,1);
            if rho_g<0
                rho_g=1e-9;
            end
            if p_g<0
                p_g=1e-9;
            end
            
            Ma_n_r=v_n_r/a_g;

            Ur=[rho_g;rho_g*u_g;rho_g*v_g;p_g/(gamma-1)+0.5*rho_g*(u_g^2+v_g^2)];



%-flux        
        if Flux==1
           F_ij=vanleerFVS(Ma_n_l,Ma_n_r,v_n_l,v_n_r,V_il,U_il,Vr,Ur,n_vector);
        elseif Flux==2
           F_ij=Roe(v_n_l,v_n_r,V_il,U_il,Vr,Ur,n_vector);
        end
%-rhs

         RHS(:,iel)=RHS(:,iel)-F_ij*tau_ij;
        if ier<=Nelement
            RHS(:,ier)=RHS(:,ier)+F_ij*tau_ij;
        end

end



end

