function RHS=CRHSP1(Ucurrent,delta_x,delta_y,Coord_center)
global Flux;
global Nface;
global Nelement;
global gamma;
global Mesh;
% Geometry
global INPOEL;
global COORD;
global BCOND;
global INTFAC;
global Striangle;
RHS=zeros(3,4,Nelement);
ts=[1/2,1/2,0;0,1/2,1/2;1/2,0,1/2];Ws=[1/3,1/3,1/3];%-Gauss in square
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

        n_vector=INTFAC(6:7,iface)';
        t_vector=[COORD(1,ip2)-COORD(1,ip1),COORD(2,ip2)-COORD(2,ip1)];
        tau_ij=INTFAC(8,iface);%face length
        Il=zeros(12,1);Ir=zeros(8,1);
        for i=1:3%Gauss
            xi=x1*0.5*(1-tf(i))+x2*0.5*(1+tf(i));
            yi=y1*0.5*(1-tf(i))+y2*0.5*(1+tf(i));

            if Mesh==1.1||Mesh==1.2||Mesh==1.3||Mesh==1.4
                 if ier>Nelement&&BCOND(3,ier-Nelement)==2
                    n_vector=[0-xi,0-yi]/norm([0-xi,0-yi]);
                    t_vector=[-yi,xi]/norm([-yi,xi]);
                 end
            end
            B2l=(xi-xcl)/delta_x(iel);B3l=(yi-ycl)/delta_y(iel);
            if ier<=Nelement
               B2r=(xi-xcr)/delta_x(ier);B3r=(yi-ycr)/delta_y(ier);
               rho_ir=Ucurrent(1,1,ier)+Ucurrent(2,1,ier)*B2r+Ucurrent(3,1,ier)*B3r;
               if rho_ir<0
                   rho_ir=1e-9;
               end
               rhou_ir=Ucurrent(1,2,ier)+Ucurrent(2,2,ier)*B2r+Ucurrent(3,2,ier)*B3r;
               rhov_ir=Ucurrent(1,3,ier)+Ucurrent(2,3,ier)*B2r+Ucurrent(3,3,ier)*B3r;
               rhoe_ir=Ucurrent(1,4,ier)+Ucurrent(2,4,ier)*B2r+Ucurrent(3,4,ier)*B3r;              
               p_ir=(gamma-1)*(rhoe_ir-0.5*rho_ir*((rhou_ir/rho_ir)^2+(rhov_ir/rho_ir)^2));
               if p_ir<0
                   p_ir=1e-9;
               end               
               a_ir=sqrt(gamma*p_ir/rho_ir);
               V_ir=[rho_ir;rhou_ir/rho_ir;rhov_ir/rho_ir;p_ir;a_ir];
            else
                V_ir=zeros(5,1);
            end

            rho_il=Ucurrent(1,1,iel)+Ucurrent(2,1,iel)*B2l+Ucurrent(3,1,iel)*B3l;
            if rho_il<0
                   rho_il=1e-9;
            end            
            rhou_il=Ucurrent(1,2,iel)+Ucurrent(2,2,iel)*B2l+Ucurrent(3,2,iel)*B3l;
            rhov_il=Ucurrent(1,3,iel)+Ucurrent(2,3,iel)*B2l+Ucurrent(3,3,iel)*B3l;
            rhoe_il=Ucurrent(1,4,iel)+Ucurrent(2,4,iel)*B2l+Ucurrent(3,4,iel)*B3l;
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

        Il(1:4,1)=Il(1:4,1)+Wf(i)*F_ij;
        Il(5:8,1)=Il(5:8,1)+Wf(i)*(F_ij*(xi-xcl)/delta_x(1,iel));
        Il(9:12,1)=Il(9:12,1)+Wf(i)*(F_ij*(yi-ycl)/delta_y(1,iel));
        if ier<=Nelement
        Ir(1:4,1)=Ir(1:4,1)+Wf(i)*(F_ij*(xi-xcr)/delta_x(1,ier));
        Ir(5:8,1)=Ir(5:8,1)+Wf(i)*(F_ij*(yi-ycr)/delta_y(1,ier));
        end

        end
        Il=Il*0.5*tau_ij;
        Ir=Ir*0.5*tau_ij;
%-rhs

        RHS(1,1,iel)=RHS(1,1,iel)-Il(1,1);
        RHS(1,2,iel)=RHS(1,2,iel)-Il(2,1);
        RHS(1,3,iel)=RHS(1,3,iel)-Il(3,1);
        RHS(1,4,iel)=RHS(1,4,iel)-Il(4,1);

        RHS(2,1,iel)=RHS(2,1,iel)-Il(5,1);
        RHS(2,2,iel)=RHS(2,2,iel)-Il(6,1);
        RHS(2,3,iel)=RHS(2,3,iel)-Il(7,1);
        RHS(2,4,iel)=RHS(2,4,iel)-Il(8,1);

        RHS(3,1,iel)=RHS(3,1,iel)-Il(9,1);
        RHS(3,2,iel)=RHS(3,2,iel)-Il(10,1);
        RHS(3,3,iel)=RHS(3,3,iel)-Il(11,1);
        RHS(3,4,iel)=RHS(3,4,iel)-Il(12,1);

        if ier<=Nelement
            RHS(1,1,ier)=RHS(1,1,ier)+Il(1,1);
            RHS(1,2,ier)=RHS(1,2,ier)+Il(2,1);
            RHS(1,3,ier)=RHS(1,3,ier)+Il(3,1);
            RHS(1,4,ier)=RHS(1,4,ier)+Il(4,1);

            RHS(2,1,ier)=RHS(2,1,ier)+Ir(1,1);
            RHS(2,2,ier)=RHS(2,2,ier)+Ir(2,1);
            RHS(2,3,ier)=RHS(2,3,ier)+Ir(3,1);
            RHS(2,4,ier)=RHS(2,4,ier)+Ir(4,1);

            RHS(3,1,ier)=RHS(3,1,ier)+Ir(5,1);
            RHS(3,2,ier)=RHS(3,2,ier)+Ir(6,1);
            RHS(3,3,ier)=RHS(3,3,ier)+Ir(7,1);
            RHS(3,4,ier)=RHS(3,4,ier)+Ir(8,1);
        end
end


%calculate RHS_domain

for ie=1:Nelement
    x_1=COORD(1,INPOEL(1,ie));
    x_2=COORD(1,INPOEL(2,ie));
    x_3=COORD(1,INPOEL(3,ie));
    y_1=COORD(2,INPOEL(1,ie));
    y_2=COORD(2,INPOEL(2,ie));
    y_3=COORD(2,INPOEL(3,ie));

    x_trale=[x_1,x_2,x_3];
    y_trale=[y_1,y_2,y_3];

    Square=Striangle(1,ie);xc=Coord_center(1,ie);yc=Coord_center(2,ie);

    RHS(1,:,ie)=RHS(1,:,ie)+0;
    I=zeros(1,8);fi=zeros(1,8);
   for i=1:3
       xi=ts(i,1)*x_trale(1)+ts(i,2)*x_trale(2)+ts(i,3)*x_trale(3);
       yi=ts(i,1)*y_trale(1)+ts(i,2)*y_trale(2)+ts(i,3)*y_trale(3);
       B2=(xi-xc)/delta_x(ie);B3=(yi-yc)/delta_y(ie);
       rho_i=Ucurrent(1,1,ie)+Ucurrent(2,1,ie)*B2+Ucurrent(3,1,ie)*B3;
               if rho_i<0
                   rho_i=1e-9;
               end
       rhou_i=Ucurrent(1,2,ie)+Ucurrent(2,2,ie)*B2+Ucurrent(3,2,ie)*B3;
       rhov_i=Ucurrent(1,3,ie)+Ucurrent(2,3,ie)*B2+Ucurrent(3,3,ie)*B3;
       rhoe_i=Ucurrent(1,4,ie)+Ucurrent(2,4,ie)*B2+Ucurrent(3,4,ie)*B3;
       p_i=(gamma-1)*(rhoe_i-0.5*rho_i*((rhou_i/rho_i)^2+(rhov_i/rho_i)^2));
                if p_i<0
                   p_i=1e-9;
                end      
       fi(1,1)=rhou_i;%rho_e
       fi(1,2)=rhov_i;%rho_e
       fi(1,3)=p_i+rho_i*((rhou_i/rho_i)^2);%rho*u^2+p
       fi(1,4)=rho_i*(rhou_i/rho_i)*(rhov_i/rho_i);%rho*u*v
       fi(1,5)=fi(1,4);
       fi(1,6)=p_i+rho_i*((rhov_i/rho_i)^2);%rho*v^2+p
       fi(1,7)=(rhou_i/rho_i)*(rhoe_i+p_i);
       fi(1,8)=(rhov_i/rho_i)*(rhoe_i+p_i);

       I(1,:)=I(1,:)+Ws(i)*fi(1,:);
   end
   I(1,:)=Square*I(1,:);

    RHS(2,1,ie)=RHS(2,1,ie)+I(1,1)/delta_x(1,ie);
    RHS(3,1,ie)=RHS(3,1,ie)+I(1,2)/delta_y(1,ie);
    RHS(2,2,ie)=RHS(2,2,ie)+I(1,3)/delta_x(1,ie);
    RHS(3,2,ie)=RHS(3,2,ie)+I(1,4)/delta_y(1,ie);
    RHS(2,3,ie)=RHS(2,3,ie)+I(1,5)/delta_x(1,ie);
    RHS(3,3,ie)=RHS(3,3,ie)+I(1,6)/delta_y(1,ie);
    RHS(2,4,ie)=RHS(2,4,ie)+I(1,7)/delta_x(1,ie);
    RHS(3,4,ie)=RHS(3,4,ie)+I(1,8)/delta_y(1,ie);
    

end







end

