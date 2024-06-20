function Vr=JudgeBC(ier,iel,Ma_n_l,v_n_l,n_vector,t_vector,Vcurrent,Ucurrent)
global Nelement;
global gamma;
global V_far;

global INPOEL;
global COORD;
global BCOND;
global ESUEL;
global INTFAC;
global ESUP1;
global ESUP2;
global Striangle;

        if ier>Nelement&&BCOND(3,ier-Nelement)==2%judge whether face is fixed
            rho_g=Vcurrent(1,iel);
            V_g=[Vcurrent(2,iel),Vcurrent(3,iel)]-2*v_n_l*n_vector;
            u_g=V_g(1,1);
            v_g=V_g(1,2);
            p_g=(gamma-1)*(Ucurrent(4,iel)-0.5*rho_g*(u_g^2+v_g^2));
            a_g=sqrt(gamma*p_g/rho_g);
        elseif ier>Nelement&&BCOND(3,ier-Nelement)==4&&Ma_n_l<=-1%judge whether face is far& judge size of M_n
                rho_g=V_far(1,1);
                u_g=V_far(2,1);
                v_g=V_far(3,1);
                p_g=V_far(4,1);
                a_g=V_far(5,1);
        elseif ier>Nelement&&BCOND(3,ier-Nelement)==4&&Ma_n_l>-1&&Ma_n_l<=0
                v_t=dot(V_far(2:3,1)',t_vector/norm(t_vector));
                v_n=0.5*(dot(V_far(2:3,1)',n_vector)+dot(Vcurrent(2:3,iel)',n_vector)+2*Vcurrent(5,iel)/(gamma-1)-2*V_far(5,1)/(gamma-1));
                v=[n_vector;t_vector/norm(t_vector)]\[v_n;v_t];
                u_g=v(1,1);v_g=v(2,1);
                a_g=(gamma-1)*(dot(Vcurrent(2:3,iel)',n_vector)-dot(V_far(2:3,1)',n_vector)+2*(V_far(5,1)+Vcurrent(5,iel))/(gamma-1))/4;
                rho_g=(a_g^2/(gamma*V_far(4,1)))^(1/(gamma-1));
                p_g=V_far(4,1)*rho_g^gamma;

        elseif ier>Nelement&&BCOND(3,ier-Nelement)==4&&Ma_n_l>0&&Ma_n_l<=1
                v_t=dot([Vcurrent(2,iel),Vcurrent(3,iel)],t_vector/norm(t_vector));
                v_n=0.5*(dot(V_far(2:3,1)',n_vector)+dot([Vcurrent(2,iel),Vcurrent(3,iel)],n_vector)+2*Vcurrent(5,iel)/(gamma-1)-2*V_far(5,1)/(gamma-1));
                v=[n_vector;t_vector/norm(t_vector)]\[v_n;v_t];
                u_g=v(1,1);v_g=v(2,1);
                a_g=(gamma-1)*(dot([Vcurrent(2,iel),Vcurrent(3,iel)],n_vector)-dot(V_far(2:3,1)',n_vector)+2*(V_far(5,1)+Vcurrent(5,iel))/(gamma-1))/4;
                rho_g=(a_g^2*Vcurrent(1,iel)^gamma/(gamma*Vcurrent(4,iel)))^(1/(gamma-1));
                p_g=a_g^2*rho_g/gamma;

        
        elseif ier>Nelement&&BCOND(3,ier-Nelement)==4&&Ma_n_l>1
               rho_g=Vcurrent(1,iel);
               u_g=Vcurrent(2,iel);
               v_g=Vcurrent(3,iel);
               p_g=Vcurrent(4,iel);
               a_g=Vcurrent(5,iel);

        else 
               rho_g=Vcurrent(1,ier);
               u_g=Vcurrent(2,ier);
               v_g=Vcurrent(3,ier);
               p_g=Vcurrent(4,ier);
               a_g=Vcurrent(5,ier);            
        end

        Vr=[rho_g;u_g;v_g;p_g;a_g];


end