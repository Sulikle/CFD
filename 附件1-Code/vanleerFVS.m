function F_ij=vanleerFVS(Ma_n_l,Ma_n_r,v_n_l,v_n_r,Vl,Ul,Vr,Ur,n_vector)
        
global gamma;
%solve fi+
        if (Ma_n_l>1)
            F_l_plus=[Vl(1,1)*v_n_l;v_n_l*Vl(1,1)*Vl(2,1)+Vl(4,1)*n_vector(1,1);v_n_l*Vl(1,1)*Vl(3,1)+Vl(4,1)*n_vector(1,2);v_n_l*(Ul(4,1)+Vl(4,1))];
        elseif(Ma_n_l<-1)
            F_l_plus=0;
        else 
            f_l_plus=Vl(1,1)*Vl(5,1)*(Ma_n_l+1)^2/4;
            F_l_plus(1,1)=f_l_plus*1;
            F_l_plus(2,1)=f_l_plus*(Vl(2,1)+(-v_n_l+2*Vl(5,1))*n_vector(1,1)/gamma);
            F_l_plus(3,1)=f_l_plus*(Vl(3,1)+(-v_n_l+2*Vl(5,1))*n_vector(1,2)/gamma);
            F_l_plus(4,1)=f_l_plus*((Vl(2,1)^2+Vl(3,1)^2-v_n_l^2)/2+((gamma-1)*v_n_l+2*Vl(5,1))^2/(2*(gamma^2-1)));
        end

%solve fj-
        if(Ma_n_r>1)
           F_r_minus=0;
        elseif(Ma_n_r<-1)
            F_r_minus=[Vr(1,1)*v_n_r;v_n_r*Vr(1,1)*Vr(2,1)+Vr(4,1)*n_vector(1,1);v_n_r*Vr(1,1)*Vr(3,1)+Vr(4,1)*n_vector(1,2);v_n_r*(Ur(4,1)+Vr(4,1))];
        else
            f_r_minus=-Vr(1,1)*Vr(5,1)*(Ma_n_r-1)^2/4;
            F_r_minus(1,1)=f_r_minus;
            F_r_minus(2,1)=f_r_minus*(Vr(2,1)+(-v_n_r-2*Vr(5,1))*n_vector(1,1)/gamma);
            F_r_minus(3,1)=f_r_minus*(Vr(3,1)+(-v_n_r-2*Vr(5,1))*n_vector(1,2)/gamma);
            F_r_minus(4,1)=f_r_minus*((Vr(2,1)^2+Vr(3,1)^2-v_n_r^2)/2+((gamma-1)*v_n_r-2*Vr(5,1))^2/(2*(gamma^2-1)));
        end
              

        F_ij=F_l_plus+F_r_minus;%flux

end