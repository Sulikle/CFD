function UHLSr=HLSr(U,delta_x,delta_y,Coord_center,ESUEL)
global Nelement;
%% Preproceeding
Urect=zeros(2,Nelement);


%% Proceeding
%------LSR
for ie=1:Nelement
    LHS=zeros(2,2);
    RHS=zeros(2,1);
    SUEL=0;
    for k=1:3
        je=ESUEL(k,ie);
        if je<=Nelement
            SUEL=SUEL+1;
            B2=(Coord_center(1,je)-Coord_center(1,ie))/delta_x(1,ie);
            B3=(Coord_center(2,je)-Coord_center(2,ie))/delta_y(1,ie);
            LHS(1,1)=LHS(1,1)+B2^2;LHS(1,2)=LHS(1,2)+B2*B3;
            LHS(2,1)=LHS(2,1)+B2*B3;LHS(2,2)=LHS(2,2)+B3^2;
            RHS(1,1)=RHS(1,1)+B2*(U(1,je)-U(1,ie));
            RHS(2,1)=RHS(2,1)+B3*(U(1,je)-U(1,ie));
        end
    end
    if SUEL==1%when boundary element isn't sufficient
%        REPLACE=0;
        for k=1:3
            je=ESUEL(k,ie);
            if je<=Nelement
                for i=1:3
                    ke=ESUEL(i,je);%search for neiber's nerber
                    if ke~=ie&&ke<=Nelement
                       B2=(Coord_center(1,ke)-Coord_center(1,ie))/delta_x(1,ie);
                       B3=(Coord_center(2,ke)-Coord_center(2,ie))/delta_y(1,ie);
                       LHS(1,1)=LHS(1,1)+B2^2;LHS(1,2)=LHS(1,2)+B2*B3;
                       LHS(2,1)=LHS(2,1)+B2*B3;LHS(2,2)=LHS(2,2)+B3^2;
                       RHS(1,1)=RHS(1,1)+B2*(U(1,ke)-U(1,ie));
                       RHS(2,1)=RHS(2,1)+B3*(U(1,ke)-U(1,ie));
                       REPLACE=REPLACE+1;
                    end
%                     if REPLACE==2
%                        break;
%                     end
                end
            end
        end
    end
    Urect(:,ie)=lsqminnorm(LHS,RHS);
    %Urect(:,ie)=LHS\RHS;
end
UHLSr=Urect;


%------LSr
%In P0P1 LSr is equal to LSR

%------HLSr
%In P0P1 HLSr is equal to LSR too

end