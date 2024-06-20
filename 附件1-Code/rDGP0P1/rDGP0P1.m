function OK=rDGP0P1
%% Preproceeding
%global paramaters
global Timeadvance;
global Ma_far;
global V_far;
global gamma;
global Ndim;
global Nelement;
global Npoint;
global tol;
global rDG;

% Geometry
global INPOEL;
global COORD;
global BCOND;
global ESUEL;
global INTFAC;
global ESUP1;
global ESUP2;
global Striangle;

%calculate delta_x,delta_y,xc,yc
delta_x=zeros(1,Nelement);
delta_y=zeros(1,Nelement);
Coord_center=zeros(2,Nelement);

for ie=1:Nelement
x_1=COORD(1,INPOEL(1,ie));
x_2=COORD(1,INPOEL(2,ie));
x_3=COORD(1,INPOEL(3,ie));
y_1=COORD(2,INPOEL(1,ie));
y_2=COORD(2,INPOEL(2,ie));
y_3=COORD(2,INPOEL(3,ie));

x_trale=[x_1,x_2,x_3];
y_trale=[y_1,y_2,y_3];
%In this subroutine we calculate delta_x&delta_y
delta_x(1,ie)=0.5*(max(x_trale)-min(x_trale));
delta_y(1,ie)=0.5*(max(y_trale)-min(y_trale));
%In this subroutine we calculate x_c&y_c
Coord_center(1,ie)=sum(x_trale)/3;
Coord_center(2,ie)=sum(y_trale)/3;
end

t=[1/2,1/2,0;0,1/2,1/2;1/2,0,1/2];W=[1/3,1/3,1/3];%Gauss square

%rDGP0P1
Ucurrent=zeros(4,Nelement);%rho rhou rhov rhoe
Urect=zeros(2,4,Nelement);

%-TVDRK3
Uhold=zeros(4,Nelement);
Unhold=zeros(4,Nelement);
afa=[1,3/4,1/3];gamma1=[0,1/4,2/3];belta=[1,1/4,2/3];
Unext=zeros(4,Nelement);


%% Proceeding
%initial condition
%远场初值
for ie=1:Nelement
    Ucurrent(1,ie)=V_far(1,1);
    Ucurrent(2,ie)=V_far(2,1);
    Ucurrent(3,ie)=V_far(3,1);
    Ucurrent(4,ie)=0.5+1/(gamma*(gamma-1)*Ma_far^2);
end
Vcurrent=convert(Ucurrent);

deltat=Deltat(ESUEL,INPOEL,INTFAC,Vcurrent,Ucurrent,Striangle,COORD,BCOND);
%deltat=0.0005;
%time iteration

endtimes=1e5;
nowtime=0;
for itimes=1:endtimes
%-reconstruct
if rDG==1
for ieq=1:4
    UHLSr=HLSr(Ucurrent(ieq,:),delta_x,delta_y,Coord_center,ESUEL);
    %UGG=GG(Ucurrent(ieq,:),delta_x,delta_y,Coord_center,ESUEL);
    for ie=1:Nelement
        Urect(:,ieq,ie)=UHLSr(:,ie);
    end
end
elseif rDG==2
for ieq=1:4
    UGG=GG(Ucurrent(ieq,:),delta_x,delta_y,Coord_center,ESUEL);
    for ie=1:Nelement
        Urect(:,ieq,ie)=UGG(:,ie);
    end
end
end

%-calculate RHS
    RHS=CRHSP0P1(Ucurrent,Urect,delta_x,delta_y,Coord_center);
%-time advance
    if Timeadvance==1
        for ie=1:Nelement
            Unext(:,ie)=Ucurrent(:,ie)+deltat*RHS(:,ie)/Striangle(1,ie);
        end
    elseif Timeadvance==2
        Uhold=Ucurrent;
    for istage=1:3
        Unhold=afa(istage).*Ucurrent+gamma1(istage).*Uhold+belta(istage)*deltat*(RHS./Striangle);
        Uhold=Unhold;
        if istage~=3
           RHS=CRHSP0P1(Ucurrent,Urect,delta_x,delta_y,Coord_center);
        end
    end
    Unext=Uhold;
    end


%-break cycle
    JUDGE=zeros(Ndim+2,1);
for ie=1:Nelement
    JUDGE=JUDGE+(Unext(:,ie)-Ucurrent(:,ie)).^2*Striangle(1,ie);
end
    judge=max(sqrt(JUDGE));
    if judge<tol
        break;
    end
    nowtime=nowtime+deltat
itimes

%-next level
    Ucurrent=Unext;
    Vcurrent=convert(Ucurrent);
%-next time
    deltat=Deltat(ESUEL,INPOEL,INTFAC,Vcurrent,Ucurrent,Striangle,COORD,BCOND);

if itimes/1000==round(itimes/1000)
rho_point=zeros(1,Npoint);
rhou_point=zeros(1,Npoint);
rhov_point=zeros(1,Npoint);
rhoe_point=zeros(1,Npoint);
v_point=zeros(2,Npoint);
P_point=zeros(1,Npoint);
a_point=zeros(1,Npoint);
for ip=1:Npoint
    x=COORD(1,ip);
    y=COORD(2,ip);
    weight=0;
    for ie2=ESUP2(ip)+1:ESUP2(ip+1)
        ie=ESUP1(ie2);
        weight=weight+Striangle(ie);
        xc=Coord_center(1,ie);
        yc=Coord_center(2,ie);
        phi=zeros(4,1);

        for eq=1:4
        %-calculate uig-
        u_i_g_minus=Urect(1,eq,ie)*(x-xc)+Urect(2,eq,ie)*(y-yc);


       %calculate uig-
        Umax=Ucurrent(eq,ie);
        Umin=Ucurrent(eq,ie);
        for j=1:3
            je=ESUEL(j,ie);
            if je<=Nelement
               if Ucurrent(eq,je)>Umax
                  Umax=Ucurrent(eq,je);
               end

               if Ucurrent(eq,je)<Umin
                  Umin=Ucurrent(eq,je);
               end
            end
        end


        if u_i_g_minus<0
           u_i_g_plus=Umin-Ucurrent(eq,ie);
        else
           u_i_g_plus=Umax-Ucurrent(eq,ie);
        end

        if u_i_g_minus==0
           phi(eq,1)=1;
        else
           phi(eq,1)=min(1,u_i_g_plus/u_i_g_minus);
        end

        end

        rho_point(1,ip)=rho_point(1,ip)+(Ucurrent(1,ie)+phi(1,1)*(Urect(1,1,ie)*(x-Coord_center(1,ie))/delta_x(1,ie)+Urect(2,1,ie)*(y-Coord_center(2,ie))/delta_y(1,ie)))*Striangle(ie);
        rhou_point(1,ip)=rhou_point(1,ip)+(Ucurrent(2,ie)+phi(2,1)*(Urect(1,2,ie)*(x-Coord_center(1,ie))/delta_x(1,ie)+Urect(2,2,ie)*(y-Coord_center(2,ie))/delta_y(1,ie)))*Striangle(ie);
        rhov_point(1,ip)=rhov_point(1,ip)+(Ucurrent(3,ie)+phi(3,1)*(Urect(1,3,ie)*(x-Coord_center(1,ie))/delta_x(1,ie)+Urect(2,3,ie)*(y-Coord_center(2,ie))/delta_y(1,ie)))*Striangle(ie);
        rhoe_point(1,ip)=rhoe_point(1,ip)+(Ucurrent(4,ie)+phi(4,1)*(Urect(1,4,ie)*(x-Coord_center(1,ie))/delta_x(1,ie)+Urect(2,4,ie)*(y-Coord_center(2,ie))/delta_y(1,ie)))*Striangle(ie);
    end
    rho_point(1,ip)=rho_point(1,ip)/weight;
    rhou_point(1,ip)=rhou_point(1,ip)/weight;
    rhov_point(1,ip)=rhov_point(1,ip)/weight;
    rhoe_point(1,ip)=rhoe_point(1,ip)/weight;

end
for ip=1:Npoint
    V=convert([rho_point(1,ip);rhou_point(1,ip);rhov_point(1,ip);rhoe_point(1,ip)]);
    v_point(1:2,ip)=V(2:3,1);
    P_point(1,ip)=V(4,1);
    a_point(1,ip)=V(5,1);    

end
%-calculate error(through S)
t=[1/2,1/2,0;0,1/2,1/2;1/2,0,1/2];W=[1/3,1/3,1/3];%Gauss square
S_far=V_far(4,1)/(V_far(1,1)^gamma);
e=0;
for ie=1:Nelement
    x_1=COORD(1,INPOEL(1,ie));
    x_2=COORD(1,INPOEL(2,ie));
    x_3=COORD(1,INPOEL(3,ie));
    y_1=COORD(2,INPOEL(1,ie));
    y_2=COORD(2,INPOEL(2,ie));
    y_3=COORD(2,INPOEL(3,ie));
    x_trale=[x_1,x_2,x_3];
    y_trale=[y_1,y_2,y_3];
    Square=Striangle(1,ie);
    for i=1:3
       xi=t(i,1)*x_trale(1)+t(i,2)*x_trale(2)+t(i,3)*x_trale(3);
       yi=t(i,1)*y_trale(1)+t(i,2)*y_trale(2)+t(i,3)*y_trale(3);
       B2=(xi-Coord_center(1,ie))/delta_x(1,ie);
       B3=(yi-Coord_center(2,ie))/delta_y(1,ie);
       xc=Coord_center(1,ie);
       yc=Coord_center(2,ie);
        phi=zeros(4,1);

        for eq=1:4
        %-calculate uig-
        u_i_g_minus=Urect(1,eq,ie)*(x-xc)+Urect(2,eq,ie)*(y-yc);


       %calculate uig-
        Umax=Ucurrent(eq,ie);
        Umin=Ucurrent(eq,ie);
        for j=1:3
            je=ESUEL(j,ie);
            if je<=Nelement
               if Ucurrent(eq,je)>Umax
                  Umax=Ucurrent(eq,je);
               end

               if Ucurrent(eq,je)<Umin
                  Umin=Ucurrent(eq,je);
               end
            end
        end


        if u_i_g_minus<0
           u_i_g_plus=Umin-Ucurrent(eq,ie);
        else
           u_i_g_plus=Umax-Ucurrent(eq,ie);
        end

        if u_i_g_minus==0
           phi(eq,1)=1;
        else
           phi(eq,1)=min(1,u_i_g_plus/u_i_g_minus);
        end

        end
       rho_i=Ucurrent(1,ie)+phi(1,1)*(Urect(1,1,ie)*B2+Urect(2,1,ie)*B3);
       rhou_i=Ucurrent(2,ie)+phi(2,1)*(Urect(1,2,ie)*B2+Urect(2,2,ie)*B3);
       rhov_i=Ucurrent(3,ie)+phi(3,1)*(Urect(1,3,ie)*B2+Urect(2,3,ie)*B3);
       p_i=(gamma-1)*(Ucurrent(4,ie)+phi(4,1)*(Urect(1,4,ie)*B2+Urect(2,4,ie)*B3)-0.5*(rhou_i^2/rho_i+rhov_i^2/rho_i));
       S_i=p_i/(rho_i^gamma);
       fi=((S_i-S_far)/S_far)^2;
       e=e+W(i)*fi*Square;
    end
end
e=sqrt(e);



%% Postproceeding
%导出网格
% 打开文件以进行写入
DATA_mix=[COORD',v_point',rho_point',P_point',a_point'];
DATA = fopen('rDGP0P1_data.dat', 'w');
% 使用 fprintf 将数据写入文件
fprintf(DATA,'Variables=x,y,u,v,rho,P,a\n');
fprintf(DATA,'Zone n=');fprintf(DATA,num2str(Npoint));fprintf(DATA,',e=');fprintf(DATA,num2str(Nelement));fprintf(DATA,',f=fepoint,et=triangle\n');
fprintf(DATA,'%f\t%f\t%f\t%f\t%f\t%f\t%f\n',DATA_mix');
fprintf(DATA,'%d\t%d\t%d\n',INPOEL);
% 关闭文件
fclose(DATA);
end

end

rho_point=zeros(1,Npoint);
rhou_point=zeros(1,Npoint);
rhov_point=zeros(1,Npoint);
rhoe_point=zeros(1,Npoint);
v_point=zeros(2,Npoint);
P_point=zeros(1,Npoint);
a_point=zeros(1,Npoint);
for ip=1:Npoint
    x=COORD(1,ip);
    y=COORD(2,ip);
    weight=0;
    for ie2=ESUP2(ip)+1:ESUP2(ip+1)
        ie=ESUP1(ie2);
        weight=weight+Striangle(ie);
        xc=Coord_center(1,ie);
        yc=Coord_center(2,ie);
        phi=zeros(4,1);

        for eq=1:4
        %-calculate uig-
        u_i_g_minus=Urect(1,eq,ie)*(x-xc)+Urect(2,eq,ie)*(y-yc);


       %calculate uig-
        Umax=Ucurrent(eq,ie);
        Umin=Ucurrent(eq,ie);
        for j=1:3
            je=ESUEL(j,ie);
            if je<=Nelement
               if Ucurrent(eq,je)>Umax
                  Umax=Ucurrent(eq,je);
               end

               if Ucurrent(eq,je)<Umin
                  Umin=Ucurrent(eq,je);
               end
            end
        end


        if u_i_g_minus<0
           u_i_g_plus=Umin-Ucurrent(eq,ie);
        else
           u_i_g_plus=Umax-Ucurrent(eq,ie);
        end

        if u_i_g_minus==0
           phi(eq,1)=1;
        else
           phi(eq,1)=min(1,u_i_g_plus/u_i_g_minus);
        end

        end

        rho_point(1,ip)=rho_point(1,ip)+(Ucurrent(1,ie)+phi(1,1)*(Urect(1,1,ie)*(x-Coord_center(1,ie))/delta_x(1,ie)+Urect(2,1,ie)*(y-Coord_center(2,ie))/delta_y(1,ie)))*Striangle(ie);
        rhou_point(1,ip)=rhou_point(1,ip)+(Ucurrent(2,ie)+phi(2,1)*(Urect(1,2,ie)*(x-Coord_center(1,ie))/delta_x(1,ie)+Urect(2,2,ie)*(y-Coord_center(2,ie))/delta_y(1,ie)))*Striangle(ie);
        rhov_point(1,ip)=rhov_point(1,ip)+(Ucurrent(3,ie)+phi(3,1)*(Urect(1,3,ie)*(x-Coord_center(1,ie))/delta_x(1,ie)+Urect(2,3,ie)*(y-Coord_center(2,ie))/delta_y(1,ie)))*Striangle(ie);
        rhoe_point(1,ip)=rhoe_point(1,ip)+(Ucurrent(4,ie)+phi(4,1)*(Urect(1,4,ie)*(x-Coord_center(1,ie))/delta_x(1,ie)+Urect(2,4,ie)*(y-Coord_center(2,ie))/delta_y(1,ie)))*Striangle(ie);
    end
    rho_point(1,ip)=rho_point(1,ip)/weight;
    rhou_point(1,ip)=rhou_point(1,ip)/weight;
    rhov_point(1,ip)=rhov_point(1,ip)/weight;
    rhoe_point(1,ip)=rhoe_point(1,ip)/weight;

end
for ip=1:Npoint
    V=convert([rho_point(1,ip);rhou_point(1,ip);rhov_point(1,ip);rhoe_point(1,ip)]);
    v_point(1:2,ip)=V(2:3,1);
    P_point(1,ip)=V(4,1);
    a_point(1,ip)=V(5,1);    

end
%-calculate error(through S)
t=[1/2,1/2,0;0,1/2,1/2;1/2,0,1/2];W=[1/3,1/3,1/3];%Gauss square
S_far=V_far(4,1)/(V_far(1,1)^gamma);
e=0;
for ie=1:Nelement
    x_1=COORD(1,INPOEL(1,ie));
    x_2=COORD(1,INPOEL(2,ie));
    x_3=COORD(1,INPOEL(3,ie));
    y_1=COORD(2,INPOEL(1,ie));
    y_2=COORD(2,INPOEL(2,ie));
    y_3=COORD(2,INPOEL(3,ie));
    x_trale=[x_1,x_2,x_3];
    y_trale=[y_1,y_2,y_3];
    Square=Striangle(1,ie);
    for i=1:3
       xi=t(i,1)*x_trale(1)+t(i,2)*x_trale(2)+t(i,3)*x_trale(3);
       yi=t(i,1)*y_trale(1)+t(i,2)*y_trale(2)+t(i,3)*y_trale(3);
       B2=(xi-Coord_center(1,ie))/delta_x(1,ie);
       B3=(yi-Coord_center(2,ie))/delta_y(1,ie);
       xc=Coord_center(1,ie);
       yc=Coord_center(2,ie);
        phi=zeros(4,1);

        for eq=1:4
        %-calculate uig-
        u_i_g_minus=Urect(1,eq,ie)*(x-xc)+Urect(2,eq,ie)*(y-yc);


       %calculate uig-
        Umax=Ucurrent(eq,ie);
        Umin=Ucurrent(eq,ie);
        for j=1:3
            je=ESUEL(j,ie);
            if je<=Nelement
               if Ucurrent(eq,je)>Umax
                  Umax=Ucurrent(eq,je);
               end

               if Ucurrent(eq,je)<Umin
                  Umin=Ucurrent(eq,je);
               end
            end
        end


        if u_i_g_minus<0
           u_i_g_plus=Umin-Ucurrent(eq,ie);
        else
           u_i_g_plus=Umax-Ucurrent(eq,ie);
        end

        if u_i_g_minus==0
           phi(eq,1)=1;
        else
           phi(eq,1)=min(1,u_i_g_plus/u_i_g_minus);
        end

        end
       rho_i=Ucurrent(1,ie)+phi(1,1)*(Urect(1,1,ie)*B2+Urect(2,1,ie)*B3);
       rhou_i=Ucurrent(2,ie)+phi(2,1)*(Urect(1,2,ie)*B2+Urect(2,2,ie)*B3);
       rhov_i=Ucurrent(3,ie)+phi(3,1)*(Urect(1,3,ie)*B2+Urect(2,3,ie)*B3);
       p_i=(gamma-1)*(Ucurrent(4,ie)+phi(4,1)*(Urect(1,4,ie)*B2+Urect(2,4,ie)*B3)-0.5*(rhou_i^2/rho_i+rhov_i^2/rho_i));
       S_i=p_i/(rho_i^gamma);
       fi=((S_i-S_far)/S_far)^2;
       e=e+W(i)*fi*Square;
    end
end
e=sqrt(e);



%% Postproceeding
%导出网格
% 打开文件以进行写入
INPOEL=INPOEL';
DATA_mix=[COORD',v_point',rho_point',P_point',a_point'];
DATA = fopen('rDGP0P1_data.dat', 'w');
% 使用 fprintf 将数据写入文件
fprintf(DATA,'Variables=x,y,u,v,rho,P,a\n');
fprintf(DATA,'Zone n=');fprintf(DATA,num2str(Npoint));fprintf(DATA,',e=');fprintf(DATA,num2str(Nelement));fprintf(DATA,',f=fepoint,et=triangle\n');
fprintf(DATA,'%f\t%f\t%f\t%f\t%f\t%f\t%f\n',DATA_mix');
fprintf(DATA,'%d\t%d\t%d\n',INPOEL');
% 关闭文件
fclose(DATA);
OK=[1,itimes,log10(e),log10(min(Striangle))];

end


