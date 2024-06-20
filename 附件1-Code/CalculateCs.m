clc
clear
close all
%% Preproceeding
%input data
%% Input data
addpath(genpath('D:/CFD课题组/研一下/Project2/CFD_Project2_王程/Mesh_2_通过NACA0012 翼型的跨音速流/DGP0_data'))
addpath(genpath('D:/CFD课题组/研一下/Project2/CFD_Project2_王程/Mesh_2_通过NACA0012 翼型的跨音速流/DGP0P1_data'))
addpath(genpath('D:/CFD课题组/研一下/Project2/CFD_Project2_王程/Mesh_2_通过NACA0012 翼型的跨音速流/DGP1_data'))
addpath(genpath('D:/CFD课题组/研一下/Project2/CFD_Project2_王程/附件1-Code/pro2_mesh'))
global Ndim;
global Nelement;
global Npoint;
global Nnode;
global NBface;
global Nface;
%Input mesh
    feflo_bump = fopen('pro2_mesh/feflo.domn.naca0012', 'rb');%open the file
    
skiprows=fscanf(feflo_bump, '%lf',1);% skip some rows
for i = 1:skiprows+2
    line=fgetl(feflo_bump);
end
Ndata = fscanf(feflo_bump, '%lf',[2,1]);
%read dimention&Nnode
Ndim=Ndata(1,1);Nnode=Ndata(2,1);
for i = 1:2
    line=fgetl(feflo_bump);
end
Ndata = fscanf(feflo_bump, '%lf',[3,1]);
Nelement=Ndata(1,1);Npoint=Ndata(2,1);NBface=Ndata(3,1);%Information of elements&points&boundary

for i = 1:3
    line=fgetl(feflo_bump);
end
%read data per row
C = textscan(line, '%f');
numbers = length(str2num(line));
INPOEL = fscanf(feflo_bump, '%lf',[numbers,Nelement]);
INPOEL=[C{1,1},INPOEL];
INPOEL=INPOEL(2:4,:);

%-store the coordinates of the points?COORD
for i = 1:1% skip some rows
    line=fgetl(feflo_bump);
end
COORD = fscanf(feflo_bump, '%lf',[Ndim+1,Npoint]);
COORD=COORD(2:3,:);
%-store Boundary conditions
for i = 1:2+Npoint+2
    line=fgetl(feflo_bump);
end
BCOND=zeros(3,NBface);
for i=1:NBface
    C = textscan(line, '%f');
    BCOND(:,i)=C{1,1}(2:4,1);
    line=fgetl(feflo_bump);
end

fclose(feflo_bump);
 
%-Elements surrounding points
ESUP2=zeros(1,Npoint+1);
for ie=1:Nelement
    for in=1:Nnode
        ESUP2(INPOEL(in,ie)+1)=ESUP2(INPOEL(in,ie)+1)+1;  
    end
end

for IPOIN=2:Npoint+1
    ESUP2(IPOIN)=ESUP2(IPOIN)+ESUP2(IPOIN-1);
end
MESUP=ESUP2(Npoint+1);
ESUP1=zeros(1,MESUP);
for ie=1:Nelement
    for in=1:Nnode
        IPOIN=INPOEL(in,ie);
        ISTOR=ESUP2(IPOIN)+1;
        ESUP2(IPOIN)=ISTOR;
        ESUP1(ISTOR)=ie;
    end
end
for IPOIN=Npoint+1:-1:2
    ESUP2(IPOIN)=ESUP2(IPOIN-1);
end
ESUP2(1)=0;

%Elements surrounding elements
ESUEL=zeros(3,Nelement);
LPOIN=zeros(1,Npoint);
%-Auxiliary array
LNOFA=zeros(2,3);
LNOFA(1,1)=2;LNOFA(2,1)=3;
LNOFA(1,2)=3;LNOFA(2,2)=1;
LNOFA(1,3)=1;LNOFA(2,3)=2;

for ie=1:Nelement
    for iface=1:Nnode
        LHELP=INPOEL(LNOFA(:,iface),ie);
        LPOIN(1,LHELP(1,1))=1;LPOIN(1,LHELP(2,1))=1;
        IPOIN=INPOEL(LNOFA(1,iface),ie);%select a point of the face
        FLAG=0;
        for ISTOR=ESUP2(IPOIN)+1:ESUP2(IPOIN+1)
            JELEM=ESUP1(ISTOR);
            if(JELEM~=ie)
                for JFAEL=1:Nnode%loop over the faces of element JELEM
                    ICOUN=0;
                    for JOFA=1:Nnode-1
                        JPOIN=INPOEL(LNOFA(JOFA,JFAEL),JELEM);
                        if(LPOIN(1,JPOIN)==1)
                            ICOUN=ICOUN+1;
                        end
                    end
                        if(ICOUN==Nnode-1)
                            ESUEL(iface,ie)=JELEM;
                            ESUEL(JFAEL,JELEM)=ie;
                            FLAG=1;
                        end
                 end
            end
        end

        if(FLAG==0)
            for IFACE=1:NBface
                if((BCOND(1,IFACE)==LHELP(1,1)&&BCOND(2,IFACE)==LHELP(2,1))||(BCOND(2,IFACE)==LHELP(1,1)&&BCOND(1,IFACE)==LHELP(2,1)))
                    ESUEL(iface,ie)=IFACE+Nelement;
                    break;
                end
            end
        end
        LPOIN=zeros(1,Npoint);
    end
end


Nface=NBface;%Numbers of all faces
for ie=1:Nelement
    for IFAEL=1:Nnode
        JELEM=ESUEL(IFAEL,ie);
        if(JELEM>ie&&JELEM<=Nelement)
            Nface=Nface+1;
        end
    end
end

NIFFA=7;%number of information of face
INTFAC=zeros(NIFFA,Nface);%interface array 
iface=1;
for ie=1:Nelement
    for IFAEL=1:Nnode
        JELEM=ESUEL(IFAEL,ie);
        if(JELEM>ie)
            INTFAC(2,iface)=JELEM;%store right element
            INTFAC(1,iface)=ie;%store left element
            INTFAC(3,iface)=INPOEL(LNOFA(1,IFAEL),ie);%store left point of face
            INTFAC(4,iface)=INPOEL(LNOFA(2,IFAEL),ie);%store right point of face
            if(JELEM>Nelement)
                INTFAC(5,iface)=BCOND(3,JELEM-Nelement);%judge whether this face is fixed or far or internal
            else
                INTFAC(5,iface)=-1;
            end
            iel=ie;ier=JELEM;
            %calculate V_n
            ip1=INPOEL(LNOFA(1,IFAEL),ie);%points in face
            ip2=INPOEL(LNOFA(2,IFAEL),ie);
             for IPL=1:Nnode
                 if(INPOEL(IPL,iel)~=ip1&&INPOEL(IPL,iel)~=ip2)
                     ip3=INPOEL(IPL,iel);%another point in iel
                 end
             end
             Judge_vector1=[COORD(1,ip2)-COORD(1,ip1),COORD(2,ip2)-COORD(2,ip1)];
             Judge_vector2=[COORD(1,ip3)-COORD(1,ip1),COORD(2,ip3)-COORD(2,ip1)];             
             V_a=COORD(1,ip2)-COORD(1,ip1);V_b=COORD(2,ip2)-COORD(2,ip1);
             if V_b==0
                  V_product1=[0,1];V_product2=[0,-1];
             else
                  V_product1=[-V_b,V_a]/norm([-V_b,V_a]);
                  V_product2=[V_b,-V_a]/norm([V_b,-V_a]);
             end

             if dot(Judge_vector2,V_product1)<0
                 V_n=V_product1;
             elseif dot(Judge_vector2,V_product2)<0
                 V_n=V_product2;
             end
             INTFAC(6,iface)=V_n(1,1);
             INTFAC(7,iface)=V_n(1,2);
             tau_ij=sqrt((COORD(1,ip2)-COORD(1,ip1))^2+(COORD(2,ip2)-COORD(2,ip1))^2);
             INTFAC(8,iface)=tau_ij;

            iface=iface+1;
        end
    end
end

%calculate squares of elements
Striangle=zeros(1,Nelement);
for ie=1:Nelement
    ipoin1=INPOEL(1,ie);ipoin2=INPOEL(2,ie);ipoin3=INPOEL(3,ie);
    HPcycle=[ipoin1,ipoin2,ipoin3,ipoin1,ipoin2];c_ie=zeros(1,3);
    for k=1:Nnode
           c_ie(k)=COORD(1,HPcycle(k+1))*COORD(2,HPcycle(k+2))-COORD(1,HPcycle(k+2))*COORD(2,HPcycle(k+1));
    end
    Striangle(1,ie)=0.5*sum(c_ie);
end

%input rho&u&v&p in point
feflo_bump = fopen('DGP0_data.dat', 'rb');%open the file
for i = 1:2
    line=fgetl(feflo_bump);
end
Data= fscanf(feflo_bump, '%lf',[7,Npoint]);
U_point=Data(3,:);
V_point=Data(4,:);
Rho_point=Data(5,:);
P_point=Data(6,:);
A_point=Data(7,:);
fclose(feflo_bump);

gamma=1.4;Ma_far=0.85;afa_attack=pi/180;
P_far=1/(gamma*Ma_far^2);
rho_far=1;
u_far=cos(afa_attack);v_far=sin(afa_attack);
S_far=P_far/rho_far^gamma;
%calculate Coefficient p
Cp=zeros(3,Npoint);
% for bface=1:NBface
%     if BCOND(3,bface)==2
% xc=COORD(1,BCOND(1,1));x_min=xc;x_max=xc;
% yc=COORD(2,BCOND(1,1));y_min=yc;y_max=yc;
% break;
%     end
% end
for bface=1:NBface
    if BCOND(3,bface)==2
    ip1=BCOND(1,bface);
    ip2=BCOND(2,bface);
%     if x_min>COORD(1,ip1)
%         x_min=COORD(1,ip1);
%         y_min=COORD(2,ip1);
%     end
%     if x_max<COORD(1,ip1)
%         x_max=COORD(1,ip1);
%         y_max=COORD(2,ip1);
%     end 
%if COORD(1,ip1)>=1&&COORD(1,ip1)<=2&&COORD(2,ip1)>=0&&COORD(2,ip1)<1
if COORD(1,ip1)>=0&&COORD(1,ip1)<=1
    Cp(1,ip1)=COORD(1,ip1);
    Cp(2,ip1)=COORD(2,ip1);
    Cp(3,ip1)=(P_point(1,ip1)/Rho_point(1,ip1)^gamma-S_far)/S_far;
end
%if COORD(1,ip2)>=1&&COORD(1,ip2)<=2&&COORD(2,ip2)>=0&&COORD(2,ip2)<1
if COORD(1,ip2)>=0&&COORD(1,ip2)<=1
    Cp(1,ip2)=COORD(1,ip2);
    Cp(2,ip2)=COORD(2,ip2);
    Cp(3,ip2)=(P_point(1,ip2)/Rho_point(1,ip2)^gamma-S_far)/S_far;
end
    end
end

% c=sqrt((x_min-x_max)^2+(y_min-y_max)^2);
c=1;
k=1;j=1;
Cp_plus=zeros(2,NBface+1);
Cp_minus=zeros(2,NBface+1);
for ip=1:Npoint
    if sum(Cp(:,ip))~=0&&Cp(2,ip)>=0
        Cp_plus(:,k)=Cp([1,3],ip);
        k=k+1;
    end
    if sum(Cp(:,ip))~=0&&Cp(2,ip)<0
        Cp_minus(:,j)=Cp([1,3],ip);
        j=j+1;
    end    
end
Cp_plus(:,k:NBface+1)=[];
Cp_minus(:,j:NBface+1)=[];
Cp_plus(1,:)=Cp_plus(1,:)/c;
Cp_minus(1,:)=Cp_minus(1,:)/c;
length_plus=size(Cp_plus,2);
length_minus=size(Cp_minus,2);
for ip1=1:length_plus-1
    for ip2=ip1+1:length_plus
    if Cp_plus(1,ip1)>Cp_plus(1,ip2)
        Cpt=Cp_plus(:,ip1);
        Cp_plus(:,ip1)=Cp_plus(:,ip2);
        Cp_plus(:,ip2)=Cpt;
    end
    end
end

for ip1=1:length_minus-1
    for ip2=ip1+1:length_minus
    if Cp_minus(1,ip1)>Cp_minus(1,ip2)
        Cpt=Cp_minus(:,ip1);
        Cp_minus(:,ip1)=Cp_minus(:,ip2);
        Cp_minus(:,ip2)=Cpt;
    end
    end
end

H1=plot(Cp_plus(1,:),Cp_plus(2,:),'-*b','Linewidth',1.5);
hold on
H2=plot(Cp_minus(1,:),Cp_minus(2,:),'-*r','Linewidth',1.5);
ylim([-0.06,0.05]);
lgd=legend([H1,H2],'Arc protrusions','Under the wing');
lgd.FontSize=12;
xlabel('x/c','FontSize',15);
ylabel('Coefficient S','FontSize',15);     
grid on
hold off




