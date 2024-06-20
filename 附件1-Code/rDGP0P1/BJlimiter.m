function phi=BJlimiter(iel,U,Urect,Coord_face,Coord_center,ESUEL)
global Nelement;
ig=2;
tf=[-sqrt(15)/5,0,sqrt(15)/5];Wf=[5/9,8/9,5/9];%-Gauss in face
x1=Coord_face(1,1);x2=Coord_face(1,2);
y1=Coord_face(2,1);y2=Coord_face(2,2);
xi=x1*0.5*(1-tf(ig))+x2*0.5*(1+tf(ig));
yi=y1*0.5*(1-tf(ig))+y2*0.5*(1+tf(ig));
xc=Coord_center(1,iel);
yc=Coord_center(2,iel);
phi=zeros(4,1);

for eq=1:4
%-calculate uig-
u_i_g_minus=Urect(1,eq,iel)*(xi-xc)+Urect(2,eq,iel)*(yi-yc);


%calculate uig-
Umax=U(eq,iel);
Umin=U(eq,iel);
for j=1:3
    je=ESUEL(j,iel);
    if je<=Nelement
        if U(eq,je)>Umax
            Umax=U(eq,je);
        end

        if U(eq,je)<Umin
            Umin=U(eq,je);
        end
    end
end


if u_i_g_minus<0
    u_i_g_plus=Umin-U(eq,iel);
else
    u_i_g_plus=Umax-U(eq,iel);
end

if u_i_g_minus==0
    phi(eq,1)=1;
else
    phi(eq,1)=min(1,u_i_g_plus/u_i_g_minus);
end

end

end