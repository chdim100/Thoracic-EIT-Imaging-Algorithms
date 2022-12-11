function [RES,llung_rec_area,rlung_rec_area]=Calculate_RES(polyelement,llcand,rlcand)

if size(polyelement,3)==1024||size(polyelement,3)==1060
    scaling=148;
else
    scaling=1.479;
end

%%%total medium area
Ao=0;
for element=1:size(polyelement,3)
    Ao=Ao+area(polyshape([polyelement(:,1,element) polyelement(:,2,element)]/scaling));
end

%%%%calculate resolution
%%%%left lung
llung_rec_area=0;
if size(polyelement,3)==1024||size(polyelement,3)==1060
    for ill=1:length(llcand)
        llung_rec_area=llung_rec_area+area(polyshape([polyelement(:,1,llcand(ill)) polyelement(:,2,llcand(ill))]/scaling));
    end
else
    for ill=1:length(llcand)
        llung_rec_area=llung_rec_area+2*area(polyshape([polyelement(:,1,2*llcand(ill)-1) polyelement(:,2,2*llcand(ill)-1)]/scaling));
    end
end

RES.Left_Lung=sqrt(llung_rec_area/Ao);

%%%%calculate resolution
%%%%right lung
rlung_rec_area=0;
if size(polyelement,3)==1024||size(polyelement,3)==1060
    for irl=1:length(rlcand)
        rlung_rec_area=rlung_rec_area+area(polyshape([polyelement(:,1,rlcand(irl)) polyelement(:,2,rlcand(irl))]/scaling));
    end
else
    for irl=1:length(rlcand)
        rlung_rec_area=rlung_rec_area+2*area(polyshape([polyelement(:,1,2*rlcand(irl)-1) polyelement(:,2,2*rlcand(irl)-1)]/scaling));
    end
end

RES.Right_Lung=sqrt(rlung_rec_area/Ao);
RES.Total=sqrt((rlung_rec_area+llung_rec_area)/Ao);


end