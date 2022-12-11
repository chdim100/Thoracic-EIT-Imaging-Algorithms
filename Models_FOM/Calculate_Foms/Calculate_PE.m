function [PE,RL_rec_centre,LL_rec_centre]=Calculate_PE(polyelement,llcand,rlcand,Left_lung_ref_centre, Right_lung_ref_centre,centrelement)

if size(polyelement,3)==1024||size(polyelement,3)==1060
    scaling=148;
else
    scaling=1.479;
end

if size(polyelement,3)==1024||size(polyelement,3)==1060
    llx3=squeeze(polyelement(:,1,llcand(:)));
    lly3=squeeze(polyelement(:,2,llcand(:)));
    llx3_vec=llx3(:);
    lly3_vec=lly3(:);
    llx=unique(llx3_vec,'stable');
    lly=unique(lly3_vec,'stable');
    LL_rec_centre=[mean(llx),mean(lly)];
else
    LL_rec_centre=[mean(centrelement(llcand,1)) mean(centrelement(llcand,2))];
end
%%%%Position Error (PE)

PE.left_lung=norm(Left_lung_ref_centre-LL_rec_centre)/scaling;
%sqrt(sum((abs(Left_lung_ref_centre-LL_rec_centre)/scaling).^2));

%%%%right lung target
if size(polyelement,3)==1024||size(polyelement,3)==1060
    rlx3=squeeze(polyelement(:,1,rlcand(:))); rlx3_vec=rlx3(:);
    rly3=squeeze(polyelement(:,2,rlcand(:))); rly3_vec=rly3(:);
    rlx=unique(rlx3_vec,'stable');
    rly=unique(rly3_vec,'stable');
    RL_rec_centre=[mean(rlx),mean(rly)];
else
    RL_rec_centre=[mean(centrelement(rlcand,1)) mean(centrelement(rlcand,2))];
end
%%%%Position Error (PE)

PE.right_lung=norm(Right_lung_ref_centre-RL_rec_centre)/scaling;
%sqrt(sum((abs(Right_lung_ref_centre-RL_rec_centre)/scaling).^2));
PE.Total=PE.right_lung+PE.left_lung;

end