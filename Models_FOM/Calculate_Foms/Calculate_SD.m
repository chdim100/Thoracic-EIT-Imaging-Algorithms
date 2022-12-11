function SD=Calculate_SD(polyelement,llcand,rlcand,llunginf,rlunginf,llung_rec_area,rlung_rec_area)

if size(polyelement,3)==1024||size(polyelement,3)==1060
    scaling=148;
    factorA=1;
    factorB=0;
else
    scaling=1.479;
    factorA=2;
    factorB=1;
end

%%%%shape deformation
%find right lung reconstructed elements that are outside the reference
%shape and their corresponding area
Aout_rl=0;
for irl=1:length(rlcand)
    pelem=polyshape([polyelement(:,1,factorA*rlcand(irl)-factorB)...
        polyelement(:,2,factorA*rlcand(irl)-factorB)]/scaling);
    porg=polyshape(rlunginf(:,1)/scaling,rlunginf(:,2)/scaling);
    wrl_out(irl)=1-area(intersect(pelem, porg))/area(pelem);
    Aout_rl=Aout_rl+factorA*area(pelem)*wrl_out(irl);
end
SD.Right_Lung=Aout_rl/rlung_rec_area;

%find left lung reconstructed elements that are outside the reference
%shape and their corresponding area
Aout_ll=0;
for ill=1:length(llcand)
    pelem=polyshape([polyelement(:,1,factorA*llcand(ill)-factorB)...
        polyelement(:,2,factorA*llcand(ill)-factorB)]/scaling);
    porg=polyshape(llunginf(:,1)/scaling,llunginf(:,2)/scaling);
    wll_out(ill)=1-area(intersect(pelem, porg))/area(pelem);
    Aout_ll=Aout_ll+factorA*area(pelem)*wll_out(ill);
end
SD.Left_Lung=Aout_ll/llung_rec_area;
SD.Total=(Aout_ll+Aout_rl)/(llung_rec_area+rlung_rec_area);

% o=[];

end