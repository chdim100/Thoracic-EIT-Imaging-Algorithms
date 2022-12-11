function [RLind,RLind_c,LLind,LLind_c]=find_lungs_indices(polyelement,llunginf,rlunginf)
%%%returns the element indices that are in the reference image lungs
if size(polyelement,3)==1024||size(polyelement,3)==1060
    scaling=148;
    factorA=1;
    factorB=0;
else
    scaling=1.479;
    factorA=2;
    factorB=1;
end

%find right lung reconstructed elements that are inside the reference
%shape and their corresponding area
for ii=1:size(polyelement,3)
    pelem=polyshape([polyelement(:,1,ii)...
        polyelement(:,2,ii)]/scaling);
    porg=polyshape(rlunginf(:,1)/scaling,rlunginf(:,2)/scaling);
    wrl(ii)=area(intersect(pelem, porg))/area(pelem);
end

RLind=find(wrl>=0.5);
if size(polyelement,3)==2052
    wrl=wrl(1:2:end);
    RLind_c=find(wrl>=0.5);
else
    RLind_c=RLind;
end

%find left lung reconstructed elements that are inside the reference
%shape and their corresponding area
for ii=1:size(polyelement,3)
    pelem=polyshape([polyelement(:,1,ii)...
        polyelement(:,2,ii)]/scaling);
    porg=polyshape(llunginf(:,1)/scaling,llunginf(:,2)/scaling);
    wll(ii)=area(intersect(pelem, porg))/area(pelem);
end
LLind=find(wll>=0.5);
if size(polyelement,3)==2052
    wll=wll(1:2:end);
    LLind_c=find(wll>=0.5);
else
    LLind_c=LLind;
end


end