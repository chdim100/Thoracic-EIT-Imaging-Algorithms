function RNG=Calculate_RNG(polyelement,llcand,rlcand,image)

if size(polyelement,3)==1060
image.elem_data=image.elem_data-1;
end

%%%Count total lung values in the reference area
SXRL=0;
for irl=1:length(rlcand)
    SXRL=SXRL+real(image.elem_data(rlcand(irl)));
end

SXLL=0;
for ill=1:length(llcand)
    SXLL=SXLL+real(image.elem_data(llcand(ill)));
end
%%%Count Ringing
SRING=0;
if size(polyelement,3)==1024||size(polyelement,3)==1060
    for ies=1:size(polyelement,3)
        if ~ismember(ies,rlcand)&&~ismember(ies,llcand)&&real(image.elem_data(ies))>0
            SRING=SRING+real(image.elem_data(ies));
        end
    end
else
    for ies=1:size(polyelement,3)/2
    if ~ismember(ies,rlcand)&&~ismember(ies,llcand)&&real(image.elem_data(ies))>0
        SRING=SRING+real(image.elem_data(ies));
    end
    end
end
RNG=-SRING/(SXRL+SXLL);

end