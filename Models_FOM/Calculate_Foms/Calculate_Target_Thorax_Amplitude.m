function [TA,image,filtered_image]=Calculate_Target_Thorax_Amplitude(image,polyelement)

%%%%%%threshold 25% of the change

image.filtered_elem_data=image.elem_data;

%%%%if M.o.M., remove backround
if min(image.filtered_elem_data)>0
    image.filtered_elem_data=image.filtered_elem_data-1;
end

image.filtered_elem_data(image.filtered_elem_data...
    >-1/4*max(abs(image.filtered_elem_data)))=0;
image.filtered_elem_data(image.filtered_elem_data...
    <=-1/4*max(abs(image.filtered_elem_data)))=1;
filtered_image=image;
filtered_image.elem_data=image.filtered_elem_data;
%%%check if GREIT pixel domain is activated
if size(polyelement,3)==length(filtered_image.elem_data)
    filtered_image=remove_bound_artefacts(filtered_image,polyelement);
end

 if min(real(image.elem_data))<0
     TA=sum(real(image.elem_data));
 else
    TA=sum(real(image.elem_data)-1);
 end


end

function filt_im=remove_bound_artefacts(filtered_image,polyelement)
filt_im=filtered_image;

if size(polyelement,1)==3
    %%%F.E.M. (GN/GREIT)
boundnodeinx=filtered_image.fwd_model.boundary(:);

for inn=1:length(boundnodeinx)
    for element=1:size(polyelement,3)
        if ismember(boundnodeinx(inn),filtered_image.fwd_model.elems(element,:))
            filt_im.elem_data(element)=0;
        end
    end
end

else
    %%%%M.o.M.
    A=mean(polyelement(:,1,:),1);
    B=mean(polyelement(:,2,:),1);
    
    kk=boundary(squeeze(A),squeeze(B));
    filt_im.elem_data(kk)=0;
end
end
