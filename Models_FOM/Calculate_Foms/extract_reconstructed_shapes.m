function [llcand,rlcand,centrelement]=extract_reconstructed_shapes(filtered_image,polyelement)
%%%%inputs: -filtered_image: a binary image filtered with the 1/4-amplitude threshold
%(see Calculate_Target_Thorax_Amplitude) and with the boundary artifacts
%removed
%%%%       -polyelement: a complete element-array (Nodes/element)xDimxNo_of_Elements

%%%%outputs: rlcand and llcand: the elements of the reconstructed image
%%%%that correspond to the right and left reconstructed lung respectively.

%%%%Categorize the elements to left and right
rightelements=[];
leftelements=[];
if size(polyelement,3)==1024||size(polyelement,3)==1060  %%%%FEM/M.o.M.
    
    for element=1:size(polyelement,3)
        [centrelement(element,1),centrelement(element,2)]=centroid(polyshape(polyelement(:,1,element),polyelement(:,2,element)));
        if centrelement(element,1)<0
            rightelements=[rightelements element];
        else
            leftelements=[leftelements element];
        end
    end
    
else  %%%%GREIT
    
    for element=1:size(polyelement,3)/2
        [x1,y1]=centroid(polyshape(polyelement(:,1,2*element-1),polyelement(:,2,2*element-1)));
        [x2,y2]=centroid(polyshape(polyelement(:,1,2*element),polyelement(:,2,2*element)));
        centrelement(element,1)=(x1+x2)/2;
        centrelement(element,2)=(y1+y2)/2;
        if centrelement(element,1)<0
            rightelements=[rightelements element];
        else
            leftelements=[leftelements element];
        end
    end
    
end

rlcand=[];
llcand=[];
for element=1:length(filtered_image.elem_data)
    if filtered_image.elem_data(element)==1
        if ismember(element,rightelements)
            rlcand=[rlcand element];
        else
            llcand=[llcand element];
        end
    end
end
o=[];
end