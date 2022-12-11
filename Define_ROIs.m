function [Elems_ROI,LROI,L]=Define_ROIs(N,ROItype)

eidors_msg('Defining ROIs');

invmodel=mk_common_model('d2T3',N);

invmodel.fwd_solve.get_all_meas = 1;

elements=invmodel.fwd_model.elems;
L=size(elements,1);
elementcentre=zeros(L,2);
for element=1:1024
    nodes=[invmodel.fwd_model.elems(element,1) invmodel.fwd_model.elems(element,2)...
        invmodel.fwd_model.elems(element,3)];
    elementcentre(element,1)=(invmodel.fwd_model.nodes(nodes(1),1)+...
        invmodel.fwd_model.nodes(nodes(2),1)+invmodel.fwd_model.nodes(nodes(3),1))/3;
    elementcentre(element,2)=(invmodel.fwd_model.nodes(nodes(1),2)+...
        invmodel.fwd_model.nodes(nodes(2),2)+invmodel.fwd_model.nodes(nodes(3),2))/3;
end
clear nodes

switch ROItype
    
    case{'Lungs and Heart'}
        %%%right lung ROI
        aRL=45; % horizontal radius
        bRL=140; % vertical radius
        x0RL=-85; % x0,y0 ellipse centre coordinates
        y0RL=0;
        t=-pi:0.01:pi;
        xRL=(x0RL+aRL*cos(t))';
        yRL=(y0RL+bRL*sin(t))';
        elements_in_right_lung = inpolygon(elementcentre(:,1), elementcentre(:,2),...
            [xRL; xRL(1)],[yRL; yRL(1)]);
        ROI_elements.elements_in_right_lung=elements_in_right_lung;
        
        %%%left lung ROI
        aLL=45; % horizontal radius
        bLL=140; % vertical radius
        x0LL=85; % x0,y0 ellipse centre coordinates
        y0LL=0;
        t=-pi:0.01:pi;
        xLL=(x0LL+aLL*cos(t))';
        yLL=(y0LL+bLL*sin(t))';
        elements_in_left_lung = inpolygon(elementcentre(:,1), elementcentre(:,2),...
            [xLL; xLL(1)],[yLL; yLL(1)]);
        ROI_elements.elements_in_left_lung=elements_in_left_lung;
        
        %%%Heart ROI
        aH=60; % horizontal radius
        bH=90; % vertical radius
        x0H=0; % x0,y0 ellipse centre coordinates
        y0H=-20;
        t=-pi:0.01:pi;
        xH=(x0H+aH*cos(t))';
        yH=(y0H+bH*sin(t))';
        elements_in_heart = inpolygon(elementcentre(:,1), elementcentre(:,2),...
            [xH; xH(1)],[yH; yH(1)]);
        ROI_elements.elements_in_heart=elements_in_heart;
        
        Elems_ROI=ROI_elements.elements_in_left_lung|ROI_elements.elements_in_right_lung|...
            ROI_elements.elements_in_heart;
        
    case{'central','Central'}
        %%%Central ROI
        aC=100; % horizontal radius
        bC=100; % vertical radius
        x0C=0; % x0,y0 ellipse centre coordinates
        y0C=0;
        t=-pi:0.01:pi;
        xC=(x0C+aC*cos(t))';
        yC=(y0C+bC*sin(t))';
        elements_R = inpolygon(elementcentre(:,1), elementcentre(:,2),...
            [xC; xC(1)],[yC; yC(1)]);
        ROI_elements.elements_R=elements_R;
        Elems_ROI=elements_R;
        
    case{'All','Omega'}
        Elems_ROI=boolean(ones(L,1));
        
    otherwise
        error('Non defined ROI selection of this type')
end

LROI=sum(Elems_ROI);

figure
scatter(elementcentre(:,1),elementcentre(:,2),85,'b','filled')
hold on
scatter(elementcentre(Elems_ROI,1),elementcentre(Elems_ROI,2),85,'r','filled')
title('ROI selection')
Leg=legend({'\Omega','ROI'});
axis off
axis square

end