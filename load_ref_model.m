function [inv_reference]=load_ref_model2(type,invmodel)

%run C:\Users\Chris\Desktop\PHD_files\eidors-v3.9-ng\eidors\startup.m

cavity=5;

invref=mk_image(invmodel,1);
elem_index=invref.fwd_model.elems;
%%%% polyelement: (NODES==3) X 2
for element=1:size(elem_index,1)
    polyelement(:,1,element)=[invref.fwd_model.nodes(elem_index(element,1),1) invref.fwd_model.nodes(elem_index(element,2),1)...
        invref.fwd_model.nodes(elem_index(element,3),1)];
    polyelement(:,2,element)=[invref.fwd_model.nodes(elem_index(element,1),2) invref.fwd_model.nodes(elem_index(element,2),2)...
        invref.fwd_model.nodes(elem_index(element,3),2)];
end


switch type
    case{'Thorax3D'}
        f=100000;
        addpath('C:\Users\Chris\Desktop\PHD_files\eidors-v3.9-ng\eidors\EIT_Circuit_Sim\Structures\Thoracic_New\')
        addpath('C:\Users\Chris\Desktop\PHD_files\eidors-v3.9-ng\eidors\EIT_Circuit_Sim\Structures\Thoracic_New\freq_elements')
        addpath('D:\EIT\absolute_imaging\Models_FOM\Thorax_3D\')
        addpath('C:\Users\Chris\Desktop\PHD_files\eidors-v3.9-ng\eidors\EIT_Circuit_Sim\')
        addpath('C:\Users\Chris\Desktop\PHD_files\eidors-v3.9-ng\eidors\EIT_Circuit_Sim\Tissue_Properties\Thorax\')
        addpath('D:\EIT\absolute_imaging\Thoracic_Cavities\Thorax1')
        addpath('D:\EIT\absolute_imaging\Thoracic_Cavities\Thorax2')
        %%%%%%parse all data
        load(['Tissues_' num2str(round(f/1000)) 'kHz.mat'])
        load('Lung_med_values.mat')
        load(['thorax' num2str(cavity) 'inf.mat'])
        load(['thorax' num2str(cavity) '.mat'])
        load(['spondilus' num2str(cavity) '.mat'])
        load(['rlung' num2str(cavity) 'inf.mat'])
        load(['rlung' num2str(cavity) '.mat'])
        load(['llung' num2str(cavity) 'inf.mat'])
        load(['llung' num2str(cavity) '.mat'])
        load(['newheart' num2str(cavity) '.mat'])
        load(['spondilus' num2str(cavity) '.mat'])
        load(['llung' num2str(cavity) '_medstate.mat'])
        load(['rlung' num2str(cavity) '_medstate.mat'])
        load(['thorax' num2str(cavity) '_medstate.mat'])
        
        %%%assign to global names (non-cavity dependent)
        thorax=eval(['thorax' num2str(cavity)]);
        llung=eval(['llung' num2str(cavity)]);
        rlung=eval(['rlung' num2str(cavity)]);
        newheart=eval(['newheart' num2str(cavity)]);
        spondilus=eval(['spondilus' num2str(cavity)]);
        llungstate2=squeeze(llung_medstate(2,:,:));
        rlungstate2=squeeze(rlung_medstate(2,:,:));
        thoraxstate2=squeeze(thorax_medstate(2,:,:));
        llungstate3=squeeze(llung_medstate(3,:,:));
        rlungstate3=squeeze(rlung_medstate(3,:,:));
        thoraxstate3=squeeze(thorax_medstate(3,:,:));
        llungstate4=squeeze(llung_medstate(4,:,:));
        rlungstate4=squeeze(rlung_medstate(4,:,:));
        thoraxstate4=squeeze(thorax_medstate(4,:,:));
        
        
        
        %%%%%fix orientations
        thorax=[thorax(:,1) -thorax(:,2)];
        thoraxstate2=[thoraxstate2(:,1) -thoraxstate2(:,2)];
        thoraxstate3=[thoraxstate3(:,1) -thoraxstate3(:,2)];
        thoraxstate4=[thoraxstate4(:,1) -thoraxstate4(:,2)];
        thoraxinf=[thoraxinf(:,1) -thoraxinf(:,2)];
        
        rlung=[rlung(:,1) -rlung(:,2)];
        rlungstate2=[rlungstate2(:,1) -rlungstate2(:,2)];
        rlungstate3=[rlungstate3(:,1) -rlungstate3(:,2)];
        rlungstate4=[rlungstate4(:,1) -rlungstate4(:,2)];
        rlunginf=[rlunginf(:,1) -rlunginf(:,2)];
        
        
        llung=[llung(:,1) -llung(:,2)];
        llungstate2=[llungstate2(:,1) -llungstate2(:,2)];
        llungstate3=[llungstate3(:,1) -llungstate3(:,2)];
        llungstate4=[llungstate4(:,1) -llungstate4(:,2)];
        llunginf=[llunginf(:,1) -llunginf(:,2)];
        
        newheart=[newheart(:,1) -newheart(:,2)];
        
        spondilus=[spondilus(:,1) -spondilus(:,2)];
        
        
        factor=148;   %%%%GN
        if strcmp(invref.fwd_model.name,'Grid model 36 x 36')==1
            factor=1.475;   %%%%%GREIT
        end
        
        thorax_scaled(:,:,1)=thoraxinf*factor;
        thorax_scaled(:,:,2)=thoraxstate4*factor;
        thorax_scaled(:,:,3)=thoraxstate3*factor;
        thorax_scaled(:,:,4)=thoraxstate2*factor; 
        thorax_scaled(:,:,5)=thorax*factor; 
        
        
                
        right_lung_scaled(:,:,1)=rlunginf*factor;  
        right_lung_scaled(:,:,2)=rlungstate4*factor;
        right_lung_scaled(:,:,3)=rlungstate3*factor;
        right_lung_scaled(:,:,4)=rlungstate2*factor;
        right_lung_scaled(:,:,5)=rlung*factor;
        
        left_lung_scaled(:,:,1)=llunginf*factor;  
        left_lung_scaled(:,:,2)=llungstate4*factor;
        left_lung_scaled(:,:,3)=llungstate3*factor;
        left_lung_scaled(:,:,4)=llungstate2*factor;
        left_lung_scaled(:,:,5)=llung*factor;
        
        spondilus_scaled(:,:,1)=spondilus*factor; 
        spondilus_scaled(:,:,2)=spondilus*factor; 
        spondilus_scaled(:,:,3)=spondilus*factor; 
        spondilus_scaled(:,:,4)=spondilus*factor; 
        spondilus_scaled(:,:,5)=spondilus*factor;
        
        heart_scaled(:,:,1)=newheart*factor;
        heart_scaled(:,:,2)=newheart*factor;
        heart_scaled(:,:,3)=newheart*factor;
        heart_scaled(:,:,4)=newheart*factor;
        heart_scaled(:,:,5)=newheart*factor;
               
        w_skin=zeros(size(elem_index,1),2);
        for frame=1:5
            
            %%%%right lung
            porg=polyshape([right_lung_scaled(:,1,frame) right_lung_scaled(:,2,frame)]);
            for element=1:size(elem_index,1)
                pelem=polyshape([polyelement(:,1,element) polyelement(:,2,element) ]);
                polyout(element) = intersect(pelem, porg);
                w_right_lung(element,frame)=area(polyout(element))/area(pelem);
            end
            
            %%%%left lung
            porg=polyshape([left_lung_scaled(:,1,frame) left_lung_scaled(:,2,frame)]);
            for element=1:size(elem_index,1)
                pelem=polyshape([polyelement(:,1,element) polyelement(:,2,element) ]);
                polyout(element) = intersect(pelem, porg);
                w_left_lung(element,frame)=area(polyout(element))/area(pelem);
            end
            
            %%%%spondilus
            porg=polyshape([spondilus_scaled(:,1,frame) spondilus_scaled(:,2,frame)]);
            for element=1:size(elem_index,1)
                pelem=polyshape([polyelement(:,1,element) polyelement(:,2,element) ]);
                polyout(element) = intersect(pelem, porg);
                w_spondilus(element,frame)=area(polyout(element))/area(pelem);
            end
            
            %%%%heart
            porg=polyshape([heart_scaled(:,1,frame) heart_scaled(:,2,frame)]);
            for element=1:size(elem_index,1)
                pelem=polyshape([polyelement(:,1,element) polyelement(:,2,element) ]);
                polyout(element) = intersect(pelem, porg);
                w_heart(element,frame)=area(polyout(element))/area(pelem);
            end
            
            %%%%skin
            porg=polyshape([thorax_scaled(:,1,frame) thorax_scaled(:,2,frame)]);
            for element=1:size(elem_index,1)
                pelem=polyshape([polyelement(:,1,element) polyelement(:,2,element) ]);
                polyout(element) = intersect(pelem, porg);
                if area(polyout(element))/area(pelem)<0.75&&area(polyout(element))/area(pelem)>0.25
                    w_skin(element,frame)=1;
                end
            end
        end
        w_muscle=1-(w_heart+w_skin+w_spondilus+w_right_lung+w_left_lung);
        
        s_lung=[Tissues_admittance(2,1) Lung_med_values(3,1)...
            Lung_med_values(2,1) Lung_med_values(1,1) Tissues_admittance(3,1)];
        s_heart=[Tissues_admittance(1,1) Tissues_admittance(1,1)...
            Tissues_admittance(1,1) Tissues_admittance(1,1) Tissues_admittance(1,1)];
        s_spondilus=[Tissues_admittance(5,1) Tissues_admittance(5,1)...
            Tissues_admittance(5,1) Tissues_admittance(5,1) Tissues_admittance(5,1)];
        s_muscles=[Tissues_admittance(4,1) Tissues_admittance(4,1)...
            Tissues_admittance(4,1) Tissues_admittance(4,1) Tissues_admittance(4,1)];
        s_skin=[Tissues_admittance(4,1) Tissues_admittance(4,1)...
            Tissues_admittance(4,1) Tissues_admittance(4,1) Tissues_admittance(4,1)];
        
        selem_ref=w_left_lung.*repmat(s_lung,size(elem_index,1),1)+...
            w_right_lung.*repmat(s_lung,size(elem_index,1),1)+...
            w_heart.*repmat(s_heart,size(elem_index,1),1)+...
            w_spondilus.*repmat(s_spondilus,size(elem_index,1),1)+...
            w_muscle.*repmat(s_muscles,size(elem_index,1),1)+...
            w_skin.*repmat(s_skin,size(elem_index,1),1);
        
        inv_reference(1)=invref;
        inv_reference(2)=invref;
        inv_reference(3)=invref;
        inv_reference(4)=invref;
        
        inv_reference(1).elem_data=selem_ref(:,4)-selem_ref(:,5);
        inv_reference(2).elem_data=selem_ref(:,3)-selem_ref(:,5);
        inv_reference(3).elem_data=selem_ref(:,2)-selem_ref(:,5);
        inv_reference(4).elem_data=selem_ref(:,1)-selem_ref(:,5);
        figure
        subplot(1,4,1)
        show_fem(inv_reference(1),1)
        subplot(1,4,2)
        show_fem(inv_reference(2),1)
        subplot(1,4,3)
        show_fem(inv_reference(3),1)
        subplot(1,4,4)
        show_fem(inv_reference(4),1)
        o=[];
%         invv1=invref; invv1.elem_data=selem_ref(:,1);
%         invv2=invref; invv2.elem_data=selem_ref(:,2);
        %%%%%done(?) : FLIP VERTICALLY EACH TISSUE AND BOUNDARY POLYGON!!!
        
        
        
        
        
    case{'Dynamic'}
        %addpath('D:\EIT\absolute_imaging\Models_FOM\Dynamic\Trauma\')
        %         load('Tissues.mat')
        addpath('D:\EIT\absolute_imaging\Models_FOM\Dynamic\Healthy\')
        
        load('Inputs_dat.mat')
        %         load('initial_total_positions.mat')
        %         load('initial_times.mat')
        %         load('electrodes.mat')
        %         load('Zelectrodes.mat')
        load('Electrodes_History.mat')
        load('History_of_Tissues.mat')
        factor=5.8;
        %factor=0.0578;
        figure
        for frame=1:30
            electrode_coords(:,1:2,frame)=mean(electrodes.positions(:,1:2,1,(frame-1)*16+1:frame*16),4);
            electrode_coords(:,3:4,frame)=mean(electrodes.positions(:,1:2,2,(frame-1)*16+1:frame*16),4);
            boundary_coords(:,1,frame)=reshape(electrode_coords(:,1:2,frame).',1,[]);
            boundary_coords(:,2,frame)=reshape(electrode_coords(:,3:4,frame).',1,[]);
            boundary_coords(:,:,frame)=boundary_coords(:,:,frame)*factor;
            %
            plot([boundary_coords(:,1,frame); boundary_coords(1,1,frame)],[boundary_coords(:,2,frame);boundary_coords(1,2,frame)],'m','LineWidth',2)
            hold on
            %
            %%%% Tissue_polygons: TISSUE_NODES X 2
            
            Left_Lung_coords(:,1,frame)=mean(History_of_Tissues.Left_lung(:,(frame-1)*16+1:2:2*frame*16),2);
            Left_Lung_coords(:,2,frame)=mean(History_of_Tissues.Left_lung(:,(frame-1)*16+2:2:2*frame*16),2);
            Left_Lung_coords(:,:,frame)=Left_Lung_coords(:,:,frame)*factor;
            
            Right_Lung_coords(:,1,frame)=mean(History_of_Tissues.Right_lung(:,(frame-1)*16+1:2:2*frame*16),2);
            Right_Lung_coords(:,2,frame)=mean(History_of_Tissues.Right_lung(:,(frame-1)*16+2:2:2*frame*16),2);
            Right_Lung_coords(:,:,frame)=Right_Lung_coords(:,:,frame)*factor;
            
            Heart_int_coords(:,1,frame)=mean(History_of_Tissues.Heart.interior(:,(frame-1)*16+1:2:2*frame*16),2);
            Heart_int_coords(:,2,frame)=mean(History_of_Tissues.Heart.interior(:,(frame-1)*16+2:2:2*frame*16),2);
            Heart_int_coords(:,:,frame)=Heart_int_coords(:,:,frame)*factor;
            
            Heart_ext_coords(:,1,frame)=mean(History_of_Tissues.Heart.exterior(:,(frame-1)*16+1:2:2*frame*16),2);
            Heart_ext_coords(:,2,frame)=mean(History_of_Tissues.Heart.exterior(:,(frame-1)*16+2:2:2*frame*16),2);
            Heart_ext_coords(:,:,frame)=Heart_ext_coords(:,:,frame)*factor;
            
            Bone1_coords(:,1,frame)=mean(History_of_Tissues.Bones.Bone1(:,(frame-1)*16+1:2:2*frame*16),2);
            Bone1_coords(:,2,frame)=mean(History_of_Tissues.Bones.Bone1(:,(frame-1)*16+2:2:2*frame*16),2);
            Bone1_coords(:,:,frame)=Bone1_coords(:,:,frame)*factor;
            
            Bone2_coords(:,1,frame)=mean(History_of_Tissues.Bones.Bone2(:,(frame-1)*16+1:2:2*frame*16),2);
            Bone2_coords(:,2,frame)=mean(History_of_Tissues.Bones.Bone2(:,(frame-1)*16+2:2:2*frame*16),2);
            Bone2_coords(:,:,frame)=Bone2_coords(:,:,frame)*factor;
            
            Bone3_coords(:,1,frame)=mean(History_of_Tissues.Bones.Bone3(:,(frame-1)*16+1:2:2*frame*16),2);
            Bone3_coords(:,2,frame)=mean(History_of_Tissues.Bones.Bone3(:,(frame-1)*16+2:2:2*frame*16),2);
            Bone3_coords(:,:,frame)=Bone3_coords(:,:,frame)*factor;
            
            Bone4_coords(:,1,frame)=mean(History_of_Tissues.Bones.Bone4(:,(frame-1)*16+1:2:2*frame*16),2);
            Bone4_coords(:,2,frame)=mean(History_of_Tissues.Bones.Bone4(:,(frame-1)*16+2:2:2*frame*16),2);
            Bone4_coords(:,:,frame)=Bone4_coords(:,:,frame)*factor;
            
            Bone5_coords(:,1,frame)=mean(History_of_Tissues.Bones.Bone5(:,(frame-1)*16+1:2:2*frame*16),2);
            Bone5_coords(:,2,frame)=mean(History_of_Tissues.Bones.Bone5(:,(frame-1)*16+2:2:2*frame*16),2);
            Bone5_coords(:,:,frame)=Bone5_coords(:,:,frame)*factor;
            
            Bone6_coords(:,1,frame)=mean(History_of_Tissues.Bones.Bone6(:,(frame-1)*16+1:2:2*frame*16),2);
            Bone6_coords(:,2,frame)=mean(History_of_Tissues.Bones.Bone6(:,(frame-1)*16+2:2:2*frame*16),2);
            Bone6_coords(:,:,frame)=Bone6_coords(:,:,frame)*factor;
            
            Bone7_coords(:,1,frame)=mean(History_of_Tissues.Bones.Bone7(:,(frame-1)*16+1:2:2*frame*16),2);
            Bone7_coords(:,2,frame)=mean(History_of_Tissues.Bones.Bone7(:,(frame-1)*16+2:2:2*frame*16),2);
            Bone7_coords(:,:,frame)=Bone7_coords(:,:,frame)*factor;
            
            Bone8_coords(:,1,frame)=mean(History_of_Tissues.Bones.Bone8(:,(frame-1)*16+1:2:2*frame*16),2);
            Bone8_coords(:,2,frame)=mean(History_of_Tissues.Bones.Bone8(:,(frame-1)*16+2:2:2*frame*16),2);
            Bone8_coords(:,:,frame)=Bone8_coords(:,:,frame)*factor;
            
            fat_coords(:,1,frame)=[Bone1_coords(3,1,frame);  Bone2_coords(1,1,frame); Bone2_coords(4,1,frame);...
                Bone3_coords(1,1,frame); Bone3_coords(2,1,frame); Bone3_coords(3,1,frame); Bone3_coords(4,1,frame);...
                Bone3_coords(5,1,frame); Bone3_coords(6,1,frame); Bone3_coords(7,1,frame); Bone3_coords(8,1,frame);...
                Bone3_coords(9,1,frame);...
                Bone4_coords(1,1,frame); Bone4_coords(4,1,frame);...
                Bone5_coords(2,1,frame); Bone5_coords(1,1,frame); Bone6_coords(1,1,frame); Bone6_coords(4,1,frame);...
                Bone7_coords(1,1,frame); Bone7_coords(4,1,frame); Bone8_coords(1,1,frame); Bone8_coords(4,1,frame);...
                Bone1_coords(4,1,frame)];
            fat_coords(:,2,frame)=[Bone1_coords(3,2,frame);  Bone2_coords(1,2,frame); Bone2_coords(4,2,frame);...
                Bone3_coords(1,2,frame); Bone3_coords(2,2,frame); Bone3_coords(3,2,frame); Bone3_coords(4,2,frame);...
                Bone3_coords(5,2,frame); Bone3_coords(6,2,frame); Bone3_coords(7,2,frame); Bone3_coords(8,2,frame);...
                Bone3_coords(9,2,frame);...
                Bone4_coords(1,2,frame); Bone4_coords(4,2,frame);...
                Bone5_coords(2,2,frame); Bone5_coords(1,2,frame); Bone6_coords(1,2,frame); Bone6_coords(4,2,frame);...
                Bone7_coords(1,2,frame); Bone7_coords(4,2,frame); Bone8_coords(1,2,frame); Bone8_coords(4,2,frame);...
                Bone1_coords(4,2,frame)];
            
            %             %%%get difference approximations
            %             if frame>1
            %             Left_Lung_diff_coords(:,:,frame-1)=(Left_Lung_coords(:,:,frame)+Left_Lung_coords(:,:,1))/2;
            %             Right_Lung_diff_coords(:,:,frame-1)=(Right_Lung_coords(:,:,frame)+Right_Lung_coords(:,:,1))/2;
            %             Bone1_diff_coords(:,:,frame-1)=(Bone1_coords(:,:,frame)+Bone1_coords(:,:,1))/2;
            %             Bone2_diff_coords(:,:,frame-1)=(Bone2_coords(:,:,frame)+Bone2_coords(:,:,1))/2;
            %             Bone3_diff_coords(:,:,frame-1)=(Bone3_coords(:,:,frame)+Bone3_coords(:,:,1))/2;
            %             Bone4_diff_coords(:,:,frame-1)=(Bone4_coords(:,:,frame)+Bone4_coords(:,:,1))/2;
            %             Bone5_diff_coords(:,:,frame-1)=(Bone5_coords(:,:,frame)+Bone5_coords(:,:,1))/2;
            %             Bone6_diff_coords(:,:,frame-1)=(Bone6_coords(:,:,frame)+Bone6_coords(:,:,1))/2;
            %             Bone7_diff_coords(:,:,frame-1)=(Bone7_coords(:,:,frame)+Bone7_coords(:,:,1))/2;
            %             Bone8_diff_coords(:,:,frame-1)=(Bone8_coords(:,:,frame)+Bone8_coords(:,:,1))/2;
            %             Heart_int_diff_coords(:,:,frame-1)=(Heart_int_coords(:,:,frame-1)+Heart_int_coords(:,:,frame-1))/2;
            %             Heart_ext_diff_coords(:,:,frame-1)=(Heart_ext_coords(:,:,frame-1)+Heart_ext_coords(:,:,frame-1))/2;
            %             fat_diff_coords(:,:,frame-1)=(fat_coords(:,:,frame-1)+fat_coords(:,:,frame-1))/2;
            %             end
            
            
            plot([Left_Lung_coords(:,1,frame); Left_Lung_coords(1,1,frame)],[Left_Lung_coords(:,2,frame);Left_Lung_coords(1,2,frame)],'b','LineWidth',2)
            hold on
            
            plot([Right_Lung_coords(:,1,frame); Right_Lung_coords(1,1,frame)],[Right_Lung_coords(:,2,frame);Right_Lung_coords(1,2,frame)],'b','LineWidth',2)
            hold on
            
            plot([Heart_int_coords(:,1,frame); Heart_int_coords(1,1,frame)],[Heart_int_coords(:,2,frame);Heart_int_coords(1,2,frame)],'r','LineWidth',2)
            hold on
            
            plot([Heart_ext_coords(:,1,frame); Heart_ext_coords(1,1,frame)],[Heart_ext_coords(:,2,frame);Heart_ext_coords(1,2,frame)],'r','LineWidth',2)
            hold on
            
            plot([Bone1_coords(:,1,frame); Bone1_coords(1,1,frame)],[Bone1_coords(:,2,frame); Bone1_coords(1,2,frame)],'k','LineWidth',2)
            hold on
            
            plot([Bone2_coords(:,1,frame); Bone2_coords(1,1,frame)],[Bone2_coords(:,2,frame); Bone2_coords(1,2,frame)],'k','LineWidth',2)
            hold on
            
            plot([Bone3_coords(:,1,frame); Bone3_coords(1,1,frame)],[Bone3_coords(:,2,frame); Bone3_coords(1,2,frame)],'k','LineWidth',2)
            hold on
            
            plot([Bone4_coords(:,1,frame); Bone4_coords(1,1,frame)],[Bone4_coords(:,2,frame); Bone4_coords(1,2,frame)],'k','LineWidth',2)
            hold on
            
            plot([Bone5_coords(:,1,frame); Bone5_coords(1,1,frame)],[Bone5_coords(:,2,frame); Bone5_coords(1,2,frame)],'k','LineWidth',2)
            hold on
            
            plot([Bone6_coords(:,1,frame); Bone6_coords(1,1,frame)],[Bone6_coords(:,2,frame); Bone6_coords(1,2,frame)],'k','LineWidth',2)
            hold on
            
            plot([Bone7_coords(:,1,frame); Bone7_coords(1,1,frame)],[Bone7_coords(:,2,frame); Bone7_coords(1,2,frame)],'k','LineWidth',2)
            hold on
            
            plot([Bone8_coords(:,1,frame); Bone8_coords(1,1,frame)],[Bone8_coords(:,2,frame); Bone8_coords(1,2,frame)],'k','LineWidth',2)
            hold on
            
            plot([fat_coords(:,1,frame); fat_coords(1,1,frame)],[fat_coords(:,2,frame); fat_coords(1,2,frame)],'c','LineWidth',2)
            hold on
            
        end
        
        if exist('weights.mat')==0
            
            %%%%% find percentage (w) of each element's area in each tissue's
            %%%%% ROI
            %%%heart interior
            for frame=1:30
                porg=polyshape([Heart_int_coords(:,1,frame) Heart_int_coords(:,2,frame)]);
                for element=1:size(elem_index,1)
                    pelem=polyshape([polyelement(:,1,element) polyelement(:,2,element) ]);
                    polyout(element) = intersect(pelem, porg);
                    w_heart_int(element,frame)=area(polyout(element))/area(pelem);
                end
                
                %%%heart exterior
                porg=polyshape([Heart_ext_coords(:,1,frame) Heart_ext_coords(:,2,frame)]);
                for element=1:size(elem_index,1)
                    pelem=polyshape([polyelement(:,1,element) polyelement(:,2,element) ]);
                    polyout(element) = intersect(pelem, porg);
                    w_heart_ext(element,frame)=area(polyout(element))/area(pelem)-w_heart_int(element);
                end
                
                %%%left lung
                porg=polyshape([Left_Lung_coords(:,1,frame) Left_Lung_coords(:,2,frame)]);
                for element=1:size(elem_index,1)
                    pelem=polyshape([polyelement(:,1,element) polyelement(:,2,element) ]);
                    polyout(element) = intersect(pelem, porg);
                    w_left_lung(element,frame)=area(polyout(element))/area(pelem);
                end
                
                %%%right lung
                porg=polyshape([Right_Lung_coords(:,1,frame) Right_Lung_coords(:,2,frame)]);
                for element=1:size(elem_index,1)
                    pelem=polyshape([polyelement(:,1,element) polyelement(:,2,element) ]);
                    polyout(element) = intersect(pelem, porg);
                    w_right_lung(element,frame)=area(polyout(element))/area(pelem);
                end
                
                %%%bone 1
                porg=polyshape([Bone1_coords(:,1,frame) Bone1_coords(:,2,frame)]);
                for element=1:size(elem_index,1)
                    pelem=polyshape([polyelement(:,1,element) polyelement(:,2,element) ]);
                    polyout(element) = intersect(pelem, porg);
                    w_bone1(element,frame)=area(polyout(element))/area(pelem);
                end
                
                %%%bone 2
                porg=polyshape([Bone2_coords(:,1,frame) Bone2_coords(:,2,frame)]);
                for element=1:size(elem_index,1)
                    pelem=polyshape([polyelement(:,1,element) polyelement(:,2,element) ]);
                    polyout(element) = intersect(pelem, porg);
                    w_bone2(element,frame)=area(polyout(element))/area(pelem);
                end
                
                %%%bone 3
                porg=polyshape([Bone3_coords(:,1,frame) Bone3_coords(:,2,frame)]);
                for element=1:size(elem_index,1)
                    pelem=polyshape([polyelement(:,1,element) polyelement(:,2,element) ]);
                    polyout(element) = intersect(pelem, porg);
                    w_bone3(element,frame)=area(polyout(element))/area(pelem);
                end
                
                %%%bone 4
                porg=polyshape([Bone4_coords(:,1,frame) Bone4_coords(:,2,frame)]);
                for element=1:size(elem_index,1)
                    pelem=polyshape([polyelement(:,1,element) polyelement(:,2,element) ]);
                    polyout(element) = intersect(pelem, porg);
                    w_bone4(element,frame)=area(polyout(element))/area(pelem);
                end
                
                %%%bone 5
                porg=polyshape([Bone5_coords(:,1,frame) Bone5_coords(:,2,frame)]);
                for element=1:size(elem_index,1)
                    pelem=polyshape([polyelement(:,1,element) polyelement(:,2,element) ]);
                    polyout(element) = intersect(pelem, porg);
                    w_bone5(element,frame)=area(polyout(element))/area(pelem);
                end
                
                %%%bone 6
                porg=polyshape([Bone6_coords(:,1,frame) Bone6_coords(:,2,frame)]);
                for element=1:size(elem_index,1)
                    pelem=polyshape([polyelement(:,1,element) polyelement(:,2,element) ]);
                    polyout(element) = intersect(pelem, porg);
                    w_bone6(element,frame)=area(polyout(element))/area(pelem);
                end
                
                %%%bone 7
                porg=polyshape([Bone7_coords(:,1,frame) Bone7_coords(:,2,frame)]);
                for element=1:size(elem_index,1)
                    pelem=polyshape([polyelement(:,1,element) polyelement(:,2,element) ]);
                    polyout(element) = intersect(pelem, porg);
                    w_bone7(element,frame)=area(polyout(element))/area(pelem);
                end
                
                %%%bone 8
                porg=polyshape([Bone8_coords(:,1,frame) Bone8_coords(:,2,frame)]);
                for element=1:size(elem_index,1)
                    pelem=polyshape([polyelement(:,1,element) polyelement(:,2,element) ]);
                    polyout(element) = intersect(pelem, porg);
                    w_bone8(element,frame)=area(polyout(element))/area(pelem);
                end
                
                %%%fat
                %for diff_frame=1:29
                porg=polyshape([fat_coords(:,1,frame) fat_coords(:,2,frame)]);
                for element=1:size(elem_index,1)
                    pelem=polyshape([polyelement(:,1,element) polyelement(:,2,element) ]);
                    polyout(element) = intersect(pelem, porg);
                    w_fat(element,frame)=(area(pelem)-area(polyout(element)))/area(pelem)...
                        -w_bone1(element,frame)-w_bone2(element,frame)...
                        -w_bone4(element,frame)-w_bone5(element,frame)-w_bone6(element,frame)...
                        -w_bone7(element,frame)-w_bone8(element,frame);
                end
                %end
                disp(frame)
            end
            
        end
        
        %%%% muscle-plasma weights
        wtissues=w_fat+w_bone1+w_bone2+w_bone3+w_bone4+w_bone5+w_bone6+w_bone7-...
            w_bone8+w_heart_ext+w_heart_int+w_left_lung+w_right_lung;
        w_muscle=1-wtissues;
        
        
        s_heart_ext=0.186;
        for frame=1:30
            s_heart_int(frame)=mean(Inputs_dat.Graph.heart.conductivity((frame-1)*16+1:1:frame*16),2);
            s_left_lung(frame)=mean(Inputs_dat.Graph.left_lung.conductivity((frame-1)*16+1:1:frame*16),2);
            s_right_lung(frame)=mean(Inputs_dat.Graph.left_lung.conductivity((frame-1)*16+1:1:frame*16),2);
            s_muscles(frame)=mean(Inputs_dat.Graph.muscles.conductivity((frame-1)*16+1:1:frame*16),2);
        end
        s_bones=0.07;
        s_fat=0.03;
        
        selem_ref=w_heart_int.*repmat(s_heart_int,size(elem_index,1),1)...
            +w_heart_ext.*repmat(s_heart_ext,size(elem_index,1),1)...
            +w_left_lung.*repmat(s_left_lung,size(elem_index,1),1)...
            +w_right_lung.*repmat(s_right_lung,size(elem_index,1),1)...
            +w_fat.*repmat(s_fat,size(elem_index,1),1)...
            +w_bone1.*repmat(s_bones,size(elem_index,1),1)...
            +w_bone2.*repmat(s_bones,size(elem_index,1),1)...
            +w_bone3.*repmat(s_bones,size(elem_index,1),1)...
            +w_bone4.*repmat(s_bones,size(elem_index,1),1)...
            +w_bone5.*repmat(s_bones,size(elem_index,1),1)...
            +w_bone6.*repmat(s_bones,size(elem_index,1),1)...
            +w_bone7.*repmat(s_bones,size(elem_index,1),1)...
            +w_bone8.*repmat(s_bones,size(elem_index,1),1)...
            +w_muscle.*repmat(s_muscles,size(elem_index,1),1);
        
        for diff_frame=1:29
            inv_reference(diff_frame)=invref;
            selem_ref_diff(:,diff_frame)=selem_ref(:,diff_frame+1)-selem_ref(:,1);
            inv_reference(diff_frame).elem_data=selem_ref_diff(:,diff_frame);
        end
        
        for frame=1:30
            inv_reference_stat(frame)=invref;
            inv_reference_stat(frame).elem_data=selem_ref(:,frame);
        end
        
        show_slices(inv_reference_stat)
        
        
end



end
