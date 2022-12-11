path='C:\Users\Chris\Desktop\PHD_files\eidors-v3.9-ng\eidors\EIT_Circuit_Sim\';
cavity=5;
addpath(['D:\EIT\absolute_imaging\Thoracic_Cavities\Thorax' num2str(cavity) '\'])
addpath(['D:\EIT\absolute_imaging\'])

addpath(path)
run C:\Users\Chris\Desktop\PHD_files\eidors-v3.9-ng\eidors\startup.m

load(['llung' num2str(cavity) '.mat'])
load(['llung' num2str(cavity) 'inf.mat'])
load(['newheart' num2str(cavity) '.mat'])
load(['rlung' num2str(cavity) '.mat'])
load(['rlung' num2str(cavity) 'inf.mat'])
load(['spondilus' num2str(cavity) '.mat'])
load(['thorax' num2str(cavity) '.mat'])
load(['thorax' num2str(cavity) 'inf.mat'])

llungdef=eval(['llung' num2str(cavity)]);
rlungdef=eval(['rlung' num2str(cavity)]);
thoraxdef=eval(['thorax' num2str(cavity)]);
spondilus=eval(['spondilus' num2str(cavity)]);
newheart=eval(['newheart' num2str(cavity)]);

load([path 'Structures\elec_pos_dev1.mat'])
freq=100000;
conductivity_values='variational';
load(['Thorax3Ddef_' num2str(cavity) '_Assym.mat'])
img{1}=thorax_model;

for state=2:4
    if state==4&&cavity==2
        thorax_medstate(4,:,:)=thoraxdef+(3.6-1)*(thoraxinf-thoraxdef)/4;
        llung_medstate(4,:,:)=llungdef+(3.6-1)*(llunginf-llungdef)/4;
        rlung_medstate(4,:,:)=rlungdef+(3.6-1)*(rlunginf-rlungdef)/4;
    elseif state==2&&cavity==4
        thorax_medstate(2,:,:)=thoraxdef+(2.6-1)*(thoraxinf-thoraxdef)/4;
        llung_medstate(2,:,:)=llungdef+(2.6-1)*(llunginf-llungdef)/4;
        rlung_medstate(2,:,:)=rlungdef+(2.6-1)*(rlunginf-rlungdef)/4;
    elseif state==2&&cavity==5
        thorax_medstate(2,:,:)=thoraxdef+(1.6-1)*(thoraxinf-thoraxdef)/4;
        llung_medstate(2,:,:)=llungdef+(1.6-1)*(llunginf-llungdef)/4;
        rlung_medstate(2,:,:)=rlungdef+(1.6-1)*(rlunginf-rlungdef)/4;
        
    else
        llung_medstate(state,:,:)=llungdef+(state-1)*(llunginf-llungdef)/4;
        rlung_medstate(state,:,:)=rlungdef+(state-1)*(rlunginf-rlungdef)/4;
        thorax_medstate(state,:,:)=thoraxdef+(state-1)*(thoraxinf-thoraxdef)/4;
    end
    %     if cavity==1
    %         %%%%%TO CHECK: WHY CAVITY 1 DOES NOT WORK EVEN IN DEFLATED CASE?
    %         thorax_medstate(state,:,:)=thoraxdef;
    %         thorax_medstate(state,:,:)=thoraxdef;
    %         rlung_medstate(state,:,:)=rlungdef;
    %         llung_medstate(state,:,:)=llungdef;
    %     end
    if cavity==2
        objects={squeeze(thorax_medstate(state,:,:)), squeeze(rlung_medstate(state,:,:)),...
            squeeze(llung_medstate(state,:,:)), spondilus, newheart};
    elseif cavity==5
        objects={squeeze(thorax_medstate(state,:,:)), squeeze(rlung_medstate(state,:,:)),...
        squeeze(llung_medstate(state,:,:))};
    
    else
        objects={squeeze(thorax_medstate(state,:,:)), squeeze(rlung_medstate(state,:,:)),...
            squeeze(llung_medstate(state,:,:)), newheart};
    end
    
    [imgs,elec_pos]=create_Forw_model(16, objects,'Asymmetric',0.05,elec_pos);
    [imgs,elems_data,Tissues_admittance_mod]=set_elem_conducts_v2(imgs,state,freq,conductivity_values,path);
    img{state}=imgs;
    
    if state>1&&state<5
        Lung_med_values(state-1,:)=Tissues_admittance_mod(2,:);
    end
    
end
load(['Thorax3Dinf_' num2str(cavity) '_Assym.mat'])
img{5}=thorax_model;


% figure
% plot(abs(img.elem_data),'*')
% hold on
thorax_model=img;

