function thorax_model=Set_Thorax_Structure(state,freq,symmetry,Rel,conductivity_values,Device)
clc
%thorax_model=Set_Thorax_Structure('Deflated',100000,'Asymmetric',0.05,'variational','Desktop');
if strcmp(Device,'Laptop')
    path='C:\Users\user\Desktop\eidors-v3.9-ng\eidors\EIT_Circuit_Sim\';
else
    path='C:\Users\Chris\Desktop\PHD_files\eidors-v3.9-ng\eidors\EIT_Circuit_Sim\';
end
% switch Rel
%     case 0.05
%         electrode_size='Large_AgCl';
%     otherwise
%         electrode_size='Small_AgCl';
% end
%%%%set geometry
%initial_struct=['Thorax_' state '_' symmetry '_Rel=' num2str(Rel)];
initial_struct=['Thorax2_' state '_' symmetry '_Rel=' num2str(Rel)];
structure_path=[path 'Structures\Thoracic_New\'];
structure_geometry=[structure_path initial_struct '.mat'];
if exist(structure_geometry,'file')~=2
    [img,elec_pos]=create_thorax_geometry(structure_geometry,state,symmetry,Rel,path);
else
    load(structure_geometry)  
end
%%%%set element conductivities
thorax_model=set_elem_conducts(img,state,freq,conductivity_values,path);
%img.fwd_model.nodes=12*img.fwd_model.nodes;

%%%%create the equivalent N-port circuit .cir file
% name=['Thorax2020_' state '_' symmetry '_' num2str(round(freq/1000)) 'k'];
% element_path=[structure_path 'freq_elements\'];
% element_data=thorax_model.elem_data;
% elem_file_name=[element_path name '_Rel=' num2str(Rel) '.mat'];
% save(elem_file_name,'element_data')
% eit_spiceRLC(thorax_model,name,freq,'Thorax',Device,electrode_size);
end

function [img,elec_pos]=create_thorax_geometry(structure_geometry,state,symmetry,Rel,path)

thorax = shape_library('get','adult_male','boundary');
rlung  = shape_library('get','adult_male','right_lung');
llung  = shape_library('get','adult_male','left_lung');
%%%%DEFAULT MODEL IS CONSIDERATED AS DEFLATED CASE
%WE BUILD THE INFLATED CASE
thoraxinf=thorax;
%thoraxinf(1:8,2)=thoraxinf(1:8,2)+0.08;
thoraxinf(1:8,2)=thoraxinf(1:8,2)+thoraxinf(1:8,2)*0.08;
%thoraxinf(end-7:end,2)=thoraxinf(end-7:end,2)+0.08;
thoraxinf(end-7:end,2)=thoraxinf(end-7:end,2)+thoraxinf(end-7:end,2)*0.08;
rlunginf=rlung;
rlunginf(end-6:end,2)=rlunginf(end-6:end,2)+0.15;
llunginf=llung;
llunginf(2:6,2)=llunginf(2:6,2)+0.14;
% one could also run:
% shape_library('get','adult_male');
% to get all the info at once in a struct
load([path '\Boundaries\Spondilus1.mat'])
load([path '\Boundaries\NewHeart.mat'])
% show the library image
shape_library('show','adult_male');
hold on
plot([spondilus(:,1); spondilus(1,1)],[spondilus(:,2); spondilus(1,2)],'-o','LineWidth',2)
hold on
plot([newheart(:,1); newheart(1,1)], [newheart(:,2); newheart(1,2)],'-o','LineWidth',2)
figure
switch state
    case 'Deflated'
        o=[];
    case 'Inflated'
        thorax=thoraxinf;
        llung=llunginf;
        rlung=rlunginf;
end

objects={thorax, rlung, llung, flipud(spondilus), newheart};

%objects={thorax2, rlung2, llung2, spondilus2, newheart2};
%objects={thorax3, rlung3, llung3, spondilus3, newheart3};
%objects={thoraxinf, rlunginf, llunginf, spondilus2, newheart2};
%objects={thoraxinf, rlunginf, llunginf, spondilus3, newheart3};
%objects={thoraxinf, rlunginf, llunginf, spondilus4, newheart4};

objects={thorax5, rlung5, llung5};

objects={thoraxinf, rlunginf, llunginf};

if strcmp(symmetry,'Asymmetric')
    load([path 'Structures\elec_pos_dev1.mat'])
    [img,elec_pos]=create_Forw_model(16, objects,symmetry,Rel,elec_pos);
else
    [img,elec_pos]=create_Forw_model(16, objects,symmetry,Rel);
end

save(structure_geometry,'img')
fprintf('saved new structure geometry: %s\n',structure_geometry)

end

function thorax_model=set_elem_conducts(img,state,freq,conductivity_values,path)

%clear all previous values
img.elem_data(:)=1;
%select tissuelist depended on the state
switch state
    case 'Deflated'
        listofinterest=[1 3 4 5 8, 10];
    case 'Inflated'
        listofinterest=[1 2 4 5 8, 10];
end

%load Tissues_admittance!
checktissuefile=[path 'Tissue_Properties\Thorax\Tissues_',num2str(round(freq/1000)),'kHz.mat'];
while exist(checktissuefile,'file')~=2
    fprintf('RUN Python script for f=%4.0f Hz!\n',freq)
    fprintf('Then press enter\n')
    pause()
end
load(checktissuefile)

%set admittances to elements
img.elem_data(:)=Tissues_admittance(listofinterest(3),1);
img.elem_data(img.fwd_model.mat_idx{2})= ...
    Tissues_admittance(listofinterest(2),1)+...
    1i*Tissues_admittance(listofinterest(2),2); % rlung
img.elem_data(img.fwd_model.mat_idx{3})= ...
    Tissues_admittance(listofinterest(2),1)+...
    1i*Tissues_admittance(listofinterest(2),2); % llung
img.elem_data(img.fwd_model.mat_idx{4})=...
    Tissues_admittance(listofinterest(4),1)+...
    1i*Tissues_admittance(listofinterest(4),2); % spondilus
img.elem_data(img.fwd_model.mat_idx{5})= ...
    Tissues_admittance(listofinterest(1),1)+...
    1i*Tissues_admittance(listofinterest(1),2); % heart
% skin 
[srf, idx] = find_boundary(img.fwd_model.elems);
bound_nodes_coords1=img.fwd_model.nodes(srf(:,1),:);
bound_nodes_coords2=img.fwd_model.nodes(srf(:,2),:);
bound_nodes_coords3=img.fwd_model.nodes(srf(:,3),:);
bound_element_coords(:,1)=(bound_nodes_coords1(:,1)+bound_nodes_coords2(:,1)+bound_nodes_coords3(:,1))/3;
bound_element_coords(:,2)=(bound_nodes_coords1(:,2)+bound_nodes_coords2(:,2)+bound_nodes_coords3(:,2))/3;
bound_element_coords(:,3)=(bound_nodes_coords1(:,3)+bound_nodes_coords2(:,3)+bound_nodes_coords3(:,3))/3;

bound_peripheral_local_ind=find(bound_element_coords(:,3)<=0.98&bound_element_coords(:,3)>=0.02);
bound_peripheral_global_ind=idx(bound_peripheral_local_ind);
img.elem_data(bound_peripheral_global_ind)=(Tissues_admittance(listofinterest(6),1)+...
    1i*Tissues_admittance(listofinterest(6),2)+Tissues_admittance(listofinterest(5),1)+...
    1i*Tissues_admittance(listofinterest(5),2))/2;
figure
show_fem(img); view(0,70);

switch conductivity_values
        case 'constant'
            %do nothing!
        case 'variational'
            %muscle
           L1=length(img.elem_data(img.fwd_model.mat_idx{1}));
           img.elem_data(img.fwd_model.mat_idx{1})= ...
               img.elem_data(img.fwd_model.mat_idx{1}).*(1+randn(L1,1)*0.02);
           %right lung
           L2=length(img.elem_data(img.fwd_model.mat_idx{2}));
           img.elem_data(img.fwd_model.mat_idx{2})= ...
               img.elem_data(img.fwd_model.mat_idx{2}).*(1+randn(L2,1)*0.03);
           %left lung
           L3=length(img.elem_data(img.fwd_model.mat_idx{3}));
           img.elem_data(img.fwd_model.mat_idx{3})= ...
               img.elem_data(img.fwd_model.mat_idx{3}).*(1+randn(L3,1)*0.03);
           %spondilus
           L4=length(img.elem_data(img.fwd_model.mat_idx{4}));
           img.elem_data(img.fwd_model.mat_idx{4})= ...
               img.elem_data(img.fwd_model.mat_idx{4}).*(1+randn(L4,1)*0.01);
           %heart
           L5=length(img.elem_data(img.fwd_model.mat_idx{5}));
           img.elem_data(img.fwd_model.mat_idx{5})= ...
               img.elem_data(img.fwd_model.mat_idx{5}).*(1+randn(L5,1)*0.02);
end
figure
plot(abs(img.elem_data),'*')
hold on
thorax_model=img;
end