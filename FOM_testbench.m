%%%%%%  Calculate FIGURES OF MERIT---FOMS for quantitative results and
%%%%%%  comparison
%%%%%% Uses Obtained Results in folder Results\
%%%%% Only a single thoracic cavity (for 5 states) supported in this demo.
%%%%% For more, please contact chdim@central.ntua.gr

%%%%%%%path to EIDORS
Eidorspath='C:\Users\Chris\Desktop\PHD_files\'; %%%%%Set path to eidors!!
startuped=0;
if ~exist('startuped')||startuped==0
    try
        run ([Eidorspath,'eidors-v3.9-ng\eidors\startup.m'])
        startuped=1;
    catch
        error('This code needs the EIDORS package to run properly, please download from https://eidors3d.sourceforge.net/download.shtml')
    end
end
%%%%%% Number of Electrodes (current v. supports only 16)
N=16;
skipcurr=0;
skipvolt=0;
cavity=1;
if cavity>1
    %%%%%% works in original code, memory issues in Github!
    error('This demo works only for the 1st cavity! Please select cavity=1')
end
%%%%scaling factor, do not change!!
GN_factor=148; 
%%% load model(s)
%sim_model=load_sim_model(type);
%%%% obtain references
path='Models_FOM\Thorax_3D\';
addpath(path)
addpath('Models_FOM\Calculate_Foms\')
addpath('Results\')
addpath('Thoracic_Cavities\Thorax1\')
load(['invreference' num2str(cavity) '_inf.mat'])
%load(['invreferenceGREIT' num2str(cavity) '_inf.mat'])
load(['invreference' num2str(cavity) '_MoM.mat'])
%%%% obtain images
%%%%
load(['diff_smooth' num2str(cavity) '.mat'])
load(['diff_TV' num2str(cavity) '.mat'])
%load(['diff_GREIT' num2str(cavity) '.mat'])
load(['diff_move' num2str(cavity) '.mat'])
load(['abs_smooth' num2str(cavity) '.mat'])
load(['Non_linear_All' num2str(cavity) '.mat'])
load(['Non_linear_ROI_Tissues' num2str(cavity) '.mat'])
%load(['Non_linear_ROI_Central' num2str(cavity) '.mat'])
load(['images_MoM' num2str(cavity) '.mat'])
load(['images_MoM_SBL' num2str(cavity) '.mat'])
%%%load ROI elements
%load('Elems_ROI_Central.mat')
%load('Elems_ROI_Tissues.mat')

%%%%load lung curves

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
%%%scale and flip


%%%for GN
rlung=rlung*GN_factor;
rlungstate2=rlungstate2*GN_factor;
rlungstate3=rlungstate3*GN_factor;
rlungstate4=rlungstate4*GN_factor;
rlunginf=rlunginf*GN_factor;

llung=llung*GN_factor;
llungstate2=llungstate2*GN_factor;
llungstate3=llungstate3*GN_factor;
llungstate4=llungstate4*GN_factor;
llunginf=llunginf*GN_factor;

RLUNGGN(1,:,:)=rlung;
RLUNGGN(2,:,:)=rlungstate2;
RLUNGGN(3,:,:)=rlungstate3;
RLUNGGN(4,:,:)=rlungstate4;
RLUNGGN(5,:,:)=rlunginf;

LLUNGGN(1,:,:)=llung;
LLUNGGN(2,:,:)=llungstate2;
LLUNGGN(3,:,:)=llungstate3;
LLUNGGN(4,:,:)=llungstate4;
LLUNGGN(5,:,:)=llunginf;



%%%find centroids
%%%%for GN
[Xr(1),Yr(1)]=centroid(polyshape(rlungstate2(:,1), rlungstate2(:,2)));
[Xl(1),Yl(1)]=centroid(polyshape(llungstate2(:,1), llungstate2(:,2)));
[Xr(2),Yr(2)]=centroid(polyshape(rlungstate3(:,1), rlungstate3(:,2)));
[Xl(2),Yl(2)]=centroid(polyshape(llungstate3(:,1), llungstate3(:,2)));
[Xr(3),Yr(3)]=centroid(polyshape(rlungstate4(:,1), rlungstate4(:,2)));
[Xl(3),Yl(3)]=centroid(polyshape(llungstate4(:,1), llungstate4(:,2)));
[Xr(4),Yr(4)]=centroid(polyshape(rlunginf(:,1), rlunginf(:,2)));
[Xl(4),Yl(4)]=centroid(polyshape(llunginf(:,1), llunginf(:,2)));


%%%load element-to-node-map
load('polyelement.mat')
load('polyelement_MoM.mat')
%images_diff_smooth=flip(images_diff_smooth);
%%%%%% smooth linear diff
%%%target amplitude
for imm=1:4
    [TA_diff_smooth(imm),imagesmooth,imagesmoothfilt]=...
        Calculate_Target_Thorax_Amplitude(images_diff_smooth(imm),polyelement);
    TA_diff_smooth(imm)=TA_diff_smooth(imm)/max(abs(images_diff_smooth(imm).elem_data));
    %%%%%extract reconstructed lung shapes
    images_diff_smooth_filtered{imm}=imagesmoothfilt;
    [llcand_diff_smooth{imm},rlcand_diff_smooth{imm},centrelements]=extract_reconstructed_shapes...
        (images_diff_smooth_filtered{imm},polyelement);
end
%%%%Position Error
for imm=1:4
    [PE_diff_smooth{imm},RL_rec_centre,LL_rec_centre]=Calculate_PE(polyelement,llcand_diff_smooth{imm},...
        rlcand_diff_smooth{imm},[Xl(imm),Yl(imm)], [Xr(imm),Yr(imm)],centrelements);
    %%%%RESOLUTION
    [RES_diff_smooth{imm},llung_diff_smooth_rec_area(imm),rlung_diff_smooth_rec_area(imm)]=...
        Calculate_RES(polyelement,llcand_diff_smooth{imm},rlcand_diff_smooth{imm});
    %%%%Shape Deformation (for each lung particularly and overall SD)
    SD_diff_smooth{imm}=Calculate_SD(polyelement,llcand_diff_smooth{imm},rlcand_diff_smooth{imm},...
        squeeze(LLUNGGN(imm+1,:,:)),squeeze(RLUNGGN(imm+1,:,:)),llung_diff_smooth_rec_area(imm),rlung_diff_smooth_rec_area(imm));
    RNG_diff_smooth(imm)=Calculate_RNG(polyelement,llcand_diff_smooth{imm},rlcand_diff_smooth{imm},images_diff_smooth(imm));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
imageframe=1;
subplot(1,2,1)
show_fem(images_diff_smooth(imageframe))
hold on
plot([squeeze(RLUNGGN(imageframe,:,1))'; squeeze(RLUNGGN(imageframe,1,1))],...
    [squeeze(RLUNGGN(imageframe,:,2))'; squeeze(RLUNGGN(imageframe,1,2))],'r','LineWidth',2)
hold on
plot([squeeze(LLUNGGN(imageframe,:,1))'; squeeze(LLUNGGN(imageframe,1,1))],...
    [squeeze(LLUNGGN(imageframe,:,2))'; squeeze(LLUNGGN(imageframe,1,2))],'r','LineWidth',2)
axis off
subplot(1,2,2)
show_fem(images_diff_smooth_filtered{imageframe})
hold on
plot([squeeze(RLUNGGN(imageframe,:,1))'; squeeze(RLUNGGN(imageframe,1,1))],...
    [squeeze(RLUNGGN(imageframe,:,2))'; squeeze(RLUNGGN(imageframe,1,2))],'c','LineWidth',2)
hold on
plot([squeeze(LLUNGGN(imageframe,:,1))'; squeeze(LLUNGGN(imageframe,1,1))],...
    [squeeze(LLUNGGN(imageframe,:,2))'; squeeze(LLUNGGN(imageframe,1,2))],'c','LineWidth',2)
hold on
plot([0 0],[-120 120],'m--','LineWidth',3)
hold on
plot(Xl(1),Yl(1),'-p','MarkerSize',15,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor','c')
hold on
plot(Xr(1),Yr(1),'-p','MarkerSize',15,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor','c')
hold on
plot(RL_rec_centre(1),RL_rec_centre(2),'-o','MarkerSize',12,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','c')
hold on
plot(LL_rec_centre(1),LL_rec_centre(2),'-o','MarkerSize',12,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','c')
axis off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%======================================================
%%%%%% TV linear diff
for imm=1:4
    %%%target amplitude
    [TA_diff_TV(imm),imagesTV,imagesTVfiltered]=...
        Calculate_Target_Thorax_Amplitude(images_diff_TV(imm),polyelement);
    TA_diff_TV(imm)=TA_diff_TV(imm)/max(abs(images_diff_TV(imm).elem_data));
    images_diff_TV_filtered{imm}=imagesmoothfilt;
    %%%%%extract reconstructed lung shapes
    [llcand_diff_TV{imm},rlcand_diff_TV{imm},centrelements]=extract_reconstructed_shapes...
        (images_diff_TV_filtered{imm},polyelement);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%Position Error
    PE_diff_TV{imm}=Calculate_PE(polyelement,llcand_diff_TV{imm},rlcand_diff_TV{imm},...
        [Xl(imm),Yl(imm)], [Xr(imm),Yr(imm)],centrelements);
    %%%%RESOLUTION
    [RES_diff_TV{imm},llung_diff_TV_rec_area(imm),rlung_diff_TV_rec_area(imm)]=...
        Calculate_RES(polyelement,llcand_diff_TV{imm},rlcand_diff_TV{imm});
    %%%%Shape Deformation (for each lung particularly and overall SD)
    SD_diff_TV{imm}=Calculate_SD(polyelement,llcand_diff_TV{imm},rlcand_diff_TV{imm},...
        squeeze(LLUNGGN(imm+1,:,:)),squeeze(RLUNGGN(imm+1,:,:)),llung_diff_TV_rec_area(imm),rlung_diff_TV_rec_area(imm));
    RNG_diff_TV(imm)=Calculate_RNG(polyelement,llcand_diff_TV{imm},rlcand_diff_TV{imm},images_diff_TV(imm));
    
end

%%%%======================================================
%%%%%% Absolute imaging differences
for imm=1:4
    %%%target amplitude
    [TA_ABS(imm),imagABS,imagABSfiltered]=...
        Calculate_Target_Thorax_Amplitude(imagesABS(imm),polyelement);
    TA_ABS(imm)=TA_ABS(imm)/max(abs(imagesABS(imm).elem_data));
    %%%%%extract reconstructed lung shapes
    imagesABS_filtered{imm}=imagABSfiltered;
    [llcand_ABS{imm},rlcand_ABS{imm},centrelements]=extract_reconstructed_shapes...
        (imagesABS_filtered{imm},polyelement);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%Position Error
    PE_ABS{imm}=Calculate_PE(polyelement,llcand_ABS{imm},rlcand_ABS{imm},...
        [Xl(imm),Yl(imm)], [Xr(imm),Yr(imm)],centrelements);
    %%%%RESOLUTION
    [RES_ABS{imm},llung_ABS_rec_area(imm),rlung_ABS_rec_area(imm)]=...
        Calculate_RES(polyelement,llcand_ABS{imm},rlcand_ABS{imm});
    %%%%Shape Deformation (for each lung particularly and overall SD)
    SD_ABS{imm}=Calculate_SD(polyelement,llcand_ABS{imm},rlcand_ABS{imm},...
        squeeze(LLUNGGN(imm+1,:,:)),squeeze(RLUNGGN(imm+1,:,:)),llung_ABS_rec_area(imm),rlung_ABS_rec_area(imm));
    RNG_ABS(imm)=Calculate_RNG(polyelement,llcand_ABS{imm},rlcand_ABS{imm},imagesABS(imm));
    
end


%%%%=======================================
%%%%%%Linear GN imaging with Movement Prior
images_diff_movep=images_diff_move;
for imm=1:4
    images_diff_move(imm).elem_data=images_diff_move(imm).elem_data(1:length(imagesABS(imm).elem_data));
    %%%target amplitude
    [TA_movement(imm),imgmove,imagmove_filtered]=...
        Calculate_Target_Thorax_Amplitude(images_diff_move(imm),polyelement);
    TA_movement(imm)=TA_movement(imm)/max(abs(images_diff_move(imm).elem_data));
    %%%%%extract reconstructed lung shapes
    images_movement_filtered{imm}=imagmove_filtered;
    [llcand_movement{imm},rlcand_movement{imm},centrelements]=extract_reconstructed_shapes...
        (images_movement_filtered{imm},polyelement);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%Position Error
    PE_movement{imm}=Calculate_PE(polyelement,llcand_movement{imm},rlcand_movement{imm},...
        [Xl(imm),Yl(imm)], [Xr(imm),Yr(imm)],centrelements);
    %%%%RESOLUTION
    [RES_movement{imm},llung_movement_rec_area(imm),rlung_movement_rec_area(imm)]=...
        Calculate_RES(polyelement,llcand_movement{imm},rlcand_movement{imm});
    %%%%Shape Deformation (for each lung particularly and overall SD)
    SD_movement{imm}=Calculate_SD(polyelement,llcand_movement{imm},rlcand_movement{imm},...
        squeeze(LLUNGGN(imm+1,:,:)),squeeze(RLUNGGN(imm+1,:,:)),llung_movement_rec_area(imm),rlung_movement_rec_area(imm));
    RNG_movement(imm)=Calculate_RNG(polyelement,llcand_movement{imm},rlcand_movement{imm},images_diff_move(imm));
end

%%%========================================================================
%%%%%%Non-linear All_ROI imaging
for imm=1:4
    %%%target amplitude
    [TA_Non_Linear_All(imm),imgAll,imagAllfiltered]=...
        Calculate_Target_Thorax_Amplitude(imagesAll(imm),polyelement);
    TA_Non_Linear_All(imm)=TA_Non_Linear_All(imm)/max(abs(imagesAll(imm).elem_data));
    images_Non_Linear_All_filtered{imm}=imagAllfiltered;
    %%%%%extract reconstructed lung shapes
    [llcand_Non_Linear_All{imm},rlcand_Non_Linear_All{imm},centrelements]=extract_reconstructed_shapes...
        (images_Non_Linear_All_filtered{imm},polyelement);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%Position Error
    PE_Non_Linear_All{imm}=Calculate_PE(polyelement,llcand_Non_Linear_All{imm},rlcand_Non_Linear_All{imm},...
        [Xl(imm),Yl(imm)], [Xr(imm),Yr(imm)],centrelements);
    %%%%RESOLUTION
    [RES_Non_Linear_All{imm},llung_Non_Linear_All_rec_area(imm),rlung_Non_Linear_All_rec_area(imm)]=...
        Calculate_RES(polyelement,llcand_Non_Linear_All{imm},rlcand_Non_Linear_All{imm});
    %%%%Shape Deformation (for each lung particularly and overall SD)
    SD_Non_Linear_All{imm}=Calculate_SD(polyelement,llcand_Non_Linear_All{imm},rlcand_Non_Linear_All{imm},...
        squeeze(LLUNGGN(imm+1,:,:)),squeeze(RLUNGGN(imm+1,:,:)),llung_Non_Linear_All_rec_area(imm),rlung_Non_Linear_All_rec_area(imm));
    RNG_Non_Linear_All(imm)=Calculate_RNG(polyelement,llcand_Non_Linear_All{imm},rlcand_Non_Linear_All{imm},imagesAll(imm));
end




%%%========================================================================
%%%%%% MoM imaging_Laplace
for imm=1:4
    images_MoM_structed(imm).elem_data=imagesMoM(imm,:);
    %%%target amplitude
    [TA_MoM(imm),imgMoM,images_MoM_filtered]=...
        Calculate_Target_Thorax_Amplitude(images_MoM_structed(imm),polyelement_MoM);
    TA_MoM(imm)=TA_MoM(imm)/max(abs(images_MoM_structed(imm).elem_data-1));
    %%%%%extract reconstructed lung shapes
    images_MoM_filtered_seq{imm}=images_MoM_filtered;
    [llcand_MoM{imm},rlcand_MoM{imm},centrelements_MoM]=extract_reconstructed_shapes...
        (images_MoM_filtered_seq{imm},polyelement_MoM);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%Position Error
    PE_MoM{imm}=Calculate_PE(polyelement_MoM,llcand_MoM{imm},rlcand_MoM{imm},...
        [Xl(imm),Yl(imm)], [Xr(imm),Yr(imm)],centrelements_MoM);
    %%%%RESOLUTION
    [RES_MoM{imm},llung_MoM_rec_area(imm),rlung_MoM_rec_area(imm)]=...
        Calculate_RES(polyelement_MoM,llcand_MoM{imm},rlcand_MoM{imm});
    %%%%Shape Deformation (for each lung particularly and overall SD)
    SD_MoM{imm}=Calculate_SD(polyelement_MoM,llcand_MoM{imm},rlcand_MoM{imm},...
        squeeze(LLUNGGN(imm+1,:,:)),squeeze(RLUNGGN(imm+1,:,:)),llung_MoM_rec_area(imm),rlung_MoM_rec_area(imm));
    RNG_MoM(imm)=Calculate_RNG(polyelement_MoM,llcand_MoM{imm},rlcand_MoM{imm},images_MoM_structed(imm));
end


%%%========================================================================
%%%%%% MoM imaging_SBL
for imm=1:4
    images_MoM_SBL_structed(imm).elem_data=images_MoM_SBL(imm,:);
    %%%target amplitude
    [TA_MoM_SBL(imm),imgMoMSBL,images_MoM_SBL_filtered]=...
        Calculate_Target_Thorax_Amplitude(images_MoM_SBL_structed(imm),polyelement_MoM);
    TA_MoM_SBL(imm)=TA_MoM_SBL(imm)/max(abs(images_MoM_SBL_structed(imm).elem_data-1));
    %%%%%extract reconstructed lung shapes
    images_MoM_SBL_filtered_seq{imm}=images_MoM_SBL_filtered;
    [llcand_MoM_SBL{imm},rlcand_MoM_SBL{imm},centrelements_MoM_SBL]=extract_reconstructed_shapes...
        (images_MoM_SBL_filtered_seq{imm},polyelement_MoM);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%Position Error
    PE_MoM_SBL{imm}=Calculate_PE(polyelement_MoM,llcand_MoM_SBL{imm},rlcand_MoM_SBL{imm},...
        [Xl(imm),Yl(imm)], [Xr(imm),Yr(imm)],centrelements_MoM);
    %%%%RESOLUTION
    [RES_MoM_SBL{imm},llung_MoM_SBL_rec_area(imm),rlung_MoM_SBL_rec_area(imm)]=...
        Calculate_RES(polyelement_MoM,llcand_MoM_SBL{imm},rlcand_MoM_SBL{imm});
    %%%%Shape Deformation (for each lung particularly and overall SD)
    SD_MoM_SBL{imm}=Calculate_SD(polyelement_MoM,llcand_MoM_SBL{imm},rlcand_MoM_SBL{imm},...
        squeeze(LLUNGGN(imm+1,:,:)),squeeze(RLUNGGN(imm+1,:,:)),llung_MoM_SBL_rec_area(imm),rlung_MoM_SBL_rec_area(imm));
    RNG_MoM_SBL(imm)=Calculate_RNG(polyelement_MoM,llcand_MoM_SBL{imm},rlcand_MoM_SBL{imm},images_MoM_SBL_structed(imm));
end

for imm=1:4
    %%%%calculate CCs
    RR=corrcoef(real(images_diff_smooth(imm).elem_data),inv_reference(imm).elem_data);
    CC_Diff_smooth(imm)=RR(1,2);
    
    RR=corrcoef(real(images_diff_TV(imm).elem_data),inv_reference(imm).elem_data);
    CC_Diff_TV(imm)=RR(1,2);
    
    RR=corrcoef(real(images_diff_move(imm).elem_data),inv_reference(imm).elem_data);
    CC_move(imm)=RR(1,2);
    
    RR=corrcoef(real(imagesABS(imm).elem_data),inv_reference(imm).elem_data);
    CC_abs(imm)=RR(1,2);
    
    RR=corrcoef(real(imagesAll(imm).elem_data),inv_reference(imm).elem_data);
    CC_Non_Linear_All(imm)=RR(1,2);
    
    RR=corrcoef(real(images_MoM_structed(imm).elem_data),...
        real(inv_reference_MoM(imm,:)));
    CC_MoM(imm)=RR(1,2);
    
    RR=corrcoef(real(images_MoM_SBL_structed(imm).elem_data),...
        real(inv_reference_MoM(imm,:)));
    CC_MoM_SBL(imm)=RR(1,2);
end

for imm=1:4
    %%%%calculate GFR
    
    %%%%GN
    DMAX=max(abs(images_diff_smooth(imm).elem_data));
    imgelemsnorm=(images_diff_smooth(imm).elem_data+DMAX)*2/(2*DMAX)-1;
    DMAX1=max(abs(inv_reference(imm).elem_data));
    refelemsnorm=(inv_reference(imm).elem_data+DMAX1)*2/(2*DMAX1)-1;
    GFR_GN_smooth(imm)=calculate_FR(imgelemsnorm,refelemsnorm);
    RMSE_Diff_smooth(imm)=sqrt(immse(imgelemsnorm,refelemsnorm));
    %%%TV
    DMAX=max(abs(images_diff_TV(imm).elem_data));
    imgelemsnorm=(images_diff_TV(imm).elem_data+DMAX)*2/(2*DMAX)-1;
    DMAX1=max(abs(inv_reference(imm).elem_data));
    refelemsnorm=(inv_reference(imm).elem_data+DMAX1)*2/(2*DMAX1)-1;
    GFR_diff_TV(imm)=calculate_FR(imgelemsnorm,refelemsnorm);
    RMSE_Diff_TV(imm)=sqrt(immse(imgelemsnorm,refelemsnorm));
    
    
    %%%ABS
    DMAX=max(abs(imagesABS(imm).elem_data));
    imgelemsnorm=(imagesABS(imm).elem_data+DMAX)*2/(2*DMAX)-1;
    DMAX1=max(abs(inv_reference(imm).elem_data));
    refelemsnorm=(inv_reference(imm).elem_data+DMAX1)*2/(2*DMAX1)-1;
    GFR_ABS(imm)=calculate_FR(imgelemsnorm,refelemsnorm);
    RMSE_abs(imm)=sqrt(immse(imgelemsnorm,refelemsnorm));
    
    %%%ABS
    DMAX=max(abs(images_diff_move(imm).elem_data));
    imgelemsnorm=(images_diff_move(imm).elem_data+DMAX)*2/(2*DMAX)-1;
    DMAX1=max(abs(inv_reference(imm).elem_data));
    refelemsnorm=(inv_reference(imm).elem_data+DMAX1)*2/(2*DMAX1)-1;
    GFR_move(imm)=calculate_FR(imgelemsnorm,refelemsnorm);
    RMSE_move(imm)=sqrt(immse(imgelemsnorm,refelemsnorm));
    
    %%%ImagesAll
    DMAX=max(abs(imagesAll(imm).elem_data));
    imgelemsnorm=(imagesAll(imm).elem_data+DMAX)*2/(2*DMAX)-1;
    DMAX1=max(abs(inv_reference(imm).elem_data));
    refelemsnorm=(inv_reference(imm).elem_data+DMAX1)*2/(2*DMAX1)-1;
    GFR_ROIAll(imm)=calculate_FR(imgelemsnorm,refelemsnorm);
    RMSE_Non_Linear_All(imm)=sqrt(immse(imgelemsnorm,refelemsnorm));
    
    
    %%%PMM-MoM
    images_MoM_structed2=images_MoM_structed(imm);
    images_MoM_structed2.elem_data=images_MoM_structed2.elem_data-1;
    DMAX=max(abs(images_MoM_structed2.elem_data));
    imgelemsnorm=(images_MoM_structed2.elem_data+DMAX)*2/(2*DMAX)-1;
    DMAX1=max(abs(real(inv_reference_MoM(imm,:))));
    refelemsnorm=(real(inv_reference_MoM(imm,:))+DMAX1)*2/(2*DMAX1)-1;
    GFR_MoM(imm)=calculate_FR(imgelemsnorm,refelemsnorm);
    RMSE_MoM(imm)=sqrt(immse(imgelemsnorm,refelemsnorm));
    
    %%%PMM-MoM_SBL
    images_MoM_SBL_structed2=images_MoM_SBL_structed(imm);
    images_MoM_SBL_structed2.elem_data=images_MoM_SBL_structed2.elem_data-1;
    DMAX=max(abs(images_MoM_SBL_structed2.elem_data));
    imgelemsnorm=(images_MoM_SBL_structed2.elem_data+DMAX)*2/(2*DMAX)-1;
    DMAX1=max(abs(real(inv_reference_MoM(imm,:))));
    refelemsnorm=(real(inv_reference_MoM(imm,:))+DMAX1)*2/(2*DMAX1)-1;
    GFR_MoM_SBL(imm)=calculate_FR(imgelemsnorm,refelemsnorm);
    RMSE_MoM_SBL(imm)=sqrt(immse(imgelemsnorm,refelemsnorm));
    
    
end
%%%LOCAL FRS
for imm=1:4
    [RLind,LLind]=find_lungs_indices(polyelement,squeeze(LLUNGGN(imm+1,:,:)),squeeze(RLUNGGN(imm+1,:,:)));
    [RLind_MoM,LLind_MoM]=find_lungs_indices(polyelement_MoM,squeeze(LLUNGGN(imm+1,:,:)),squeeze(RLUNGGN(imm+1,:,:)));
    
    %%%FRLL
    
    %%%%GN
    DMAX=max(abs(images_diff_smooth(imm).elem_data(LLind)));
    imgelemsnorm=(images_diff_smooth(imm).elem_data(LLind)+DMAX)*2/(2*DMAX)-1;
    DMAX1=max(abs(inv_reference(imm).elem_data(LLind)));
    refelemsnorm=(inv_reference(imm).elem_data(LLind)+DMAX1)*2/(2*DMAX1)-1;
    FRLL_GN_smooth(imm)=calculate_FR(imgelemsnorm,refelemsnorm);
    
    %%%TV
    DMAX=max(abs(images_diff_TV(imm).elem_data(LLind)));
    imgelemsnorm=(images_diff_TV(imm).elem_data(LLind)+DMAX)*2/(2*DMAX)-1;
    DMAX1=max(abs(inv_reference(imm).elem_data(LLind)));
    refelemsnorm=(inv_reference(imm).elem_data(LLind)+DMAX1)*2/(2*DMAX1)-1;
    FRLL_diff_TV(imm)=calculate_FR(imgelemsnorm,refelemsnorm);
    
    
    %%%move
    DMAX=max(abs(images_diff_move(imm).elem_data(LLind)));
    imgelemsnorm=(images_diff_move(imm).elem_data(LLind)+DMAX)*2/(2*DMAX)-1;
    DMAX1=max(abs(inv_reference(imm).elem_data(LLind)));
    refelemsnorm=(inv_reference(imm).elem_data(LLind)+DMAX1)*2/(2*DMAX1)-1;
    FRLL_move(imm)=calculate_FR(imgelemsnorm,refelemsnorm);
    
    %%%ABS
    DMAX=max(abs(imagesABS(imm).elem_data(LLind)));
    imgelemsnorm=(imagesABS(imm).elem_data(LLind)+DMAX)*2/(2*DMAX)-1;
    DMAX1=max(abs(inv_reference(imm).elem_data(LLind)));
    refelemsnorm=(inv_reference(imm).elem_data(LLind)+DMAX1)*2/(2*DMAX1)-1;
    FRLL_ABS(imm)=calculate_FR(imgelemsnorm,refelemsnorm);
    
    
    %%%ImagesAll
    DMAX=max(abs(imagesAll(imm).elem_data(LLind)));
    imgelemsnorm=(imagesAll(imm).elem_data(LLind)+DMAX)*2/(2*DMAX)-1;
    DMAX1=max(abs(inv_reference(imm).elem_data(LLind)));
    refelemsnorm=(inv_reference(imm).elem_data(LLind)+DMAX1)*2/(2*DMAX1)-1;
    FRLL_ROIAll(imm)=calculate_FR(imgelemsnorm,refelemsnorm);
    
    %%%PMM-MoM
    images_MoM_structed2=images_MoM_structed(imm);
    images_MoM_structed2.elem_data=images_MoM_structed2.elem_data-1;
    DMAX=max(abs(images_MoM_structed2.elem_data(LLind_MoM)));
    imgelemsnorm=(images_MoM_structed2.elem_data(LLind_MoM)+DMAX)*2/(2*DMAX)-1;
    DMAX1=max(abs(real(inv_reference_MoM(imm,LLind_MoM))));
    refelemsnorm=(real(inv_reference_MoM(imm,LLind_MoM))+DMAX1)*2/(2*DMAX1)-1;
    FRLL_MoM(imm)=calculate_FR(imgelemsnorm,refelemsnorm);
    
    %%%PMM-MoM_SBL
    images_MoM_SBL_structed2=images_MoM_SBL_structed(imm);
    images_MoM_SBL_structed2.elem_data=images_MoM_SBL_structed2.elem_data-1;
    DMAX=max(abs(images_MoM_SBL_structed2.elem_data(LLind_MoM)));
    imgelemsnorm=(images_MoM_SBL_structed2.elem_data(LLind_MoM)+DMAX)*2/(2*DMAX)-1;
    DMAX1=max(abs(real(inv_reference_MoM(imm,LLind_MoM))));
    refelemsnorm=(real(inv_reference_MoM(imm,LLind_MoM))+DMAX1)*2/(2*DMAX1)-1;
    FRLL_MoM_SBL(imm)=calculate_FR(imgelemsnorm,refelemsnorm);
    
    
    %%%FRRL
    
    %%%%GN
    DMAX=max(abs(images_diff_smooth(imm).elem_data(RLind)));
    imgelemsnorm=(images_diff_smooth(imm).elem_data(RLind)+DMAX)*2/(2*DMAX)-1;
    DMAX1=max(abs(inv_reference(imm).elem_data(RLind)));
    refelemsnorm=(inv_reference(imm).elem_data(RLind)+DMAX1)*2/(2*DMAX1)-1;
    FRRL_GN_smooth(imm)=calculate_FR(imgelemsnorm,refelemsnorm);
    
    %%%TV
    DMAX=max(abs(images_diff_TV(imm).elem_data(RLind)));
    imgelemsnorm=(images_diff_TV(imm).elem_data(RLind)+DMAX)*2/(2*DMAX)-1;
    DMAX1=max(abs(inv_reference(imm).elem_data(RLind)));
    refelemsnorm=(inv_reference(imm).elem_data(RLind)+DMAX1)*2/(2*DMAX1)-1;
    FRRL_diff_TV(imm)=calculate_FR(imgelemsnorm,refelemsnorm);
    
    
    %%%move
    DMAX=max(abs(images_diff_move(imm).elem_data(RLind)));
    imgelemsnorm=(images_diff_move(imm).elem_data(RLind)+DMAX)*2/(2*DMAX)-1;
    DMAX1=max(abs(inv_reference(imm).elem_data(RLind)));
    refelemsnorm=(inv_reference(imm).elem_data(RLind)+DMAX1)*2/(2*DMAX1)-1;
    FRRL_move(imm)=calculate_FR(imgelemsnorm,refelemsnorm);
    
    %%%ABS
    DMAX=max(abs(imagesABS(imm).elem_data(RLind)));
    imgelemsnorm=(imagesABS(imm).elem_data(RLind)+DMAX)*2/(2*DMAX)-1;
    DMAX1=max(abs(inv_reference(imm).elem_data(RLind)));
    refelemsnorm=(inv_reference(imm).elem_data(RLind)+DMAX1)*2/(2*DMAX1)-1;
    FRRL_ABS(imm)=calculate_FR(imgelemsnorm,refelemsnorm);
    
    
    %%%ImagesAll
    DMAX=max(abs(imagesAll(imm).elem_data(RLind)));
    imgelemsnorm=(imagesAll(imm).elem_data(RLind)+DMAX)*2/(2*DMAX)-1;
    DMAX1=max(abs(inv_reference(imm).elem_data(RLind)));
    refelemsnorm=(inv_reference(imm).elem_data(RLind)+DMAX1)*2/(2*DMAX1)-1;
    FRRL_ROIAll(imm)=calculate_FR(imgelemsnorm,refelemsnorm);
    
    
    %%%PMM-MoM
    images_MoM_structed2=images_MoM_structed(imm);
    images_MoM_structed2.elem_data=images_MoM_structed2.elem_data-1;
    DMAX=max(abs(images_MoM_structed2.elem_data(RLind_MoM)));
    imgelemsnorm=(images_MoM_structed2.elem_data(RLind_MoM)+DMAX)*2/(2*DMAX)-1;
    DMAX1=max(abs(real(inv_reference_MoM(imm,RLind_MoM))));
    refelemsnorm=(real(inv_reference_MoM(imm,RLind_MoM))+DMAX1)*2/(2*DMAX1)-1;
    FRRL_MoM(imm)=calculate_FR(imgelemsnorm,refelemsnorm);
    
    %%%PMM-MoM-SBL
    images_MoM_SBL_structed2=images_MoM_SBL_structed(imm);
    images_MoM_structed2_SBL.elem_data=images_MoM_SBL_structed2.elem_data-1;
    DMAX=max(abs(images_MoM_structed2_SBL.elem_data(RLind_MoM)));
    imgelemsnorm=(images_MoM_structed2_SBL.elem_data(RLind_MoM)+DMAX)*2/(2*DMAX)-1;
    DMAX1=max(abs(real(inv_reference_MoM(imm,RLind_MoM))));
    refelemsnorm=(real(inv_reference_MoM(imm,RLind_MoM))+DMAX1)*2/(2*DMAX1)-1;
    FRRL_MoM_SBL(imm)=calculate_FR(imgelemsnorm,refelemsnorm);
    
end


o=[];


FR=size(images_diff_move,2);

%%%%%%%image TA
figure
plot(TA_diff_smooth,'-o','LineWidth',2)
hold on
plot(TA_diff_TV,'-o','LineWidth',2)
hold on
plot(TA_movement,'-+','LineWidth',2)
hold on
plot(TA_ABS,'-+','LineWidth',2)
hold on
plot(TA_Non_Linear_All,'-*','LineWidth',2)
hold on
plot(TA_MoM,'-d','LineWidth',2)
hold on
plot(TA_MoM_SBL,'-bd','LineWidth',2)

L1=legend({'GN','TV','Move Laplace','Diff. of Abs','N.L.D.',...
    'M.o.M. Laplace','M.o.M. SBL'},'Interpreter','Latex');

xlim([1 FR])
xticks(1:1:FR);
MaxTA=max([max(abs(TA_diff_smooth)),max(abs(TA_diff_TV)),max(abs(TA_movement)),...
max(abs(TA_ABS)),max(abs(TA_Non_Linear_All))...
,max(abs(TA_MoM)),max(abs(TA_MoM_SBL))]);

MinTA=min([min(abs(TA_diff_smooth)),min(abs(TA_diff_TV)),min(abs(TA_movement)),...
min(abs(TA_ABS)),min(abs(TA_Non_Linear_All))...
,min(abs(TA_MoM)),min(abs(TA_MoM_SBL))]);

ylim([-MaxTA-50 -MinTA+50])

xl=xlabel('Image Frame','Interpreter','Latex');
ax = gca;
ax.FontSize = 16;
yl=ylabel('Target Amplitude ($TA$)','Interpreter','Latex');
yl.FontSize=16;
ay = gca;
ay.FontSize = 16;
set(gca,'TickLabelInterpreter', 'latex');
grid on
TT=title(['$TA$, Thoracic Cavity ' num2str(cavity)],'Interpreter','Latex');
TT.FontSize=18;


%%%image PE
for frame=1:FR
    PE_diff_smooth_Total(frame)=PE_diff_smooth{frame}.Total;
    PE_diff_TV_Total(frame)=PE_diff_TV{frame}.Total;
    PE_ABS_Total(frame)=PE_ABS{frame}.Total;
    PE_movement_Total(frame)=PE_movement{frame}.Total;
    PE_Non_Linear_All_Total(frame)=PE_Non_Linear_All{frame}.Total;
    PE_MoM_Total(frame)=PE_MoM{frame}.Total;
    PE_MoM_SBL_Total(frame)=PE_MoM_SBL{frame}.Total;
end

figure
plot(PE_diff_smooth_Total,'-o','LineWidth',2)
hold on
plot(PE_diff_TV_Total,'-o','LineWidth',2)
hold on
plot(PE_movement_Total,'-+','LineWidth',2)
hold on
plot(PE_ABS_Total,'-+','LineWidth',2)
hold on
plot(PE_Non_Linear_All_Total,'-*','LineWidth',2)
hold on
plot(PE_MoM_Total,'-d','LineWidth',2)
hold on
plot(PE_MoM_SBL_Total,'-bd','LineWidth',2)

L1=legend({'GN','TV','Move Laplace','Diff. of Abs','N.L.D.',...
    'M.o.M. Laplace','M.o.M. SBL'},'Interpreter','Latex');

xlim([1 FR])
xticks(1:1:FR);
xl=xlabel('Image Frame','Interpreter','Latex');
ax = gca;
ax.FontSize = 16;
yl=ylabel('Position Error ($PE$)','Interpreter','Latex');
yl.FontSize=16;
ay = gca;
ay.FontSize = 16;
set(gca,'TickLabelInterpreter', 'latex');
grid on
TT=title(['$PE$, Thoracic Cavity ' num2str(cavity)],'Interpreter','Latex');
TT.FontSize=18;

MaxPE=max([max(abs(PE_diff_smooth_Total)),max(abs(PE_diff_TV_Total)),max(abs(PE_movement_Total)),...
max(abs(PE_ABS_Total)),max(abs(PE_Non_Linear_All_Total))...
,max(abs(PE_MoM_Total)),max(abs(PE_MoM_SBL_Total))]);


ylim([0 MaxPE+MaxPE/5])


%%%image RES
for frame=1:FR
    RES_diff_smooth_Total(frame)=RES_diff_smooth{frame}.Total;
    RES_diff_TV_Total(frame)=RES_diff_TV{frame}.Total;
    RES_movement_Total(frame)=RES_movement{frame}.Total;
    RES_ABS_Total(frame)=RES_ABS{frame}.Total;
    RES_Non_Linear_All_Total(frame)=RES_Non_Linear_All{frame}.Total;
    RES_MoM_Total(frame)=RES_MoM{frame}.Total;
    RES_MoM_SBL_Total(frame)=RES_MoM_SBL{frame}.Total;
end

figure
plot(RES_diff_smooth_Total,'-o','LineWidth',2)
hold on
plot(RES_diff_TV_Total,'-o','LineWidth',2)
hold on
plot(RES_movement_Total,'-+','LineWidth',2)
hold on
plot(RES_ABS_Total,'-+','LineWidth',2)
hold on

plot(RES_Non_Linear_All_Total,'-*','LineWidth',2)
hold on
plot(RES_MoM_Total,'-d','LineWidth',2)
hold on
plot(RES_MoM_SBL_Total,'-bd','LineWidth',2)

L1=legend({'GN','TV','Move Laplace','Diff. of Abs','N.L.D.',...
    'M.o.M. Laplace','M.o.M. SBL'},'Interpreter','Latex');

xlim([1 FR])
xticks(1:1:FR);
xl=xlabel('Image Frame','Interpreter','Latex');
ax = gca;
ax.FontSize = 16;
yl=ylabel('Resolution ($RES$)','Interpreter','Latex');
yl.FontSize=16;
ay = gca;
ay.FontSize = 16;
set(gca,'TickLabelInterpreter', 'latex');
grid on
TT=title(['$RES$, Thoracic Cavity ' num2str(cavity)],'Interpreter','Latex');
TT.FontSize=18;

MaxRES=max([max(abs(RES_diff_smooth_Total)),max(abs(RES_diff_TV_Total)),max(abs(RES_movement_Total)),...
max(abs(RES_ABS_Total)),max(abs(RES_Non_Linear_All_Total))...
,max(abs(RES_MoM_Total)),max(abs(RES_MoM_SBL_Total))]);


ylim([0 MaxRES+MaxRES/5])


%%%image SD
for frame=1:FR
    SD_diff_smooth_Total(frame)=SD_diff_smooth{frame}.Total;
    SD_diff_TV_Total(frame)=SD_diff_TV{frame}.Total;
    SD_movement_Total(frame)=SD_movement{frame}.Total;
    SD_ABS_Total(frame)=SD_ABS{frame}.Total;
    SD_Non_Linear_All_Total(frame)=SD_Non_Linear_All{frame}.Total;
    SD_MoM_Total(frame)=SD_MoM{frame}.Total;
    SD_MoM_SBL_Total(frame)=SD_MoM_SBL{frame}.Total;
end

figure
plot(SD_diff_smooth_Total,'-o','LineWidth',2)
hold on
plot(SD_diff_TV_Total,'-o','LineWidth',2)
hold on
plot(SD_movement_Total,'-+','LineWidth',2)
hold on
plot(SD_ABS_Total,'-+','LineWidth',2)
hold on
plot(SD_Non_Linear_All_Total,'-*','LineWidth',2)
hold on
plot(SD_MoM_Total,'-d','LineWidth',2)
hold on
plot(SD_MoM_SBL_Total,'-bd','LineWidth',2)

L1=legend({'GN','TV','Move Laplace','Diff. of Abs','N.L.D.',...
    'M.o.M. Laplace','M.o.M. SBL'},'Interpreter','Latex');

xlim([1 FR])
xticks(1:1:FR);
xl=xlabel('Image Frame','Interpreter','Latex');
ax = gca;
ax.FontSize = 16;
yl=ylabel('Shape Deformation ($SD$)','Interpreter','Latex');
yl.FontSize=16;
ay = gca;
ay.FontSize = 16;
set(gca,'TickLabelInterpreter', 'latex');
grid on
TT=title(['$SD$, Thoracic Cavity ' num2str(cavity)],'Interpreter','Latex');
TT.FontSize=18;

MaxSD=max([max(abs(SD_diff_smooth_Total)),max(abs(SD_diff_TV_Total)),max(abs(SD_movement_Total)),...
max(abs(SD_ABS_Total)),max(abs(SD_Non_Linear_All_Total))...
,max(abs(SD_MoM_Total)),max(abs(SD_MoM_SBL_Total))]);


ylim([0 MaxSD+MaxSD/5])

%%%image RNG

figure
plot(RNG_diff_smooth,'-o','LineWidth',2)
hold on
plot(RNG_diff_TV,'-o','LineWidth',2)
hold on
plot(RNG_movement,'-+','LineWidth',2)
hold on
plot(RNG_ABS,'-+','LineWidth',2)
hold on
plot(RNG_Non_Linear_All,'-*','LineWidth',2)
hold on
plot(RNG_MoM,'-d','LineWidth',2)
hold on
plot(RNG_MoM_SBL,'-bd','LineWidth',2)

L1=legend({'GN','TV','Move Laplace','Diff. of Abs','N.L.D.',...
    'M.o.M. Laplace','M.o.M. SBL'},'Interpreter','Latex');

xlim([1 FR])
xticks(1:1:FR);
xl=xlabel('Image Frame','Interpreter','Latex');
ax = gca;
ax.FontSize = 16;
yl=ylabel('Ringing ($RNG$)','Interpreter','Latex');
yl.FontSize=16;
ay = gca;
ay.FontSize = 16;
set(gca,'TickLabelInterpreter', 'latex');
grid on
TT=title(['$RNG$, Thoracic Cavity ' num2str(cavity)],'Interpreter','Latex');
TT.FontSize=18;

MaxRNG=max([max(abs(RNG_diff_smooth)),max(abs(RNG_diff_TV)),max(abs(RNG_movement)),...
max(abs(RNG_ABS)),max(abs(RNG_Non_Linear_All))...
,max(abs(RNG_MoM)),max(abs(RNG_MoM_SBL))]);


ylim([0 MaxRNG+MaxRNG/5])

figure

plot(CC_Diff_smooth,'-o','LineWidth',2)
hold on
plot(CC_Diff_TV,'-o','LineWidth',2)
hold on
plot(CC_move,'-+','LineWidth',2)
hold on
plot(CC_abs,'-+','LineWidth',2)
hold on
plot(CC_Non_Linear_All,'-*','LineWidth',2)
hold on
plot(CC_MoM,'-d','LineWidth',2)
hold on
plot(CC_MoM_SBL,'-bd','LineWidth',2)

L1=legend({'GN','TV','Move Laplace','Diff. of Abs','N.L.D.',...
    'M.o.M. Laplace','M.o.M. SBL'},'Interpreter','Latex');
L1.FontSize=14;

xlim([1 FR])
xticks(1:1:FR);
xl=xlabel('Image Frame','Interpreter','Latex');
ax = gca;
ax.FontSize = 16;
yl=ylabel('Correlation Coeff. ($CC$)','Interpreter','Latex');
yl.FontSize=16;
ay = gca;
ay.FontSize = 16;
set(gca,'TickLabelInterpreter', 'latex');
grid on
TT=title(['$CC$, Thoracic Cavity ' num2str(cavity)],'Interpreter','Latex');
TT.FontSize=18;
ylim([0 1])


figure
%%%RMSE
plot(RMSE_Diff_smooth,'-o','LineWidth',2)
hold on
plot(RMSE_Diff_TV,'-o','LineWidth',2)
hold on
plot(RMSE_move,'-+','LineWidth',2)
hold on
plot(RMSE_abs,'-+','LineWidth',2)
hold on
plot(RMSE_Non_Linear_All,'-*','LineWidth',2)
hold on
plot(RMSE_MoM,'-d','LineWidth',2)
hold on
plot(RMSE_MoM_SBL,'-bd','LineWidth',2)

L1=legend({'GN','TV','Move Laplace','Diff. of Abs','N.L.D.',...
    'M.o.M. Laplace','M.o.M. SBL'},'Interpreter','Latex');

xlim([1 FR])

xticks(1:1:FR);

xl=xlabel('Image Frame','Interpreter','Latex');
ax = gca;
ax.FontSize = 16;
yl=ylabel('$RMSE$','Interpreter','Latex');
yl.FontSize=16;
ay = gca;
ay.FontSize = 16;
set(gca,'TickLabelInterpreter', 'latex');
grid on
TT=title(['$RMSE$, Thoracic Cavity ' num2str(cavity)],'Interpreter','Latex');
TT.FontSize=18;

ylim([0 1])

figure
%%%%GFR
plot(GFR_GN_smooth,'-o','LineWidth',2)
hold on
plot(GFR_diff_TV,'-o','LineWidth',2)
hold on
plot(GFR_move,'-+','LineWidth',2)
hold on
plot(GFR_ABS,'-+','LineWidth',2)
hold on
plot(GFR_ROIAll,'-*','LineWidth',2)
hold on
plot(GFR_MoM,'-d','LineWidth',2)
hold on
plot(GFR_MoM_SBL,'-bd','LineWidth',2)

L1=legend({'GN','TV','Move Laplace','Diff. of Abs','N.L.D.',...
    'M.o.M. Laplace','M.o.M. SBL'},'Interpreter','Latex');

xlim([1 FR])
xticks(1:1:FR);

xl=xlabel('Image Frame','Interpreter','Latex');
ax = gca;
ax.FontSize = 16;
yl=ylabel('$GFR$','Interpreter','Latex');
yl.FontSize=16;
ay = gca;
ay.FontSize = 16;
set(gca,'TickLabelInterpreter', 'latex');
grid on
TT=title(['$GFR$, Thoracic Cavity ' num2str(cavity)],'Interpreter','Latex');
TT.FontSize=18;


%%%%FR_LL
figure
%%%%GFR
plot(FRLL_GN_smooth,'-o','LineWidth',2)
hold on
plot(FRLL_diff_TV,'-o','LineWidth',2)
hold on
plot(FRLL_move,'-+','LineWidth',2)
hold on
plot(FRLL_ABS,'-+','LineWidth',2)
hold on
plot(FRLL_ROIAll,'-*','LineWidth',2)
hold on
plot(FRLL_MoM,'-d','LineWidth',2)
hold on
plot(FRLL_MoM_SBL,'-bd','LineWidth',2)

L1=legend({'GN','TV','Move Laplace','Diff. of Abs','N.L.D.',...
    'M.o.M. Laplace','M.o.M. SBL'},'Interpreter','Latex');

xlim([1 FR])
xticks(1:1:FR);

xl=xlabel('Image Frame','Interpreter','Latex');
ax = gca;
ax.FontSize = 16;
yl=ylabel('$FR_{LL}$','Interpreter','Latex');
yl.FontSize=16;
ay = gca;
ay.FontSize = 16;
set(gca,'TickLabelInterpreter', 'latex');
grid on
TT=title(['$FR_{LL}$, Thoracic Cavity ' num2str(cavity)],'Interpreter','Latex');
TT.FontSize=18;

%%%%FRRL
figure
%%%%GFR
plot(FRRL_GN_smooth,'-o','LineWidth',2)
hold on
plot(FRRL_diff_TV,'-o','LineWidth',2)
hold on
plot(FRRL_move,'-+','LineWidth',2)
hold on
plot(FRRL_ABS,'-+','LineWidth',2)
hold on
plot(FRRL_ROIAll,'-*','LineWidth',2)
hold on
plot(FRRL_MoM,'-d','LineWidth',2)
hold on
plot(FRRL_MoM_SBL,'-bd','LineWidth',2)

L1=legend({'GN','TV','Move Laplace','Diff. of Abs','N.L.D.',...
    'M.o.M. Laplace','M.o.M. SBL'},'Interpreter','Latex');

xlim([1 FR])
xticks(1:1:FR);

xl=xlabel('Image Frame','Interpreter','Latex');
ax = gca;
ax.FontSize = 16;
yl=ylabel('$FR_{RL}$','Interpreter','Latex');
yl.FontSize=16;
ay = gca;
ay.FontSize = 16;
set(gca,'TickLabelInterpreter', 'latex');
grid on
TT=title(['$FR_{RL}$, Thoracic Cavity ' num2str(cavity)],'Interpreter','Latex');
TT.FontSize=18;

%%%%%%mean+/-STD

%%%%mean TA
TAMoM_m=mean(TA_MoM);
TAMoM_SBL_m=mean(TA_MoM_SBL);
TAabs_m=mean(TA_ABS);
TAGN_m=mean(TA_diff_smooth);
TATV_m=mean(TA_diff_TV);
TAmove_m=mean(TA_movement);
TANonLinearAllm=mean(TA_Non_Linear_All);


%%%%std TA
TAMoM_std=std(TA_MoM);
TAMoM_SBL_std=std(TA_MoM_SBL);
TAabs_std=std(TA_ABS);
TAGN_std=std(TA_diff_smooth);
TATV_std=std(TA_diff_TV);
TAmove_std=std(TA_movement);
TANonLinearAll_std=std(TA_Non_Linear_All);


%%%%mean PE
PEMoM_m=mean(PE_MoM_Total);
PEMoM_SBL_m=mean(PE_MoM_SBL_Total);
PEabs_m=mean(PE_ABS_Total);
PEGN_m=mean(PE_diff_smooth_Total);
PETV_m=mean(PE_diff_TV_Total);
PEmove_m=mean(PE_movement_Total);
PENonLinearAllm=mean(PE_Non_Linear_All_Total);


%%%%std PE
PEMoM_std=std(PE_MoM_Total);
PEMoM_SBL_std=std(PE_MoM_SBL_Total);
PEabs_std=std(PE_ABS_Total);
PEGN_std=std(PE_diff_smooth_Total);
PETV_std=std(PE_diff_TV_Total);
PEmove_std=std(PE_movement_Total);
PENonLinearAll_std=std(PE_Non_Linear_All_Total);

%%%%mean RES
RESMoM_m=mean(RES_MoM_Total);
RESabs_m=mean(RES_ABS_Total);
RESGN_m=mean(RES_diff_smooth_Total);
RESTV_m=mean(RES_diff_TV_Total);
RESmove_m=mean(RES_movement_Total);
RESNonLinearAllm=mean(RES_Non_Linear_All_Total);

%%%%std RES
RESMoM_std=std(RES_MoM_Total);
RESabs_std=std(RES_ABS_Total);
RESGN_std=std(RES_diff_smooth_Total);
RESTV_std=std(RES_diff_TV_Total);
RESmove_std=std(RES_movement_Total);
RESNonLinearAll_std=std(RES_Non_Linear_All_Total);


%%%%mean SD
SDMoM_m=mean(SD_MoM_Total);
SDabs_m=mean(SD_ABS_Total);
SDGN_m=mean(SD_diff_smooth_Total);
SDTV_m=mean(SD_diff_TV_Total);
SDmove_m=mean(SD_movement_Total);
SDNonLinearAllm=mean(SD_Non_Linear_All_Total);

%%%%std SD
SDMoM_std=std(SD_MoM_Total);
SDabs_std=std(SD_ABS_Total);
SDGN_std=std(SD_diff_smooth_Total);
SDTV_std=std(SD_diff_TV_Total);
SDmove_std=std(SD_movement_Total);
SDNonLinearAll_std=std(SD_Non_Linear_All_Total);



%%%%mean RNG
RNGMoM_m=mean(RNG_MoM);
RNGabs_m=mean(RNG_ABS);
RNGGN_m=mean(RNG_diff_smooth);
RNGTV_m=mean(RNG_diff_TV);
RNGmove_m=mean(RNG_movement);
RNGNonLinearAllm=mean(RNG_Non_Linear_All);

%%%%std RNG
RNGMoM_std=std(RNG_MoM);
RNGabs_std=std(RNG_ABS);
RNGGN_std=std(RNG_diff_smooth);
RNGTV_std=std(RNG_diff_TV);
RNGmove_std=std(RNG_movement);
RNGNonLinearAll_std=std(RNG_Non_Linear_All);


%%%%%%mean CC
CCMoM_m=mean(CC_MoM);
CCabs_m=mean(CC_abs);
CCGN_m=mean(CC_Diff_smooth);
CCTV_m=mean(CC_Diff_TV);
CCmove_m=mean(CC_move);
CCNonLinearAllm=mean(CC_Non_Linear_All);

%%%%%%std CC
CCMoM_std=std(CC_MoM);
CCabs_std=std(CC_abs);
CCGN_std=std(CC_Diff_smooth);
CCTV_std=std(CC_Diff_TV);
CCmove_std=std(CC_move);
CCNonLinearAll_std=std(CC_Non_Linear_All);
