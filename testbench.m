%%%%%%%TESTBENCH script for dynamic thoracic imaging using several
%%%%%%%algorithms

%%%%%%Please Cite:
%%%% Dimas, C., Alimisis, V., Uzunoglu, N., & Sotiriadis, P. P. (2021). 
%%%% A point-matching method of moment with sparse Bayesian learning applied 
%%%% and evaluated in dynamic lung electrical impedance tomography. 
%%%% Bioengineering, 8(12), 191.

%%%%% Only a single thoracic cavity (for 5 breathing states) supported in this demo.
%%%%% For more, please contact chdim@central.ntua.gr
%%%%%%% number of electrodes
N=16;
%%%%%%only 16 supported
%%%%%%% current and voltage skip patterns
%%%%%%% (measurements protocol)
skipcurr=0;
skipvolt=0;
%%%%%%% desired SNR (noise deliberately added to
%%%%%%% check each algorithm's robustness)
SNR=50;
%%%%%% type of data: Thorax3D for dynamic 5-state 3D thoracic simulation
%%%%%% model.
%%%%%% In-Vivo for in-vivo data (resource:EIDORS)
type='Thorax3D';
%type='In-Vivo';
%%%%%% if type==Thorax3D select cavity (thoracic models) from 1-5, else
%%%%%% else indifferent
cavity=1;
if cavity>1
    %%%%%% works in original code, memory issues in Github!
    error('This demo works only for the 1st cavity! Please select cavity=1')
end
%%%%% create new model, default is set to zero (See Regularized_Diff_EIT)
tp5D=0;

%%% collect the data
[vhomo,vinhomo,sim_model]=data_collection(N,skipcurr,skipvolt,SNR,type,cavity);


if strcmp(type,'Thorax3D')
    sim_model_real=sim_model;
    for i=1:5
        sim_model_real{i}.elem_data=abs(sim_model_real{i}.elem_data);
    end
    figure
    %subplot(1,5,1)
    sim_model_real{1}.calc_colours.ref_level= 0.25;
    sim_model_real{1}.calc_colours.cb_shrink_move = [0.5,0.8,-.10];
    H1=show_fem(sim_model_real{1},1);
    set(H1, 'edgecolor', [0.5 0.5 0.5]);
    set(H1, 'edgealpha', 0.2);
    axis off
    view([0 80])
    title('Deflated')
    
    %subplot(1,5,2)
    figure
    sim_model_real{2}.calc_colours.ref_level= 0.25;
    sim_model_real{2}.calc_colours.cb_shrink_move = [0.5,0.8,-.10];
    H2=show_fem(sim_model_real{2},1);
    set(H2, 'edgecolor', [0.5 0.5 0.5]);
    set(H2, 'edgealpha', 0.25);
    axis off
    view([0 80])
    title('1/4 fill')
    
    
    
    %subplot(1,5,3)
    figure
    sim_model_real{3}.calc_colours.ref_level= 0.25;
    sim_model_real{3}.calc_colours.cb_shrink_move = [0.5,0.8,-.10];
    H3=show_fem(sim_model_real{3},1);
    set(H3, 'edgecolor', [0.5 0.5 0.5]);
    set(H3, 'edgealpha', 0.25);
    axis off
    view([0 80])
    title('2/4 fill')
    
    
    %subplot(1,5,4)
    figure
    sim_model_real{4}.calc_colours.ref_level= 0.25;
    sim_model_real{4}.calc_colours.cb_shrink_move = [0.5,0.8,-.10];
    H4=show_fem(sim_model_real{4},1);
    set(H4, 'edgecolor', [0.5 0.5 0.5]);
    set(H4, 'edgealpha', 0.25);
    axis off
    view([0 80])
    title('3/4 fill')
    
    %subplot(1,5,5)
    figure
    sim_model_real{5}.calc_colours.ref_level= 0.25;
    sim_model_real{5}.calc_colours.cb_shrink_move = [0.5,0.8,-.10];
    H5=show_fem(sim_model_real{5},1);
    set(H5, 'edgecolor', [0.5 0.5 0.5]);
    set(H5, 'edgealpha', 0.25);
    axis off
    view([0 80])
    title('Inflated')
    
end

%%%%%%PARAMETERIZATION
if strcmp(type,'Thorax3D')
    %%%%Thoracic3D Laplace
    lambda_smooth=0.005;
    %%%%
    %%%%Thoracic3D TV
    lambdaTV=1e-06;
    %%%%
elseif strcmp(type,'In-Vivo')
    %%%%In-Vivo Laplace
    lambda_smooth=0.005;
    lambdaTV=1e-06;
end
%%%%

[images_diff_smooth,Q,imdl_diff_smooth]=Regularized_Diff_EIT(vhomo,vinhomo,...
    N,skipcurr,skipvolt,'GN','Laplace',lambda_smooth,tp5D);
[images_diff_TV]=Regularized_Diff_EIT(vhomo,vinhomo,N,skipcurr,skipvolt,'TV','TV',lambdaTV,tp5D);
tic
[images_GREIT,imdlgr,invmodel3D]=GREIT_imaging(vhomo,vinhomo,N,skipcurr,skipvolt,[35 35]);
toc
%images_GREIT.elem_data=real(images_GREIT.elem_data);
%[rr,prms]=test_performance(imdlgr,invmodel3D);

figure
subplot(1,3,1)
show_slices(images_diff_smooth)
title('Iterative GTR')
subplot(1,3,2)
show_slices(images_diff_TV)
title('Total Variation')
subplot(1,3,3)
show_slices(images_GREIT)
title('GREIT')

if strcmp(type,'Thorax3D')
    %%%Thoracic3D ABS
    lambdaABS=0.01;
    %%%
elseif strcmp(type,'In-Vivo')
    %%%In-Vivo ABS
    %%%
    lambdaABS=900;
end

imagesABS=Absolute_EIT(vhomo,vinhomo,N,skipcurr,skipvolt,Q,lambdaABS);
figure
show_slices(imagesABS)
title('Absolute EIT Imaging')

%%% Define ROIs
ROItype='Central';
%ROItype='Lungs and Heart';
%ROItype='All';
[Elems_ROI,LROI,L]=Define_ROIs(N,ROItype);

if strcmp(type,'Thorax3D')
    %%%%Thoracic3D, Laplace
    lambda1=1e-03;
    lambda2=1e-04;
    %%%Thoracic3D beta
    beta=1e-04;
    
elseif strcmp(type,'In-Vivo')
    %%%%In-Vivo, Laplace
    lambda1=50;
    lambda2=10;
    %%%In-Vivo beta
    beta=1e-04;
end

imagesROICentral=Non_Linear_Diff_Imaging(vhomo,vinhomo,N,L,LROI,Elems_ROI,skipcurr,skipvolt,Q,lambda1,lambda2,beta);
figure
show_slices(imagesROICentral)
title('Non-Linear Difference Imaging, ROI=Central')

if strcmp(type,'Thorax3D')
    %%%%Thoracic3D, Laplace
    lambda1=1e-03;
    lambda2=1e-04;
    %%%%
elseif strcmp(type,'In-Vivo')
    %%%%In-Vivo, Laplace
    lambda1=50;
    lambda2=10;
    %%%%
    %%%%% beta does not need to be changed from ROI-to-ROI
end

ROItype='Lungs and Heart';
[Elems_ROI,LROI,L]=Define_ROIs(N,ROItype);
imagesROI_Tissues=Non_Linear_Diff_Imaging(vhomo,vinhomo,N,L,LROI,Elems_ROI,...
    skipcurr,skipvolt,Q,lambda1,lambda2,beta);
figure
show_slices(imagesROI_Tissues)
title('Non-Linear Difference Imaging, ROI=Tissues')

if strcmp(type,'Thorax3D')
    %%%%Thoracic3D, Laplace
    lambda1=1e-03;
    lambda2=1e-04;
    %%%%
elseif strcmp(type,'In-Vivo')
    %%%%In-Vivo, Laplace
    lambda1=50;
    lambda2=10;
    %%%%
end

ROItype='All';
[Elems_ROI,LROI,L]=Define_ROIs(N,ROItype);
imagesAll=Non_Linear_Diff_Imaging(vhomo,vinhomo,N,L,LROI,Elems_ROI,skipcurr,...
    skipvolt,Q,lambda1,lambda2,beta);
figure
show_slices(imagesAll)
title('Non-Linear Difference Imaging, ROI=\Omega')

%%%%%linear single-step with movement prior

fprintf('Movement Prior-Linear Single-Step GN\n')
mu=sqrt(2);
[images_diff_move,Q,imdl_diff_move]=Linear_imaging_movement(vhomo,vinhomo,N,skipcurr,skipvolt,lambda_smooth,mu);
figure
show_slices(images_diff_move)
title('Movement Prior-Linear Reconstruction')

%lambdaMOM=0.13;
%lambdaMOM=0.11;
lambdaMOM=0.2;

mu=0;
%sqrt(2)/100;
tic
[sigmak,xc,yc,elec_nodes,movements]=PMM_EIT(vhomo/100,vinhomo/100,N,skipcurr,skipvolt,lambdaMOM,mu);
toc
figure
for frame=1:size(vinhomo,2)
    if strcmp(type,'Thorax3D')
        subplot(1,5,frame)
        dotsize=25;
    else
        subplot(7,5,frame)
        dotsize=10;
    end
    scatter3(xc,yc,sigmak(frame,:),dotsize,sigmak(frame,:),'filled')
    view([0 90])
    axis tight
    axis off
    colormap jet
    caxis([min(min(sigmak(1:end,:))) max(max(sigmak(1:end,:)))])
end
suptitle('PM-MoM Reconstruction')
