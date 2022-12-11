%%%%%%%%%First run FOM_testbench.m
close all
%%%GN
figure('Renderer', 'painters', 'Position', [10 10 900 600])
%images_diff_smoothreal=images_diff_smooth;
for imm=1:4
    subplot(1,4,imm)
    images_diff_smooth(imm).elem_data=real(images_diff_smooth(imm).elem_data);
    %images_diff_smoothreal(imm).calc_colours.cb_shrink_move = [0.5,0.8,-.10];
    images_diff_smooth(imm).calc_colours.clim= max(abs(images_diff_smooth(4).elem_data))/2;
    %images_diff_smoothreal(imm).calc_colours.sat_adj= max(abs(images_diff_smoothreal(4).elem_data))/1.3;
    H1=show_fem(images_diff_smooth(imm));
    set(H1,'edgecolor', [0.5 0.5 0.5]);
    set(H1,'edgealpha', 0.25);
    %caxis([-1 1])
    axis off
    axis tight
    colormap jet
end

%%%%TV
figure('Renderer', 'painters', 'Position', [10 10 900 600])
for imm=1:4
    subplot(1,4,imm)
    %images_diff_TV(imm).calc_colours.cb_shrink_move = [0.5,0.8,-.10];
    %images_diff_TV(imm).calc_colours.sat_adj= max(abs(images_diff_TV(4).elem_data));
    images_diff_TV(imm).calc_colours.clim= max(abs(images_diff_TV(4).elem_data))/2;
    H1=show_fem(images_diff_TV(imm));
    set(H1,'edgecolor', [0.5 0.5 0.5]);
    set(H1,'edgealpha', 0.25);
    axis off
    axis tight
    colormap jet
end

figure('Renderer', 'painters', 'Position', [10 10 900 600])
%%%ABS
for imm=1:4
    subplot(1,4,imm)
%     imagesABS(imm).calc_colours.cb_shrink_move = [0.5,0.8,-.010];
%     imagesABS(imm).calc_colours.sat_adj= max(abs(imagesABS(4).elem_data))*50;
    imagesABS(imm).calc_colours.clim= max(abs(imagesABS(4).elem_data))/2;
    H1=show_fem(imagesABS(imm));
    set(H1,'edgecolor', [0.5 0.5 0.5]);
    set(H1,'edgealpha', 0.25);
    axis off
    axis tight
    colormap jet
end

figure('Renderer', 'painters', 'Position', [10 10 900 600])
%%%%GREIT
for imm=1:4
    subplot(1,4,imm)
%     images_GREIT(imm).calc_colours.cb_shrink_move = [0.5,0.8,-.010];
%     images_GREIT(imm).calc_colours.sat_adj= max(abs(images_GREIT(4).elem_data))/1.5e+03;
%      images_GREIT(imm).calc_colours.ref_level =  -170;
    images_GREIT(imm).calc_colours.clim= max(abs(images_GREIT(4).elem_data))/2;
    H1=show_fem(images_GREIT(imm));
    set(H1,'edgecolor', [0.5 0.5 0.5]);
    set(H1,'edgealpha', 0.25);
    axis off
    axis tight
    colormap jet
end

figure('Renderer', 'painters', 'Position', [10 10 900 600])
%Movement
for imm=1:4
    subplot(1,4,imm)
%     images_diff_move(imm).calc_colours.cb_shrink_move = [0.5,0.8,-.010];
%     images_diff_move(imm).calc_colours.sat_adj= max(abs(images_diff_move(4).elem_data))/8;
    images_diff_move(imm).calc_colours.clim= max(abs(images_diff_move(4).elem_data))/2;
    H1=show_fem(images_diff_move(imm));
    set(H1,'edgecolor', [0.5 0.5 0.5]);
    set(H1,'edgealpha', 0.25);
    axis off
    axis tight
    colormap jet
end

figure('Renderer', 'painters', 'Position', [10 10 900 600])
%ROI ALL
for imm=1:4
    subplot(1,4,imm)
%     imagesAll(imm).calc_colours.cb_shrink_move = [0.5,0.8,-.010];
%     imagesAll(imm).calc_colours.sat_adj= max(abs(imagesAll(4).elem_data))*30;
    imagesAll(imm).calc_colours.clim= max(abs(imagesAll(4).elem_data))/2;
    H1=show_fem(imagesAll(imm));
    set(H1,'edgecolor', [0.5 0.5 0.5]);
    set(H1,'edgealpha', 0.25);
    axis off
    axis tight
    colormap jet
end

figure('Renderer', 'painters', 'Position', [10 10 900 600])
%ROI Central
for imm=1:4
    subplot(1,4,imm)
%     imagesROICentral(imm).calc_colours.cb_shrink_move = [0.5,0.8,-.010];
%     imagesROICentral(imm).calc_colours.sat_adj= max(abs(imagesROICentral(4).elem_data));
    imagesROICentral(imm).calc_colours.clim= max(abs(imagesROICentral(4).elem_data))/2;
    H1=show_fem(imagesROICentral(imm));
    set(H1,'edgecolor', [0.5 0.5 0.5]);
    set(H1,'edgealpha', 0.25);
    axis off
    axis tight
    colormap jet
end

figure('Renderer', 'painters', 'Position', [10 10 900 600])
%ROI Tissues
for imm=1:4
    subplot(1,4,imm)
%     imagesROICentral(imm).calc_colours.cb_shrink_move = [0.5,0.8,-.010];
%     imagesROICentral(imm).calc_colours.sat_adj= max(abs(imagesROICentral(4).elem_data));
    imagesROI_Tissues(imm).calc_colours.clim= max(abs(imagesROI_Tissues(4).elem_data))/2;
    H1=show_fem(imagesROI_Tissues(imm));
    set(H1,'edgecolor', [0.5 0.5 0.5]);
    set(H1,'edgealpha', 0.25);
    colormap jet

    axis off
    axis tight
end

figure('Renderer', 'painters', 'Position', [10 10 900 600])
pathThoracic='D:\EIT\Thoracic_MOM\';
addpath(pathThoracic)

load([pathThoracic 'xc.mat'])
load([pathThoracic 'yc.mat'])

%%%PMM-MoM
for imm=1:4
 subplot(1,4,imm)
    scatter3(xc,yc,imagesMoM(imm,:),20,imagesMoM(imm,:),'filled')
    view([0 90])
    axis equal
    axis off
    colormap jet
    range=max(max(imagesMoM(1:end,4)))-min(min(imagesMoM(1:end,4)));
    caxis([min(min(imagesMoM(1:end,:))) 1.038*max(max(imagesMoM(1:end,:)))])
%      RR=corrcoef(real(sigmak(frame,:)),...
%      real(inv_reference_MoM(frame,:)));
%      CC_MoM(frame)=RR(1,2);
end

%%%%%INV_REFERENCE
figure('Renderer', 'painters', 'Position', [10 10 900 600])
for imm=1:4
    subplot(1,4,imm)
    inv_reference(imm).elem_data=real(inv_reference(imm).elem_data);
    %images_diff_smoothreal(imm).calc_colours.cb_shrink_move = [0.5,0.8,-.10];
    inv_reference(imm).calc_colours.clim= max(abs(inv_reference(4).elem_data))/2;
    %images_diff_smoothreal(imm).calc_colours.sat_adj= max(abs(images_diff_smoothreal(4).elem_data))/1.3;
    H1=show_fem(inv_reference(imm));
    set(H1,'edgecolor', [0.5 0.5 0.5]);
    set(H1,'edgealpha', 0.01);
    %caxis([-1 1])
    axis off
    axis tight
    colormap jet
end