Cavity=1;
path2load1=['Thoracic_Cavities\Thorax' num2str(Cavity) '\'];
addpath(path2load1)
path2load2='Models_FOM\Thorax_3D\';
addpath(path2load2)
load(['thorax' num2str(Cavity) '.mat'])
load(['thorax' num2str(Cavity) 'inf.mat'])
load(['thorax' num2str(Cavity) '_medstate.mat'])
thoraxdef=eval(['thorax' num2str(Cavity)]);
load(['rlung' num2str(Cavity) '.mat'])
load(['rlung' num2str(Cavity) 'inf.mat'])
load(['rlung' num2str(Cavity) '_medstate.mat'])
rlungdef=eval(['rlung' num2str(Cavity)]);
load(['llung' num2str(Cavity) '.mat'])
load(['llung' num2str(Cavity) 'inf.mat'])
load(['llung' num2str(Cavity) '_medstate.mat'])
llungdef=eval(['llung' num2str(Cavity)]);
load(['spondilus' num2str(Cavity) '.mat'])
spondilus=eval(['spondilus' num2str(Cavity)]);
load(['newheart' num2str(Cavity) '.mat'])
heart=eval(['newheart' num2str(Cavity)]);


figure
%%%%%boundaries
plot([thoraxdef(:,1); thoraxdef(1,1)],[thoraxdef(:,2); thoraxdef(1,2)],'k','LineWidth',2)
hold on
plot([squeeze(thorax_medstate(2,:,1))'; thorax_medstate(2,1,1)],...
    [squeeze(thorax_medstate(2,:,2))'; thorax_medstate(2,1,2)],'k--','LineWidth',1)
hold on
plot([squeeze(thorax_medstate(3,:,1))'; thorax_medstate(3,1,1)],...
    [squeeze(thorax_medstate(3,:,2))'; thorax_medstate(3,1,2)],'k--','LineWidth',1)
hold on
plot([squeeze(thorax_medstate(4,:,1))'; thorax_medstate(4,1,1)],...
    [squeeze(thorax_medstate(4,:,2))'; thorax_medstate(4,1,2)],'k--','LineWidth',1)
hold on
plot([thoraxinf(:,1); thoraxinf(1,1)],[thoraxinf(:,2); thoraxinf(1,2)],'k--','LineWidth',2)

hold on
%%%%lung1 
plot([llungdef(:,1); llungdef(1,1)],[llungdef(:,2); llungdef(1,2)],'b','LineWidth',2)
hold on
plot([squeeze(llung_medstate(2,:,1))'; llung_medstate(2,1,1)],...
    [squeeze(llung_medstate(2,:,2))'; llung_medstate(2,1,2)],'b--','LineWidth',1)
hold on
plot([squeeze(llung_medstate(3,:,1))'; llung_medstate(3,1,1)],...
    [squeeze(llung_medstate(3,:,2))'; llung_medstate(3,1,2)],'b--','LineWidth',1)
hold on
plot([squeeze(llung_medstate(4,:,1))'; llung_medstate(4,1,1)],...
    [squeeze(llung_medstate(4,:,2))'; llung_medstate(4,1,2)],'b--','LineWidth',1)
hold on
plot([llunginf(:,1); llunginf(1,1)],[llunginf(:,2); llunginf(1,2)],'b--','LineWidth',2)

hold on
%%%%lung2 
plot([rlungdef(:,1); rlungdef(1,1)],[rlungdef(:,2); rlungdef(1,2)],'b','LineWidth',2)
hold on
plot([squeeze(rlung_medstate(2,:,1))'; rlung_medstate(2,1,1)],...
    [squeeze(rlung_medstate(2,:,2))'; rlung_medstate(2,1,2)],'b--','LineWidth',1)
hold on
plot([squeeze(rlung_medstate(3,:,1))'; rlung_medstate(3,1,1)],...
    [squeeze(rlung_medstate(3,:,2))'; rlung_medstate(3,1,2)],'b--','LineWidth',1)
hold on
plot([squeeze(rlung_medstate(4,:,1))'; rlung_medstate(4,1,1)],...
    [squeeze(rlung_medstate(4,:,2))'; rlung_medstate(4,1,2)],'b--','LineWidth',1)
hold on
plot([rlunginf(:,1); rlunginf(1,1)],[rlunginf(:,2); rlunginf(1,2)],'b--','LineWidth',2)

hold on
%%%%heart 
plot([heart(:,1); heart(1,1)],[heart(:,2); heart(1,2)],'r','LineWidth',2)

hold on
%%%spondilus
plot([spondilus(:,1); spondilus(1,1)],[spondilus(:,2); spondilus(1,2)],...
    'Color',[0.4 0.4 0.4] ,'LineWidth',2)
axis square
xlim([-1 1.25])
axis off
Leg=legend({'Boundary Expiration-End','Boundary median states','',...
    '','Boundary Inspiration-End','Lungs Expiration-End','Lungs median states','',...
    '','Lungs Inspiration-End','Lungs Expiration-End','Lungs median states','',...
    '','Lungs Inspiration-End','Heart','Vertebrae'});
axis equal
axis off