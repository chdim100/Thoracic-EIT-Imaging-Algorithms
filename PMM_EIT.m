function [sigmak,xc,yc,elec_nodes,movements]=PMM_EIT(vhomo,vinhomo,N,skipcurr,skipvolt,lambda,mu)

%%%%%default beta version for 16 electrodes

pathThoracic='Thoracic_MOM\';
addpath(pathThoracic)

load([pathThoracic 'Jk.mat'])
load([pathThoracic 'RtR_Laplace.mat'])
load([pathThoracic 'Gaussian.mat'])
load([pathThoracic 'xc.mat'])
load([pathThoracic 'yc.mat'])
load([pathThoracic 'elec_nodes.mat'])

LC=size(RtR,1);

%%% set pattern mapping
[stim,els] = mk_stim_patterns(N,1,[0,skipcurr+1],[0,skipvolt+1],{},0.001);
sizes=size(stim(1).meas_pattern);
meas=sizes(1);
elsnew=(reshape(els,N,N))';
voltagel1=zeros(N,meas);
for k=1:N
    space=0;
    for l=1:meas
        while elsnew(k,l+space)==0
            space=space+1;
        end
        voltagel1(k,l)=l+space;
    end
end

voltagel2=zeros(N,meas);

for currmeas=1:N
    for volmeas=1:meas
        if voltagel1(currmeas,volmeas)+skipvolt+1<=N
            voltagel2(currmeas,volmeas)=voltagel1(currmeas,volmeas)+skipvolt+1;
        else
            voltagel2(currmeas,volmeas)=voltagel1(currmeas,volmeas)+skipvolt+1-N;
        end
    end
end

if mu==0
    fprintf('No movement prior\n')
    Q=lambda^2*RtR;
    Mplus=Jk;
else
    delta=1e-03;
    load([pathThoracic 'green_on_elecs1.mat'])
    load([pathThoracic 'green_on_elecs2.mat'])
    mm=1;
    for electrode=1:N
        for voltagechannel=1:meas
            Meas_default(mm)=green_on_elecs1(electrode,voltagechannel);
            Meas_pert(mm)=green_on_elecs2(electrode,voltagechannel);
            mm=mm+1;
        end
    end
    
    for d=1:2
        for electrode=1:N
            Idx=N*(d-1)+electrode;
            Mm(:,Idx)=1/delta*(Meas_pert(:)-Meas_default(:));
        end
    end
    
    n = 2*N;
    e = -ones(n,1);
    Qb= spdiags([e -2.1*e e],-1:1,n,n);
    Q=[lambda^2*RtR zeros(size(RtR,1),size(Qb,2)); zeros(size(Qb,1),size(RtR,2)) mu^2*Qb];
    Mplus=[Jk Mm];
end

for frame=1:size(vinhomo,2)
    %img1(frame)= inv_solve(IMDL, zc_resp(:,1), zc_resp(:,frame));
    % figure
    % show_fem(img(frame))
    c1(frame,:)=-(Mplus'*Mplus+Q)\(Mplus'*(real(vinhomo(:,frame))-real(vhomo)));
    for point=1:LC
        sigmak(frame,point)=exp(sum(c1(frame,1:LC)'.*(g(:,point))));
    end
end

if mu~=0
    movements=c1(:,LC+1:end);
else
    movements=[];
end

end