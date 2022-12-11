%%%%%%% Performs EM-Based Block Sparse Bayesian Learning combined with a
%%%%%%% Radial Basis Function Based Method-of-Moment

%%%%%%Please Cite:
%%%% Dimas, C., Alimisis, V., Uzunoglu, N., & Sotiriadis, P. P. (2021).
%%%% A point-matching method of moment with sparse Bayesian learning applied
%%%% and evaluated in dynamic lung electrical impedance tomography.
%%%% Bioengineering, 8(12), 191.

%%% Dimas, C., Uzunoglu, N., & Sotiriadis, P. P. (2021).
%%% An efficient Point-Matching Method-of-Moments for 2D and 3D Electrical
%%% Impedance Tomography Using Radial Basis functions.
%%% IEEE Transactions on Biomedical Engineering, 69(2), 783-794.

%%% and

%%%%% Liu, S., Jia, J., Zhang, Y. D., & Yang, Y. (2018). Image reconstruction
%%%%% in electrical impedance tomography based on structure-aware sparse
%%%%% Bayesian learning. IEEE transactions on medical imaging, 37(9), 2090-2102.

close all
clc
clear So Bi Bitilde gi
load('Thoracic_MOM\Jk.mat')
J=Jk;
cavity=1;
if cavity>1&&~strcmp(cavity,'in-vivo')
    %%%%%% works in original code, memory issues in Github!
    error('This demo works only for the 1st cavity! Please select cavity=1 or in-vivo')
end

Reg=[30 30 30 30];
%Reg=[1 1 1 1];
%Reg=[10 10 10 10];
if isnumeric(cavity)==1
    load(['Thoracic_MOM\vhomonoise' num2str(cavity) '.mat'])
    load(['Thoracic_MOM\vinhomonoise' num2str(cavity) '.mat'])
elseif strcmp(cavity,'in-vivo')
    load('Thoracic_MOM\vinvivo.mat')
    %%%%% scale (naive way, trial and error based)
    vhomo=real(vhomo/100);
    vinhomo=real(vinhomo/100);
end
for imm=1:size(vinhomo,2)
    %%%%% get difference measurements, m-1 frames from m states
    v=-vinhomo(:,imm)+vhomo;
    %%%%% Initialize
    M=size(J,1);
    N=size(J,2);
    beta=0;
    emin=1e-06;
    epsilon=1;
    thetamax=10;
    h=4;
    g=N-h+1;
    if isnumeric(cavity==1)
        Regs=Reg(cavity);
    else
        Regs=5;
    end
    Psi=zeros(N,h*g);
    Sx=zeros(g*h,g*h);
    theta=0;
    mux=zeros(g*h,1);
    gi=ones(g,1);
    go=0.01*sqrt(1/(N-1)*sum(abs(v-mean(v)).^2));
    Bi=zeros(g,h,h);
    Biblock=toeplitz(0.9.^(0:h-1));
    Bitilde=zeros(g,h,h);
    Binew=Bi;
    %%%%compute Psi, Phi, Bi (initial), So
    for i=1:g
        Bi(i,:,:)=Biblock;
        Psi(:,(i-1)*h+1:i*h)=sparse([zeros(i-1,h)' eye(h,h)' zeros(N-i-h+1,h)']');
        So(h*(i-1)+1:i*h,:)=sparse([zeros(size(Biblock,1),(i-1)*h) gi(i)*Biblock zeros(size(Biblock,1),h*(g-1)-(i-1)*h)]);
    end
    Phi=J*Psi;
    Im=eye(M,M);
    %%%%Initialize Bitilde
    tic
    %%%%iterations and updates
    while epsilon>emin && theta<=thetamax
        PhiSoPhiT=Phi*So*Phi';
        %%%%update mux using (15)
        Sv=go*Im+PhiSoPhiT;
        muxprev=mux;
        mux=So*Phi'/Sv*v;
        %%%%update Sx using (16)
        Sx=So-So*Phi'*(Sv\Phi*So);
        %%%%update gi using (22)
        summary=0;
        for i=1:g
            Phi_i=Phi(:,(i-1)*h+1:i*h);
            Phi_iT=Phi_i';
            Sxi=Sx((i-1)*h+1:i*h,(i-1)*h+1:i*h);
            Biblock=squeeze(Bi(i,:,:));
            muxi=mux((i-1)*h+1:i*h);
            muxiT=mux((i-1)*h+1:i*h)';
            gi(i)=norm(sqrt(Biblock)*Phi_iT*(Sv\v))/...
                sqrt(trace(Phi_iT*(Sv\(Phi_i*Biblock))));
            summary=summary+...
                trace(Sxi*Phi_iT*Phi_i);
            %%%%update Bi using (18)-(21)
            %%%%update Bitilde (18)
            Bitildei=squeeze(Bitilde(i,:,:));
            Bitilde(i,:,:)=Bitildei+1/gi(i)*...
                (Sxi+muxi*muxiT);
            Bitildei=squeeze(Bitilde(i,:,:));
            %%%compute rtildei (21)
            rtildei=mean(diag(Bitildei,1))/mean(diag(Bitildei));
            %%%compute ri (20)
            ri=sign(rtildei)*min([norm(rtildei) 0.99]);
            Biblocknew=toeplitz(ri.^(0:1:h-1));
            Binew(i,:,:)=Biblocknew;
            %%%update So
            So(h*(i-1)+1:i*h,:)=sparse([zeros(size(Biblocknew,1),(i-1)*h) gi(i)*Biblocknew zeros(size(Biblocknew,1),h*(g-1)-(i-1)*h)]);
        end
        %%%%update go using (23)
        go=1/M*(norm(v-Phi*mux)^2+summary/Regs);
        Bi=Binew;
        epsilon=norm(mux-muxprev)/norm(mux);
        theta=theta+1;
    end
    toc
    %%%%output
    load('Thoracic_MOM\Gaussian.mat')
    cc=Psi*mux/100;
    %cc=Psi*mux/100;
    sigma=cc;
    for point=1:size(Jk,2)
        sigma(point)=exp(sum(cc.*(g(:,point))));
    end
    
    load('Thoracic_MOM\xc.mat')
    load('Thoracic_MOM\yc.mat')
    
    figure
    scatter3(xc,yc,real(sigma),140,real(sigma),'filled','s')
    axis tight
    axis off
    colormap jet
    images_MoM_SBL(imm,:)=sigma;
    if isnumeric(cavity)==1
        load(['invreference' num2str(cavity) '_MoM.mat'])
        RR2=corrcoef(real(sigma'),...
            real(inv_reference_MoM(imm,:)));
        CC_MoM2(imm)=RR2(1,2);
        caxis([0.93 1.07])
        view([0 90])
    else
        caxis([0.98 1.02])
        view([0 90])
    end
    title(['MoM-BSBL-Frame: ', num2str(imm)])
    fprintf('Frame %2.0f/%2.0f Completed\n',imm,size(vinhomo,2))
end
if isnumeric(cavity)==1
    CC_MoM2
    mean(CC_MoM2)
end