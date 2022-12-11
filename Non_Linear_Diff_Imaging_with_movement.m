function [images,movements]=Non_Linear_Diff_Imaging_with_movement(vhomo,vinhomo,N,L,LROI,Elems_ROI,skipcurr,skipvolt,Qsmooth,lambda,lambdaTV,lambdamove,beta)

[stimpattern, elspattern] = mk_stim_patterns(N,1,[0 skipcurr+1],[0 skipvolt+1],{'no_meas_current'}, 1);

if skipcurr==skipvolt
    meas=N-3;
else
    meas=N-4;
end
Noframes=size(vinhomo,2);

invmodel=mk_common_model('d2T3',N);
invmodel.fwd_model.stimulation = stimpattern;
invmodel.fwd_model.meas_select=elspattern;
invmodel.fwd_solve.get_all_meas = 1;

simdatatest=fwd_solve(mk_image(invmodel,1));
sref=((simdatatest.meas'*vhomo)/(simdatatest.meas'*simdatatest.meas))^(-1);
invrec1(1)=mk_image(invmodel,1);
for frame=1:Noframes
    invrec1(frame)=invrec1(1);
    invrec1(frame).jacobian_bkgnd.value=sref;
    invrec1(frame).elem_data=sref*ones(L,1);
    invrec2(frame)=invrec1(frame);
end

%%%%Calculate Edge Matrix
invmodel_for_TV=invmodel;
invmodel_for_TV.RtR_prior=@prior_TV;
Ledges=calc_RtR_prior(invmodel_for_TV);
Ledges_ROI=Ledges(:,Elems_ROI);
Ledges_ROI(~any(Ledges_ROI,2),:) = [];
Ledges_ROI(sum(Ledges_ROI,2)~=0,:)=[];
nofedges=size(Ledges_ROI,1);

%%%%Calculate movement Prior
%lambdamove=sqrt(1e3/1); 
n = 2*N;
e = -ones(n,1);
 Qmove= spdiags([e -2.1*e e],-1:1,n,n);

%Qmove= speye(n);

W=sparse(1:length(vhomo),1:length(vhomo),vhomo.^2);
Wn=[W zeros(N*meas,N*meas); zeros(N*meas,N*meas) W];
eidors_msg('Multiple Prior-ROI Reconstruction');
threshold=0.01;
max_iters=20;
b=beta;
for frame=1:Noframes
    eidors_msg('Frame  %d',frame,2);
    data=[vhomo; vinhomo(:,frame)];
    residual_change=inf;
    residual=1e+06;
    iteration=1;
    s1=sref*ones(L,1);  deltasigmaROIs=zeros(LROI,1);
    sbar=[s1; zeros(2*N,1); deltasigmaROIs; zeros(2*N,1)];
    srefo=sbar;
    Udata=zeros(2*N*meas,1);
    d_data=data;
    norm_d_data= inf;
    beta=b;
    %%% initialize dual variable xi
    xi=zeros(size(Ledges_ROI,1),1);
    
    while residual_change>threshold&&iteration<=max_iters
        residual_prev=residual;
        %%%Update Jacobi
        Jac=calc_jacobian(invrec1(frame));
        Jam=movement_jacobian_only(N*meas,N,2,1e-06,invrec1(frame));
        %Jb=zeros(N*meas,LROI);
        Jb=zeros(N*meas,LROI+2*N);
        %Jc=calc_jacobian(invrec2(frame));
        Jcc=calc_jacobian(invrec2(frame));
        Jcm=movement_jacobian_only(N*meas,N,2,1e-06,invrec2(frame));
        Jdc=Jac(:,Elems_ROI);
        Jdm=Jcm;
        Ja=[Jac Jam];
        %Jb=[Jbc Jbm];
        Jc=[Jcc Jcm];
        Jd=[Jdc Jdm];
        Jglob=[Ja Jb; Jc Jd];
        %%%
        %%%construct Qsteep and Qglobe in TV case
        Ek=zeros(nofedges,nofedges);
        K=zeros(nofedges,nofedges);
        for edgei=1:nofedges
            Ek(edgei,edgei)=sqrt((norm(Ledges_ROI(edgei,:)*deltasigmaROIs))^2+beta);
            K(edgei,edgei)=1-xi(edgei)*Ledges_ROI(edgei,:)*deltasigmaROIs/Ek(edgei,edgei);
        end
        Einv=diag(diag(Ek).^-1);
        %%%%update priors
        Qsteepa=Ledges_ROI'*Einv*K*Ledges_ROI;
        Qsteepb=Ledges_ROI'*Einv*Ledges_ROI;
        Qglobea=[lambda^2*Qsmooth zeros(size(Qsmooth,1),2*N)...
            zeros(size(Qsmooth,1),size(Qsteepa,2)) zeros(size(Qsmooth,1),2*N);...
            zeros(2*N,size(Qsmooth,2)) lambdamove^2*Qmove zeros(2*N,size(Qsteepa,2)+2*N);...
            zeros(size(Qsteepa,1),size(Qsmooth,2)+2*N) lambdaTV^2*Qsteepa zeros(size(Qsteepa,1),2*N);...
            zeros(2*N,size(Qsmooth,2)+2*N+size(Qsteepa,2)) lambdamove^2*Qmove];
        
       Qglobeb=[lambda^2*Qsmooth zeros(size(Qsmooth,1),2*N)...
            zeros(size(Qsmooth,1),size(Qsteepb,2)) zeros(size(Qsmooth,1),2*N);...
            zeros(2*N,size(Qsmooth,2)) lambdamove^2*Qmove zeros(2*N,size(Qsteepb,2)+2*N);...
            zeros(size(Qsteepb,1),size(Qsmooth,2)+2*N) lambdaTV^2*Qsteepb zeros(size(Qsteepb,1),2*N);...
            zeros(2*N,size(Qsmooth,2)+2*N+size(Qsteepb,2)) lambdamove^2*Qmove];
        
%         Qglobea=[lambda^2*Qsmooth zeros(size(Qsmooth,1),size(Qsteepa,2));...
%             zeros(size(Qsteepa,1),size(Qsmooth,2)) lambdaTV^2*Qsteepa];
%         Qglobeb=[lambda^2*Qsmooth zeros(size(Qsmooth,1),size(Qsteepb,2));...
%             zeros(size(Qsteepb,1),size(Qsmooth,2)) lambdaTV^2*Qsteepb];
        %%%%calculate change
        ds=left_divide((Jglob'*Wn*Jglob+Qglobea),...
            (Jglob'*Wn*d_data+Qglobeb*(0*srefo-sbar)));
        xi=Einv*Ledges_ROI*deltasigmaROIs+Einv*K*Ledges_ROI*ds(L+2*N+1:L+2*N+LROI);
        xi=dual_variable_scaling_rule(xi);
        beta=beta/1.5;
        %%%%linesearch
        facts= [0, logspace(-3,0,19)];
        norms= norm_d_data;
        for f= 2:length(facts)
            snew= sbar + facts(f)*ds;
            invrec1(frame).elem_data=snew(1:L);
            simdata1=fwd_solve(invrec1(frame));
            invrec2(frame).elem_data=invrec1(frame).elem_data;
            invrec2(frame).elem_data(Elems_ROI)=...
                invrec2(frame).elem_data(Elems_ROI)+snew(L+2*N+1:L+2*N+LROI);
            simdata2=fwd_solve(invrec2(frame));
            simdata=[simdata1.meas; simdata2.meas];
            norms(f)= norm(data - simdata);
        end
        ff= find(norms==min(norms));
        factor= facts(ff(end));
        sbar = sbar + factor*ds;
        s1=sbar(1:L);
        deltasigmaROIs=sbar(L+2*N+1:L+2*N+LROI);
        movements1=sbar(L+1:L+2*N);
        movements2=sbar(L+2*N+1+LROI:end);
        invrec1(frame).elem_data= s1;
        invrec2(frame).elem_data= invrec1(frame).elem_data;
        invrec2(frame).elem_data(Elems_ROI)=...
            invrec2(frame).elem_data(Elems_ROI)+deltasigmaROIs;
        simdata1=fwd_solve(invrec1(frame));
        simdata2=fwd_solve(invrec2(frame));
        Udata(1:N*meas)=simdata1.meas;
        Udata(N*meas+1:end)=simdata2.meas;
        d_data=data-Udata;
        norm_d_data=norm(d_data);
        residual=norm_d_data^2;
        eidors_msg('Multiple Prior-Absolute EIT Testbench: iter=%d err=%f factor=%f', ...
            iteration,norm_d_data, factor, 2);
        residual_change=abs(residual-residual_prev)/residual_prev;
        iteration=iteration+1;
    end
    movements(:,frame)=movements2-movements1;
    images(frame)=invrec1(frame);
    images(frame).elem_data=...
        invrec2(frame).elem_data-invrec1(frame).elem_data;
end

end