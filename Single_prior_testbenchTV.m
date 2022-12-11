N=16;
skipcurr=0;
skipvolt=0;
if skipcurr==skipvolt
    meas=N-3;
else
    meas=N-4;
end
SNR=40;
run C:\Users\Chris\Desktop\PHD_files\eidors-v3.9-ng\eidors\startup.m
%%%set simulation model
[stim,meas_sel] = mk_stim_patterns(N,1,[0,skipcurr+1],[0,skipvolt+1],{'no_meas_current'}, 1);
addpath('C:\Users\Chris\Desktop\PHD_files\eidors-v3.9-ng\eidors\EIT_Circuit_Sim\Structures\Thoracic_New\')
addpath('C:\Users\Chris\Desktop\PHD_files\eidors-v3.9-ng\eidors\EIT_Circuit_Sim\Structures\Thoracic_New\freq_elements')
load('Thorax_Deflated_Asymmetric_Rel=0.03.mat')
Thor1=img;
Thor1.fwd_model.stimulation=stim;
load('Thorax_Inflated_Asymmetric_Rel=0.03.mat')
Thor2=img;
Thor2.fwd_model.stimulation=stim;
load('Thorax2020_Deflated_Asymmetric_15k_Rel=0.03.mat')
Thor1.elem_data=element_data;
load('Thorax2020_Inflated_Asymmetric_15k_Rel=0.03.mat')
Thor2.elem_data=element_data;
vhomo=fwd_solve(Thor1);
vhomonoise=awgn(vhomo.meas,SNR,'measured');
vinhomo=fwd_solve(Thor2);
vinhomonoise=awgn(vinhomo.meas,SNR,'measured');
figure
plot(abs(vhomonoise),'LineWidth',2)
hold on
plot(abs(vinhomonoise),'LineWidth',2)

%%%%% reconstruction model
invmodel=mk_common_model('d2T3',N);
invmodel.fwd_model.stimulation = stim;
invmodel.fwd_model.meas_select=meas_sel;
invmodel.fwd_solve.get_all_meas = 1;
invmodel.reconst_type = 'difference';
invmodel.hyperparameter.value = 1e-6;
invmodel.solve=       @inv_solve_TV_pdipm;
invmodel.R_prior=     @prior_TV;
invmodel.parameters.term_tolerance= 1e-3;
invmodel.parameters.keep_iterations= 0;

ims=inv_solve(invmodel,vhomonoise,vinhomonoise);
figure
show_fem(ims)

L=length(ims.elem_data);
invmodel.jacobian_bkgnd.value=1;
invmodel.elem_data=1*ones(L,1);
simdatatest=fwd_solve(invmodel);
sref=(simdatatest.meas'*vhomonoise/(simdatatest.meas'*simdatatest.meas))^(-1);
invmodel.jacobian_bkgnd.value=sref;
invmodel.elem_data=sref*ones(L,1);
J=calc_jacobian(invmodel);
%%%TV prior
%%%L edges Matrix
invmodel_for_TV=invmodel;
invmodel_for_TV.RtR_prior=@prior_TV;
Ledges=calc_RtR_prior(invmodel_for_TV);
% lambdaTV=calc_hyperparameter(invmodel_for_TV);
lambdaTV=5*1e-03;
betasmooth=1e-03;
%==========================================
data=vinhomonoise-vhomonoise;
nofedges=size(Ledges,1);
K=eye(nofedges,nofedges);
Ek=sqrt(betasmooth)*eye(nofedges,nofedges);
Einv=1/sqrt(betasmooth)*eye(nofedges,nofedges);
Qa=Ledges'*Einv*K*Ledges;
s1=(J'*J+lambdaTV^2*Qa)\(J'*data);
xi=Einv*K*Ledges*s1;
xi=dual_variable_scaling_rule(xi);
invrec1=ims;
invrec1.elem_data=s1;
simdata=fwd_solve(invrec1);

max_iterations=20;
tolerance=1e-03;
% show_fem(invrec1)
factor= 0; norm_d_data= inf;
%%%%
for iter= 1:max_iterations
   betasmooth=betasmooth/1.1;
   d_data= data - simdata.meas;
   prev_norm_d_data= norm_d_data; norm_d_data= norm(d_data);
   eidors_msg('Multiple Prior-Absolute EIT Testbench: iter=%d err=%f factor=%f', ...
       iter,norm_d_data, factor, 2);
   
   if prev_norm_d_data - norm_d_data < tolerance
       break;
   end
   
   %%%%%%calculate new Jacobian
   J=calc_jacobian(invrec1);
   
   %%%Calculate New Priors
   Ek=zeros(nofedges,nofedges);
   K=zeros(nofedges,nofedges);
   for edgei=1:nofedges
       Ek(edgei,edgei)=sqrt((norm(Ledges(edgei,:)*s1))^2+betasmooth);
       K(edgei,edgei)=1-xi(edgei)*Ledges(edgei,:)*s1/Ek(edgei,edgei);
   end
   Einv=diag(diag(Ek).^-1);
   %%%%update prior
   Qa=Ledges'*Einv*K*Ledges;
   Qb=Ledges'*Einv*Ledges;
  
   %%%% updates of dual variables
   Fchange=left_divide((J'*J+lambdaTV^2*Qa),...
       (J'*d_data-lambdaTV^2*Qb*s1));
   xi=Einv*Ledges*s1+Einv*K*Ledges*Fchange;
   xi=dual_variable_scaling_rule(xi);
   %%%linesearch
   facts= [0, logspace(-3,0,19)];  
   norms= norm_d_data;
   for f= 2:length(facts)
      snew= s1 + facts(f)*Fchange;
      invrec1.elem_data=snew;
      simdata=fwd_solve(invrec1);
      norms(f)= norm(data - simdata.meas);
   end
   ff= find(norms==min(norms));
   factor= facts(ff(end));
   s1 = s1 + factor*Fchange;
   invrec1.elem_data= s1;
   simdata=fwd_solve(invrec1);
end
figure
show_fem(invrec1)