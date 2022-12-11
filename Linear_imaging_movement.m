function [images,Q,invmodel]=Linear_imaging_movement(vhomo,vinhomo,N,skipcurr,skipvolt,lambda_smooth,mu)
[stimpattern, elspattern] = mk_stim_patterns(N,1,[0 skipcurr+1],[0 skipvolt+1],{'no_meas_current'}, 1);

invmodel=mk_common_model('d2T3',N);
invmodel.fwd_model.stimulation = stimpattern;
invmodel.fwd_model.meas_select=elspattern;
invmodel.fwd_solve.get_all_meas = 1;

invmodel.hyperparameter.value=lambda_smooth;

invmodel.fwd_model.jacobian = @jacobian_movement;
invmodel.RtR_prior =          @prior_movement;
invmodel.prior_movement.parameters = sqrt(mu);

Noframes=size(vinhomo,2);

for frame=1:Noframes
    images(frame)=inv_solve(invmodel,vhomo,vinhomo(:,frame));
end

Q=calc_RtR_prior(invmodel);


end