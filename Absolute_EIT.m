function images=Absolute_EIT(vhomo,vinhomo,N,skipcurr,skipvolt,Q,lambda)

[stimpattern, elspattern] = mk_stim_patterns(N,1,[0 skipcurr+1],[0 skipvolt+1],{'no_meas_current'}, 1);

invmodel=mk_common_model('d2T3',N);
invmodel.fwd_model.stimulation = stimpattern;
invmodel.fwd_model.meas_select=elspattern;
invmodel.fwd_solve.get_all_meas = 1;

Noframes=size(vinhomo,2);

for frame=1:Noframes
    images(frame)=mk_image(invmodel);
end

simdatatest=fwd_solve(images(1));
sref=((simdatatest.meas'*vhomo)/(simdatatest.meas'*simdatatest.meas))^(-1);
invmodel.jacobian_bkgnd.value=sref;
invmodel.elem_data=sref*ones(length(images(1).elem_data),1);
J=calc_jacobian(invmodel);
W=sparse(1:length(vhomo),1:length(vhomo),vhomo.^2);

sigmafirstframe=(J'*W*J+lambda^2*Q)\(J'*W*vhomo);

for frame=1:Noframes
    sigma(:,frame)=(J'*W*J+lambda^2*Q)\(J'*W*vinhomo(:,frame));
    images(frame).elem_data=sigma(:,frame)-sigmafirstframe;
end

end