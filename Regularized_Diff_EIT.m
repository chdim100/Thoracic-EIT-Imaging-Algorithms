function [images,Q,invmodel]=Regularized_Diff_EIT(vhomo,vinhomo,N,skipcurr,skipvolt,algorithm,prior,lambda,tp5D)

[stimpattern, elspattern] = mk_stim_patterns(N,1,[0 skipcurr+1],[0 skipvolt+1],{'no_meas_current'}, 1);

invmodel=mk_common_model('d2T3',N);
%invmodel.fwd_model.nodes=invmodel.fwd_model.nodes/100;
invmodel.fwd_model.stimulation = stimpattern;
invmodel.fwd_model.meas_select=elspattern;
invmodel.fwd_solve.get_all_meas = 1;

switch algorithm
    case{'GN','GN2','G-N','Gauss-Newton'}
        eidors_msg('Difference GN-Recostruction');
        invmodel.solve=@inv_solve_gn;
        invmodel.inv_solve_gn.max_iterations = 20;
    case{'One-step','direct','Direct','Single-Step'}
        eidors_msg('Difference Single-Step Recostruction');
        invmodel.solve=@inv_solve_diff_GN_one_step;
    case{'Total Variation','TV','Total variation','total variation'}
        eidors_msg('Difference TV-Recostruction');
        invmodel.solve=@inv_solve_TV_pdipm;
        invmodel.parameters.term_tolerance= 1e-3;
        invmodel.parameters.keep_iterations= 0;
    otherwise
        error('This algorithm does not exist')
end

switch prior
    case{'Gaussian','gaussian','Gauss','gauss','HPF'}
        invmodel.RtR_prior = @prior_gaussian_HPF;
    case{'Laplace','laplace','LAPLACE'}
        invmodel.RtR_prior = @prior_laplace;
    case{'NOSER','Noser','noser'}
        invmodel.RtR_prior = @prior_noser;
    case{'Tikhonov','Standard Tikhonov','tikhonov','Std'}
        invmodel.RtR_prior = @prior_tikhonov;
    case{'Total Variation','TV','Total variation','total variation'}
        invmodel.R_prior = @prior_TV;
    otherwise
        error('This prior does not exist')
end

invmodel.hyperparameter.value=lambda;
Noframes=size(vinhomo,2);

if tp5D==0
    for frame=1:Noframes
        images(frame)=inv_solve(invmodel,vhomo,vinhomo(:,frame));
    end
else
    LB=size(invmodel.fwd_model.boundary,1);
    boundary=zeros(LB,2);
    
    li=1;
    for ll=1:LB
        boundary(ll,:)=[invmodel.fwd_model.nodes(invmodel.fwd_model.boundary(li,1),1)...
            invmodel.fwd_model.nodes(invmodel.fwd_model.boundary(li,1),2)];
        boundary(ll+1,:)=[invmodel.fwd_model.nodes(invmodel.fwd_model.boundary(li,2),1)...
            invmodel.fwd_model.nodes(invmodel.fwd_model.boundary(li,2),2)];
        li=li+1;
    end
    boundary(LB+1,:)=[];
    boundary=boundary/100;
    invmodel.fwd_model.nodes=invmodel.fwd_model.nodes/100;
    for electrode=1:N
        elec_x(electrode)=mean([invmodel.fwd_model.nodes(invmodel.fwd_model.electrode(electrode).nodes(1),1)...
            invmodel.fwd_model.nodes(invmodel.fwd_model.electrode(electrode).nodes(2),1)]);
        elec_y(electrode)=mean([invmodel.fwd_model.nodes(invmodel.fwd_model.electrode(electrode).nodes(1),2)...
            invmodel.fwd_model.nodes(invmodel.fwd_model.electrode(electrode).nodes(2),2)]);
    end
    elec_pos=[360*(atan2(elec_y,elec_x)/(2*pi))' 0.5*ones(N,1)];
    
    elec_shape = [0.05,0,0.01];
    %
    %invmodel3D=ng_mk_extruded_model({1, boundary, [4 50], 0.06}, elec_pos, elec_shape);
    invmodel3D=ng_mk_extruded_model({1, boundary, [2 12], 0.2}, elec_pos, elec_shape);
    c2f = mk_coarse_fine_mapping(invmodel3D, invmodel.fwd_model);
    invmodel.rec_model = invmodel.fwd_model;
    invmodel.fwd_model=invmodel3D;
    invmodel.fwd_model.stimulation = stimpattern;
    invmodel.fwd_model.meas_select=elspattern;
    invmodel.fwd_solve.get_all_meas = 1;
    invmodel.fwd_model.coarse2fine = c2f;
    %invmodel.fwd_model.jacobian = @jacobian_adjoint_2p5d_1st_order;
    for frame=1:Noframes
        images(frame)=inv_solve(invmodel,vhomo,vinhomo(:,frame));
    end
    invmodel.rec_model.nodes=invmodel.rec_model.nodes*100;
end

Q=calc_RtR_prior(invmodel);

end