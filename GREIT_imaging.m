function [images,imdlgr,invmodel3D]=GREIT_imaging(vhomo,vinhomo,N,skipcurr,skipvolt,s_size)

[stimpattern, elspattern] = mk_stim_patterns(N,1,[0 skipcurr+1],[0 skipvolt+1],{'no_meas_current'}, 1);

invmodel=mk_common_model('d2T3',N);
invmodel.fwd_model.stimulation = stimpattern;
invmodel.fwd_model.meas_select=elspattern;
invmodel.fwd_solve.get_all_meas = 1;

opt.noise_figure = 0.5;

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
%invmodel3D=mk_common_model('j2T3',N);

invmodel3D.stimulation = stimpattern;
invmodel3D.meas_select=elspattern;
invmodel3D.fwd_solve.get_all_meas = 1;

opt.imgsz=s_size;

imdlgr = mk_GREIT_model(invmodel3D,0.2,[],opt);
Noframes=size(vinhomo,2);

for frame=1:Noframes
    images(frame)=inv_solve(imdlgr,vhomo,vinhomo(:,frame));
end

end