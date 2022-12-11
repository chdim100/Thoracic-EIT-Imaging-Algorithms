function J= movement_jacobian_only( Nh,N,DE, delta, img )
J = zeros( Nh, N*DE );

node0= img.fwd_model.nodes;
d0= fwd_solve( img );
for d= 1:DE
    for i= 1:N
        idx= img.fwd_model.electrode(i).nodes;
        
        img.fwd_model.nodes( idx, d)= node0(idx,d) + delta;
        di= fwd_solve( img );
        img.fwd_model.nodes( idx, d)= node0(idx,d);
        
        J_idx =N*(d-1) + i;
        J(:,J_idx) = (1/delta) * (di.meas - d0.meas);
    end
end