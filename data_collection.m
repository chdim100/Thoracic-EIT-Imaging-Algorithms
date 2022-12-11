function [vhomonoise,vinhomonoise,sim_model]=data_collection(N,skipcurr,skipvolt,SNR,type,cavity)
%%%%%%%path to EIDORS
Eidorspath='C:\Users\Chris\Desktop\PHD_files\'; %%%%%Set path to eidors!!
startuped=0;
if ~exist('startuped')||startuped==0
    try
        run ([Eidorspath,'eidors-v3.9-ng\eidors\startup.m'])
        startuped=1;
    catch
        error('This code needs the EIDORS package to run properly, please download from https://eidors3d.sourceforge.net/download.shtml')
    end
    
end
%%% set pattern mapping
[stimpattern, elspattern] = mk_stim_patterns(N,1,[0 skipcurr+1],[0 skipvolt+1],{'no_meas_current'}, 1);

switch type
    case{'Thorax3D'}
        
        addpath('Thoracic_Basis\')
        addpath('Thoracic_Basis\freq_elements')
        addpath(['Thoracic_Cavities\Thorax' num2str(cavity) '\'])
        load(['thorax', num2str(cavity) '_3Dsequence.mat'])
        for im=1:length(thorax_model)
            %%%%configure stimulation pattern
            thorax_model{im}.fwd_model.stimulation=stimpattern;
            %%%%reconfigure electrode sequence (according to DICOM
            %%%%orientation)
            eltemp1=thorax_model{im}.fwd_model.electrode;
            thorax_model{im}.fwd_model.electrode(1:9)=eltemp1(9:-1:1);
            thorax_model{im}.fwd_model.electrode(10:16)=eltemp1(16:-1:10);
        end
        
        %%%get measurements
        
        vhomo=fwd_solve(thorax_model{1});
        vhomonoise=awgn(vhomo.meas,SNR,'measured');
        for im=2:length(thorax_model)
            vinhomo=fwd_solve(thorax_model{im});
            vinhomonoise(:,im-1)=awgn(vinhomo.meas,SNR,'measured');
            sim_model=thorax_model;
        end
        
    case{'In-Vivo','In-Vivo 1995','Guardo','Guardo1995'}
        
        load montreal_data_1995
        Measurements=zc_resp(:);
        Measurements=double(Measurements)/1000;
        vhomonoise=(Measurements(1:N^2));
        vhomonoise=vhomonoise(elspattern);
        Noframes=size(zc_resp,2);
        for frame=2:Noframes
            vinhomonoise_frame=Measurements((frame-1)*N^2+1:N^2*frame);
            vinhomonoise(:,frame-1)=(vinhomonoise_frame(elspattern))';
        end
        
        sim_model='In-Vivo Data';
        
    otherwise
        error('Data structure non-specified')
end

end