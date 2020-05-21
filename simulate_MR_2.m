function [Att_approx Phi_Dist]= simulate_MR_2(Water,Seq,Grad,World)
%%% Simulation gradients
Phi_Dist=[];
Nb_Water=size(Water,1);

%[Grad m1]= diff_grad_gen([8000 0 0 4000 8000 0 0],[1 0 0 0 -1 0 0], 350, 500,1) ;              % mT/m (us) & T.us² /m 
    %%% Loop over Direction %%%
    for cpt_dir=1:1:size(Seq.Dir,1)
        %%% Loop over Waveform %%%
        for cpt_g=1:1:size(Grad,2)          
             %%% Intregral Diffusion + gradients  %%%           
             phi=zeros(Nb_Water,1);   
             for cpt=2:1:size(Grad,1)
                phi(:)=phi(:)+ nansum( World.Gamma *World.dT*( repmat(Grad(cpt,cpt_g),Nb_Water,3).*repmat(Seq.Dir(cpt_dir,:),Nb_Water,1).*squeeze(Water(:,:,cpt))),2);
             end

            %%% Phase Distribution
            Phi_Dist(cpt_g,cpt_dir,:)=phi(:);
            %%% MR Signal attenuation
            Att_approx(cpt_g,cpt_dir)=nanmean(cos(phi(:))+i*sin(phi(:)));  
            %%%
        end
    end
end