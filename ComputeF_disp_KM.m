% function [Eff, Err, Ecc, Ell, Jac, V_t, V_t_nodes] = ComputeF_Exp(phi_all, phi_up, phi_down, query_up, query_down, dN_all, QPmid_im, cc, rr, ll, fail_index, f, phase)
function [nodes2] = ComputeF_disp_KM(nodes_dense, nodes, ref_phase,nphase)



phi_up=nodes_dense(1).points;
phi_mid=nodes_dense(2).points;
phi_down=nodes_dense(3).points;
phi_all = [nodes_dense(1).points;nodes_dense(3).points];



nodes2=nodes;
% Select the DTI points we are going to move.  
QDTI_Points=[];
QROI_Points_endo=[];
QROI_Points_epi=[];
ref_phase_DTI=ref_phase;

% If there is only one cardiac phase for DTI
if size(nodes,3)<ref_phase
    ref_phase_DTI=1;
end

QDTI_Points=      nodes(:,1:3,ref_phase_DTI);


cpt_t1=1;
for cpt_t = nphase
    %DENSE
    
     % Intial config at the given phase from DTI (reference config)
    X_radius = [phi_up(:,:,ref_phase);phi_mid(:,:,ref_phase);phi_down(:,:,ref_phase)];
        
    % Position throught time
    x_radius = [phi_up(:,:,cpt_t);phi_mid(:,:,cpt_t);phi_down(:,:,cpt_t)];   
    
    % Interpolant of position in X Y Z as a regard of initial config. 
    Phi_x = scatteredInterpolant(X_radius(:,1), X_radius(:,2), X_radius(:,3), x_radius(:,1),'linear','linear');
    Phi_y = scatteredInterpolant(X_radius(:,1), X_radius(:,2), X_radius(:,3), x_radius(:,2),'linear','linear');
    Phi_z = scatteredInterpolant(X_radius(:,1), X_radius(:,2), X_radius(:,3), x_radius(:,3),'linear','linear');

    % Apply deformation for each dTI point
    for cpt_p = 1:size(QDTI_Points,1)
        
        % Current DTI point
        DTI_point = squeeze(QDTI_Points(cpt_p,:));
        
        % Local definition of the orientation vectors

        %Fscat = Compute_F_local (DTI_point,Phi_x,Phi_y,Phi_z);
        %Fscat_save(:,:,cpt_p,cpt_t1)=Fscat;
        
        %% Apply the deformation to the vectors

        
        nodes2(cpt_p,:,cpt_t1) = [Phi_x(DTI_point); Phi_y(DTI_point); Phi_z(DTI_point)];
        
      
    end
    cpt_t1=cpt_t1+1;
     
end

end


function F = Compute_F_local (point,Phi_x,Phi_y,Phi_z)
       %Defining Fscat
        Delta = 0.1; % was in pixel should be in mm (0.1 originally)
       % 11 21 31
        Xq_plus  = point; 
        Xq_plus(1)  = Xq_plus(1)  + Delta;
        Xq_minus = point; 
        Xq_minus(1) = Xq_minus(1) - Delta;

        % J is the colum coordinate in the reference system in which we are
        % taking the derivative
        % i is the line the deformation mapping conponement we are
        % differentiate. 
        
        F(1,1) = (Phi_x(Xq_plus) - Phi_x(Xq_minus))/(2*Delta);
        F(2,1) = (Phi_y(Xq_plus) - Phi_y(Xq_minus))/(2*Delta);
        F(3,1) = (Phi_z(Xq_plus) - Phi_z(Xq_minus))/(2*Delta);

        % 12 22 32
        Xq_plus  = point; 
        Xq_plus(2)  = Xq_plus(2)  + Delta;
        Xq_minus = point; 
        Xq_minus(2) = Xq_minus(2) - Delta;

        F(1,2) = (Phi_x(Xq_plus) - Phi_x(Xq_minus))/(2*Delta);
        F(2,2) = (Phi_y(Xq_plus) - Phi_y(Xq_minus))/(2*Delta);
        F(3,2) = (Phi_z(Xq_plus) - Phi_z(Xq_minus))/(2*Delta);

        % 13 23 33
        Xq_plus  = point; 
        Xq_plus(3)  = Xq_plus(3)  + Delta;
        Xq_minus = point; 
        Xq_minus(3) = Xq_minus(3) - Delta;

        % F for Phi
        F(1,3) = (Phi_x(Xq_plus) - Phi_x(Xq_minus))/(2*Delta);
        F(2,3) = (Phi_y(Xq_plus) - Phi_y(Xq_minus))/(2*Delta);
        F(3,3) = (Phi_z(Xq_plus) - Phi_z(Xq_minus))/(2*Delta);


end

function WD= Compute_WD_local (Points,P_Endo,P_Epi,Rot)

    Endo_Line = zeros(size(P_Endo,1),1);
    Epi_Line  = ones(size(P_Epi,1),1);
    PosRoi    = cat(1,P_Endo,P_Epi);
    LineRoi   = cat(1,Endo_Line,Epi_Line);
    
    
    %% Approximate in 2D 
    PosRoi_Rot=(inv(Rot)*(PosRoi)')';
    Points_Rot=(inv(Rot)*(Points)')';
    % griddata(X,Y,V, xq, yq)
    WD=griddata(PosRoi_Rot(:,1),PosRoi_Rot(:,2),LineRoi,Points_Rot(:,1),Points_Rot(:,2));
    
%     Idx_exist=find(~isnan(WD));
%     Idx_nan=find(isnan(WD));
%     Scatter_WD = scatteredInterpolant(Points(Idx_exist,1), Points(Idx_exist,2), Points(Idx_exist,3), WD(Idx_exist),'linear','linear');
%     WD(Idx_nan)=Scatter_WD(Points(Idx_nan,1), Points(Idx_nan,2), Points(Idx_nan,3));
    
end