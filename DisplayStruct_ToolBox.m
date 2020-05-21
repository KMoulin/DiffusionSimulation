classdef DisplayStruct_ToolBox

    methods(Static)
        
        %% Function that display the Cells inside a voxel 
        function Cells(stack_cell,voxel_dim)
            hold on
            vect=[0 0 voxel_dim(1) voxel_dim(2)];
            rectangle('Position',vect,'Curvature',[0,0], 'FaceColor','w','EdgeColor','k','LineWidth',2)
            for cpt=1:1:size(stack_cell,1) 
                      x1(cpt)=stack_cell(cpt,3)-stack_cell(cpt,1);
                      y1(cpt)=stack_cell(cpt,4)-stack_cell(cpt,1);
                      x2(cpt)=2*stack_cell(cpt,1);
                      y2(cpt)=2*stack_cell(cpt,1);  
                      vect=[x1(cpt) y1(cpt) x2(cpt) y2(cpt)];
                      rectangle('Position',vect,'Curvature',[1,1], 'FaceColor','w','EdgeColor','k','LineWidth',2)
            end
        end

        %% Function that display the Polynomes Cells inside a voxel 
        function Cells_Poly(stack_cell,voxel_dim)
            hold on
            vect=[0 0 voxel_dim(1) voxel_dim(2)];
            rectangle('Position',vect,'Curvature',[0,0], 'FaceColor','w','EdgeColor','k','LineWidth',2)
            for cpt=1:1:size(stack_cell,1)    
                      x=stack_cell(cpt,3)+stack_cell(cpt,9:end).*cos(2*pi*([1:stack_cell(cpt,8)]-1)/stack_cell(cpt,8));
                      x=[x (stack_cell(cpt,3)+stack_cell(cpt,9).*cos(0))];
                      y=stack_cell(cpt,4)+stack_cell(cpt,9:end).*sin(2*pi*([1:stack_cell(cpt,8)]-1)/stack_cell(cpt,8));
                      y=[y (stack_cell(cpt,4)+stack_cell(cpt,9).*sin(0))];
                      plot(x,y,'k','LineWidth',2)
            end
        end

         %% Function that display the Cells inside a voxel 
        function Cells_Mask(mask_cell)
            hold on
            imagesc(mask_cell);
        end
        
        %% Function that add to a gif stack and eventually print it 
        function img_stack=Add_2_Gif(name,bprint,img_stack)
                img=getframe(gcf);
                [tmpp,~]=frame2im(img);
                img_stack(:,:,:,end+1)=tmpp;

                if bprint
                     [tmp_frame,cm] = rgb2ind(squeeze(img_stack(:,:,:,end)),256);
                     for cpt=2:1:size(img_stack,4)
                         [tmp_frame,~] = rgb2ind(squeeze(img_stack(:,:,:,cpt)),cm);
                         if cpt==2
                             imwrite( tmp_frame,cm,[name '.gif'],'gif', 'Loopcount',inf,'DelayTime',0.5);   %%%% First image, delay time = 0.1s         
                         else
                            imwrite( tmp_frame,cm,[name '.gif'],'gif','WriteMode','append','DelayTime',0.5); %%%% Following images
                         end 
                     end
                end
        end

         %% Function that display the cells in 3D
        function Cells_3D(stack_cell,voxel_dim)
            hold on
            vect=[0 0 voxel_dim(1) voxel_dim(2)];
            rectangle('Position',vect,'Curvature',[0,0], 'FaceColor','w','EdgeColor','k','LineWidth',2)
            for cpt=1:1:size(stack_cell,1) 
                      [X,Y,Z] = cylinder([stack_cell(cpt,1) stack_cell(cpt,1)]);   
                      mesh(X+stack_cell(cpt,3),Y+stack_cell(cpt,4),Z*stack_cell(cpt,2)+stack_cell(cpt,5),'facecolor',[1 0 0])
            end
        end
        
         %% Function that display the motion of the water inside the cells
        function Water(Water)
            hold on
            plot(squeeze(Water(:,1,:))',squeeze(Water(:,2,:))')
            scatter(squeeze(Water(:,1,end))',squeeze(Water(:,2,end))',5,'b','filled')
        end
        
         %% Function that display the motion of the water inside the cells
        function Water_Mask(Water,Resol)
            hold on
            plot(squeeze(Water(:,1,:)./Resol(1))',squeeze(Water(:,2,:)./Resol(2))')
            scatter(squeeze(Water(:,1,end)./Resol(1))',squeeze(Water(:,2,end)./Resol(2))',5,'b','filled')
        end
        
        %% Function that display the motion of the water inside the cells
        function Water_3D(Water)
            hold on
            plot3(squeeze(Water(:,1,:))',squeeze(Water(:,2,:))',squeeze(Water(:,3,:))')
            scatter3(squeeze(Water(:,1,end))',squeeze(Water(:,2,end))',squeeze(Water(:,3,end))',5,'b','filled')
        end
    
    end  
end
