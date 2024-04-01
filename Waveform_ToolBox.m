classdef Waveform_ToolBox

    methods(Static)
        
        
        
        function [F]= Fcalc(Dn,Handle,dT,Gmax)
            % G in mT/m 
            % dT in us
            % This give a b-value  in s/mm^2
            equpu =@(d) d^3/3;
            
            
            
            equpu0=@(d) d^2;
            
            equv =@(d) 0;
            

            equpu(1000e-3)+equpu0(1000e-3)*(1000e-3)+equpu(1000e-3)
            
            
            bfact=[];
            
            GammaH= 267.513e6;  % Hz/T            
            Gmax=Gmax*1e-3*1e-3;      % mT/m -> T/mm
            
            Handle = {@equh1, @equh2, @equh3};
            
            
            
            dT=dT*1e-6;
            Tn=cumsum(Dn);
            F=[];
            for cpt=1:1:length(Dn)
               Tn=0;
                if cpt==1
                   Tn=linspace(0,Dn(cpt),Dn(cpt)/dT) 
                else
                    Tn=linspace(Dn(cpt-1),Dn(cpt),(Dn(cpt)-Dn(cpt-1))/dT)
                end
                F=[F Handle{cpt}(Tn)];
                
            end
        
        
        end
        function  [bval, M1, Grad, Slope]= AddBvalEvent(bval, M1, Grad, Slope,Sn,TimeStep)

            
            
        dTimeStep  = TimeStep;   
        dTimeStep2 = dTimeStep  * dTimeStep;
        dTimeStep3 = dTimeStep2 * dTimeStep;
        dTimeStep4 = dTimeStep2 * dTimeStep2;
        dTimeStep5 = dTimeStep2 * dTimeStep3;
       
       
        
            
         % Calculate change of b-value
         % (unit: (mT/m)^2 us^3)
          bval   =     bval +      M1     * M1                                     * dTimeStep    ...
                      + 1./2.   * (M1     * Grad  + M1 * Grad)                     * dTimeStep2   ...
                      + 1./6.   * (M1     * Slope + M1 * Slope + 2. * Grad * Grad) * dTimeStep3   ...
                      + 1./8.   * (Grad  * Slope + Grad * Slope)                   * dTimeStep4   ...
                      + 1./20.  * (Slope * Slope)                                  * dTimeStep5;

            % Calculate change of 0th moment (time-integral of gradient)
            % (unit: mT/m us)
             M1 = M1 + Grad * dTimeStep + 1./2. * Slope * dTimeStep2;

            % Calculate change of gradient amplitude 
            % (unit: mT/m)
            Grad = Grad + Slope * dTimeStep;  
            Slope=Slope+Sn;
           
        end
        function [bval] = Biter(Dn,Gn,Rt)
               tmp_time=0;
               
               bval=0;
               M1=0;
               Grad=0;
               Slope=0;
               cG=0;
                for cpt=1:1:length(Dn)   
                       if Dn(cpt)>0 && Dn(cpt)>2*Rt && Gn(cpt)~=0   % Trapezoid
                           
                           % Ramp Up
                           [bval, M1, Grad, Slope]= Waveform_ToolBox.AddBvalEvent(bval, M1, Grad, Slope,Gn(cpt)*1/Rt,0);
                           cG=[cG Grad];
                           % Flatop
                           [bval, M1, Grad, Slope]= Waveform_ToolBox.AddBvalEvent(bval, M1, Grad, Slope,-Gn(cpt)*1/Rt,Rt);
                           cG=[cG Grad];
                           % Ramp Dowm
                           [bval, M1, Grad, Slope]= Waveform_ToolBox.AddBvalEvent(bval, M1, Grad, Slope,-Gn(cpt)*1/Rt,Dn(cpt)-2*Rt);
                           cG=[cG Grad];
                           % zero
                           [bval, M1, Grad, Slope]= Waveform_ToolBox.AddBvalEvent(bval, M1, Grad, Slope,Gn(cpt)*1/Rt,Rt);
                           cG=[cG Grad];
                           
                            tmp_time=tmp_time+Dn(cpt);
                       elseif Dn(cpt)>0 && Dn(cpt)<=2*Rt  && Gn(cpt)~=0  % Triangle
                            GM(cpt).StartTime=tmp_time;
                            GM(cpt).RampTime=Dn(cpt)/2;
                            GM(cpt).FlapTop=0;
                            GM(cpt).StopTime=GM(cpt).StartTime+GM(cpt).FlapTop+2*GM(cpt).RampTime;
                            GM(cpt).EffectiveTime=GM(cpt).FlapTop+GM(cpt).RampTime;
                            GM(cpt).Ampl=Gn(cpt)*GM(cpt).RampTime/Rt;
                            tmp_time=GM(cpt).StopTime; 
                       else                             % Blank
                           % zero
                           [bval, M1, Grad, Slope]= Waveform_ToolBox.AddBvalEvent(bval, M1, Grad, Slope,0,Dn(cpt));
                           cG=[cG Grad];
                           tmp_time=tmp_time+Dn(cpt);
                       end
                end
                GammaH= 267.513e6;  % Hz/T 
                bval=GammaH*GammaH*bval;
        end
        function [bval, bfact]= Bcalc(G,dT)
            % G in mT/m 
            % dT in us
            % This give a b-value in  in s/mm^2
            
            bfact=[];
            
            GammaH= 267.513e6;  % Hz/T            
            G=G*1e-3*1e-3;      % mT/m -> T/mm
            dT=dT*1e-6;         % us -> s
%             btmp=0;
%             bfunc=zeros(1,length(G));
%             
%             for cpt=1:1:length(G)
%                      bfunc(cpt)=btmp+G(cpt)*G(cpt)*dT;
%                      btmp=bfunc(cpt);
%             end
% 
%             TimeFact=(GammaH*GammaH * dT * dT ); %1e-18 * 1e-6
%             bfact = bfunc*TimeFact;
%             bval = btmp*TimeFact;
%             
            Integration1 = cumsum(G.*dT);     % First integral
            bfact=Integration1.*Integration1; % Square of the First integral
            bval = (GammaH^2)*sum(bfact* dT); % Second integral
            bfact=(GammaH^2)*cumsum(bfact*dT);
        end
         
        
         function [bval]= BValue_Mono_ref(delta,DELTA,Gmax,rT)
            %
            % Calculate reference b-value for a monopolar waveform
            % delta (Ramptime +Plateau)in s
            % DELTA in S
            % Gmax could be a 1D vector waveform or an amplitude both in mT/m 
            % ramptime rT in s
            % This give a b-value in s/mm^2
            
            
            Gmax=Gmax*1e-3*1e-3; % Convert mT/m to T/mm
            GammaH= 267.513e6; %% Hz/T

            bval = (GammaH^2)* (max(Gmax)^2) * ( delta* delta*(DELTA - delta/3) + rT*rT*rT/30 - delta *rT*rT/6); % s/mm^2
            
         end
        
        function [M]= M0(G,dT)                  
            t=(0:1:(length(G)-1))*dT;
            M=cumsum(G)*dT;

        end
        
        function [M]= M1(G,dT)
            t=(0:1:(length(G)-1))*dT;
            M=cumsum(G.*t')*dT;
            
        end
        
        function [M]= M2(G,dT)      
            t=(0:1:(length(G)-1))*dT;
            M=cumsum(G.*t'.*t')*dT;
            
        end
        
        function [M]= M3(G,dT)
            t=(0:1:(length(G)-1))*dT;
            M=cumsum(G.*t'.*t'.*t')*dT;
            
        end
        
        function [M]= MomentN(G,dT,N)
            t=(0:1:(length(G)-1))*dT;
            M=cumsum(G.*(t.^N)')*dT;
        end
        
        function [M0,M1,M2,M3]= Moment(G,dT)
            
           M0=Waveform_ToolBox.M0(G,dT);
           M1=Waveform_ToolBox.M1(G,dT);
           M2=Waveform_ToolBox.M2(G,dT);
           M3=Waveform_ToolBox.M3(G,dT);
            
        end
        
        
       
        
%       function [M0]=M0iter( G,  G2, t0,  t1, mode)
% 
%         t=t1-t0;
%         if(mode==1) %'RampUp'
% 
%             M0=(G2*t/2);
% 
%         elseif(mode==2) %'FlapTop'
% 
%             M0=(G*t);
%             
%         elseif(mode==3) %'RampDown'
%             
%             M0= G*t/2;
%             
%         elseif(mode==4) %'BridGe'
%             
%             M0= t*(G+G2)/2;
%             
%         else
%             
%             M0= 0.0;
%             
%         end
%       end
%         
%         function [M1]=M1iter( G,  G2, t0,  t1, mode, scale)
% 
%              t=t1-t0;
%             if(mode==1) %'RampUp'
% 
%                 M1=(G2*t*t/3+scale*t0*Waveform_ToolBox.M0iter(G,G2,0,t,mode));
% 
%             elseif(mode==2) %'FlapTop'
% 
%                 M1=(G*t*t/2+scale*t0*Waveform_ToolBox.M0iter(G,G2,0,t,mode));
% 
%             elseif(mode==3) %'RampDown'
% 
%                 M1=(G*t*t/6+scale*t0*Waveform_ToolBox.M0iter(G,G2,0,t,mode));
% 
%             elseif(mode==4) %'BridGe'
% 
%                 M1=(t*t*(G/6+G2/3)+scale*t0*Waveform_ToolBox.M0iter(G,G2,0,t,mode));
% 
%             else
%                 M1=0.0;
%             end
%         end
% 
%         function [M2]=M2iter(G, G2,t0, t1,mode,scale)
% 
%             t=t1-t0;
%             if(mode==1) %'RampUp'
% 
%                M2=(G2*t*t*t/4+scale*2*t0*Waveform_ToolBox.M1iter(G,G2,0,t,mode,0)+scale*t0*t0*Waveform_ToolBox.M0iter(G,G2,0,t,mode));
% 
%             elseif(mode==2) %'FlapTop'
% 
%                 M2=(G*t*t*t/3+scale*2*t0*Waveform_ToolBox.M1iter(G,G2,0,t,mode,0)+scale*t0*t0*Waveform_ToolBox.M0iter(G,G2,0,t,mode));
% 
%             elseif(mode==3) %'RampDown'
% 
%                 M2=(G*t*t*t/12+scale*2*t0*Waveform_ToolBox.M1iter(G,G2,0,t,mode,0)+scale*t0*t0*Waveform_ToolBox.M0iter(G,G2,0,t,mode));
% 
%             elseif(mode==4) %'BridGe'
% 
%                 M2=(t*t*t*(G/12+G2/4)+scale*2*t0*Waveform_ToolBox.M1iter(G,G2,0,t,mode,1/2)+scale*t0*t0*Waveform_ToolBox.M0iter(G,G2,0,t,mode));
% 
%             else
%                 M2=0;
%             end
%         end
        function [M0]=M0iter(G1, G2,t0,t1)
              t=t1-t0;
              M0=Waveform_ToolBox.Mniter_0(G1, G2,t,0);
        end
        function [M1]=M1iter(G1, G2,t0,t1)
              t=t1-t0;
              M1=Waveform_ToolBox.Mniter_0(G1, G2,t,0)*t0+Waveform_ToolBox.Mniter_0(G1, G2,t,1);
        end
        function [M2]=M2iter(G1, G2,t0,t1)
              t=t1-t0;
              M2=Waveform_ToolBox.Mniter_0(G1, G2,t,0)*t0^2+2*Waveform_ToolBox.Mniter_0(G1, G2,t,1)*t0+Waveform_ToolBox.Mniter_0(G1, G2,t,2);
        end
        function [M3]=M3iter(G1, G2,t0,t1)
              t=t1-t0;
              M3=Waveform_ToolBox.Mniter_0(G1, G2,t,0)*t0^3+3*Waveform_ToolBox.Mniter_0(G1, G2,t,1)*t0^2+3*Waveform_ToolBox.Mniter_0(G1, G2,t,2)*t0+Waveform_ToolBox.Mniter_0(G1, G2,t,3);
        end
         function [M4]=M4iter(G1, G2,t0,t1)
              t=t1-t0;
              M4=Waveform_ToolBox.Mniter_0(G1, G2,t,0)*t0^4+4*Waveform_ToolBox.Mniter_0(G1, G2,t,1)*t0^3+6*Waveform_ToolBox.Mniter_0(G1, G2,t,2)*t0^2+4*Waveform_ToolBox.Mniter_0(G1, G2,t,3)*t0^3+Waveform_ToolBox.Mniter_0(G1, G2,t,4);
        end
        function [Mn]=Mniter_0(G1, G2,t,n)
            Mn=t^(n+1)*(G1/(n+1)+G2)/(n+2);
        end
        function [M0,M1,M2,M3,M4] = Momentiter(Dn,Gn,Rt)
               
               t0=0;
               t1=0;
               M0=0;
               M1=0;
               M2=0;
               M3=0;
               M4=0;
                for cpt=1:1:length(Dn)   
                       if Dn(cpt)>0 && Dn(cpt)>2*Rt && Gn(cpt)~=0   % Trapezoid    
                            t0=t1;
                            t1=t0+Rt; % RampUp
                            G0=0;
                            G1=Gn(cpt);
                            M0=M0+Waveform_ToolBox.M0iter(G0,G1,t0,t1);
                            M1=M1+Waveform_ToolBox.M1iter(G0,G1,t0,t1);
                            M2=M2+Waveform_ToolBox.M2iter(G0,G1,t0,t1);
                            M3=M3+Waveform_ToolBox.M3iter(G0,G1,t0,t1);
                            M4=M4+Waveform_ToolBox.M4iter(G0,G1,t0,t1);

                            t0=t1;
                            t1=t0+Dn(cpt)-2*Rt; % Flaptop
                            G0=Gn(cpt);
                            G1=Gn(cpt);
                           M0=M0+Waveform_ToolBox.M0iter(G0,G1,t0,t1);
                            M1=M1+Waveform_ToolBox.M1iter(G0,G1,t0,t1);
                            M2=M2+Waveform_ToolBox.M2iter(G0,G1,t0,t1);
                            M3=M3+Waveform_ToolBox.M3iter(G0,G1,t0,t1);
                            M4=M4+Waveform_ToolBox.M4iter(G0,G1,t0,t1);

                            t0=t1;
                            t1=t0+Rt;  % RampDown
                            G0=Gn(cpt);
                            G1=0;
                            M0=M0+Waveform_ToolBox.M0iter(G0,G1,t0,t1);
                            M1=M1+Waveform_ToolBox.M1iter(G0,G1,t0,t1);
                            M2=M2+Waveform_ToolBox.M2iter(G0,G1,t0,t1);
                            M3=M3+Waveform_ToolBox.M3iter(G0,G1,t0,t1);
                            M4=M4+Waveform_ToolBox.M4iter(G0,G1,t0,t1);
                           
                       elseif Dn(cpt)>0 && Dn(cpt)<=2*Rt  && Gn(cpt)~=0  % Triangle
                            
                            Ampl=Gn(cpt)*Dn(cpt)/Rt/2;
                            
                            t0=t1;
                            t1=t0+Dn(cpt)/2; % RampUp
                             G0=0;
                            G1=Gn(cpt);
                           M0=M0+Waveform_ToolBox.M0iter(G0,G1,t0,t1);
                            M1=M1+Waveform_ToolBox.M1iter(G0,G1,t0,t1);
                            M2=M2+Waveform_ToolBox.M2iter(G0,G1,t0,t1);
                            M3=M3+Waveform_ToolBox.M3iter(G0,G1,t0,t1);
                            M4=M4+Waveform_ToolBox.M4iter(G0,G1,t0,t1);

                            t0=t1;
                            t1=t0+Dn(cpt)/2; % RampDown
                             G0=Gn(cpt);
                            G1=0;
                            M0=M0+Waveform_ToolBox.M0iter(G0,G1,t0,t1);
                            M1=M1+Waveform_ToolBox.M1iter(G0,G1,t0,t1);
                            M2=M2+Waveform_ToolBox.M2iter(G0,G1,t0,t1);
                            M3=M3+Waveform_ToolBox.M3iter(G0,G1,t0,t1);
                            M4=M4+Waveform_ToolBox.M4iter(G0,G1,t0,t1);
                       else                             % Blank
                           % zero
                            t0=t1;
                            t1=t0+Dn(cpt);  % Flaptop
                            G0=0;
                            G1=0;
                            M0=M0+Waveform_ToolBox.M0iter(G0,G1,t0,t1);
                            M1=M1+Waveform_ToolBox.M1iter(G0,G1,t0,t1);
                            M2=M2+Waveform_ToolBox.M2iter(G0,G1,t0,t1);
                            M3=M3+Waveform_ToolBox.M3iter(G0,G1,t0,t1);
                            M4=M4+Waveform_ToolBox.M4iter(G0,G1,t0,t1);
                       end
                end
        end
%            function [Ed_rsp]= EddyCurrent2(G,dT,lambda)
%                 Ed_rsp=[];
%                 N_lam = 200;
%                 all_lam = linspace(1e-4,lambda,N_lam);
%                 all_e = zeros(1, N_lam);
%                 ii = 1;
%                 for lam = all_lam
%                     lam_s = lam * 1.0e-3;
%                     r = diff(exp(-[1:numel(G)+1].*dT./lam_s));
%                     r = r(end:-1:1);
%                     e = dot(r,G);
%                     all_e(ii) = e;
%                     ii = ii + 1;
%                 end
% 
%                 if all_e(2) < 0
%                     all_e = -all_e;
%                 end
%          end
        function [Ed_rsp]= EddyCurrent(G,dT,lambda)
                  t=(0:1:length(G)-1)*dT;
                  dG=diff(G);
                  dG=[0; dG];
                  
%                   calculate time course of eddy currents for chosen lambda
%                     eddy_response = exp(-(1:length(dGdt))/(LAMBDA0*(1e-3)/dt));
%                     eddy_plot = conv(dGdt,eddy_response);
%                     eddy_plot = eddy_plot(1:n);
                  
                  dExpL=exp(-t./lambda)';
                  
                  Ed_rsp=conv(dG,dExpL);
                  Ed_rsp=Ed_rsp(1:length(G));
         end
          function [Ed_rsp]= EddySprectrum(G,dT,lam_vect) 
              Ed_rsp=[];
                for cpt=1:1:length(lam_vect)
                    Ed_rsp(:,cpt)=Waveform_ToolBox.EddyCurrent(G,dT,lam_vect(cpt)) ;
                end
                
         end
          function [Stm]= PNS(G,dT)
                    alpha = 0.333;
                    r = 23.4;
                    c = 334e-6;
                    Smin = 60;
                    coeff = zeros(numel(G), 1);
                    for cpt = 1:(numel(G))
                        coeff(cpt) = ( c / ((c + dT*(numel(G)-1) - dT*(cpt-1))^2.0) / Smin );
                    end

                    Stm = zeros(numel(G), 1);
                    for cpt2 = 1:(numel(G)-2)
                        ss = 0;
                        for cpt = 1:(cpt2+1)
                            ss = ss + coeff( numel(coeff)-cpt2+cpt-1 ) * (G(cpt+1)-G(cpt));
                        end
                        Stm(cpt2) = ss;
                    end

          end
          
          function [G,t]= Balanced_Factory(Dn,Gn,Rt,dT,Bval)
                
              
                %%% Build a Structure with all the information we need
                [GM] = Waveform_ToolBox.BuildGStruct(Dn,Gn,Rt)
                [GM, M0b, M0a] =Waveform_ToolBox.Balancer(GM);

                RTfunction=@(t)t/Rt;
                
                %%% Generate a vector of gradient
                [G_tmp] = Waveform_ToolBox.GStruct2GVect(GM,RTfunction);
                G=G_tmp(1:dT:end); % T/m
                t=1:dT: GM(end).StopTime;
                
                if ~(nargin<5||isempty(Bval))
                    b_fact=Waveform_ToolBox.Bcalc(G,dT);                
                    Ampl=sqrt(Bval./b_fact); % mT/m
                    G=repmat(G,1,length(Ampl)).*Ampl;
                end

          end
          function [Dn,Gn]= Revert_Factory(G,Rt,dT)
          % This function is a copy paste (converted to matlab) to the code
          % implemented in the sequence to estimate the gradient from the waveform vector.
          % If it buggy here the sequence is too. 
          Dn=[];
          Gn=[];
          
              cpt_gr=1;
              Save_Grad=[];
              Grad_tmp=[];
              Grad_tmp_n=[];
              Grad_tmp_z=[];
               
              trig=false;
              trig_n=false;
              trig_z=false;
              
              dMaxG=0;
              MaxAmpl=max(abs(G));
              
              for cpt=1:1:length(G)

                      dNegG=0;
                      dPosG=0;
                      if G(cpt)>0
                          dPosG=1;
                      elseif G(cpt)<0
                          dNegG=1;
                      end

                       if abs(G(cpt))>dMaxG
                          dMaxG=abs(G(cpt));
                      end



    %          ////////// Positive Gradient ////////////
                    if (~trig)&&(dPosG~=0)   %%% Rampup from 0
                        Grad_tmp.dStart=(cpt-1)*dT;
                        trig=true;
                    end

                    if (trig&&(dPosG==0||cpt==length(G)))  %%% End of pulse

                        Grad_tmp.dStop=(cpt-1)*dT;
                        Grad_tmp.lSign=1;
                        Grad_tmp.dAmpl=1;
                        Grad_tmp.dDuration=Grad_tmp.dStop-Grad_tmp.dStart;
                       % Grad_tmp.lDuration_us=fSDSRoundDownGRT(static_cast<long>(Grad_tmp.dDuration*1e6));
                        Grad_tmp.lRampTime=Rt;
                        if(Grad_tmp.dDuration<(2*Rt)) % Triangle case
                            Grad_tmp.dAmpl=dMaxG/MaxAmpl;
                            %Grad_tmp.lRampTime=fSDSRoundDownGRT(Grad_tmp.lDuration_us/2);
                        end


                        Save_Grad(cpt_gr).grad=Grad_tmp;
                        cpt_gr=cpt_gr+1;

                        Grad_tmp=[];
                        trig=false;
                        dMaxG=0;

                    end

    %         ////////// Negative Gradient ////////////
                if ((~trig_n)&&(dNegG~=0))   %%% Rampup from 0       
                    Grad_tmp_n.dStart=(cpt-1)*dT;
                    trig_n=true;
                end

                if (trig_n&&(dNegG==0||cpt==length(G)))  %%% End of pulse

                    Grad_tmp_n.dStop=(cpt-1)*dT;
                    Grad_tmp_n.lSign=-1;
                    Grad_tmp_n.dAmpl=-1;
                    Grad_tmp_n.dDuration=Grad_tmp_n.dStop-Grad_tmp_n.dStart;
                    %Grad_tmp_n.lDuration_us=fSDSRoundDownGRT(static_cast<long>(Grad_tmp_n.dDuration*1e6));
                    Grad_tmp_n.lRampTime=Rt;
                    if(Grad_tmp_n.dDuration<(2*Rt))

                        Grad_tmp_n.dAmpl=-dMaxG/MaxAmpl;
                        %Grad_tmp_n.lRampTime=fSDSRoundDownGRT(Grad_tmp_n.lDuration_us/2);
                    end

                    Save_Grad(cpt_gr).grad=Grad_tmp_n;
                    cpt_gr=cpt_gr+1;

                    Grad_tmp_n=[];
                    trig_n=false;
                    dMaxG=0;
                end
            
                
             %         ////////// Fill Time ////////////
                 if (~trig_z)&&(dNegG==0)&&(dPosG==0)  %%% Trig 0
                    Grad_tmp_z.dStart=(cpt-1)*dT;
                    trig_z=true;
                 end   
                 
                  if trig_z&&(trig_n||trig||cpt==length(G))  %%% End of fill time
                    Grad_tmp_z.dStop=(cpt-1)*dT;
                    Grad_tmp_z.lSign=0;
                    Grad_tmp_z.dAmpl=0;
                    Grad_tmp_z.dDuration=Grad_tmp_z.dStop-Grad_tmp_z.dStart;
                    
                    Save_Grad(cpt_gr).grad=Grad_tmp_z;
                    cpt_gr=cpt_gr+1;

                    Grad_tmp_z=[];
                    trig_z=false;
                   
                  end
    
              end
            
              %%%% Save into the struct format
              for cpt_g=1:1:length(Save_Grad)
                        
                  Dn=[Dn; Save_Grad(cpt_g).grad.dDuration];
                  Gn=[Gn; Save_Grad(cpt_g).grad.dAmpl];
              end             
          end
          
          
          function [G,t]= Factory(Dn,Gn,Rt,dT,Bval,RTfunction)
                
              
                % This function build a waveform base on the timing and in mT/m 
                % dT in us
                % This give a b-value in  in s/mm^2
                %
                %
                %
                % Example :
                % 
                % Dn=[7000 273000 7000];        % us
                % Gn=[1  0 -1];                 % Normalized or mT/m
                % RampTime=1520;                % us
                % dT=10;                        % us
                % Bval=[200 350];               % s/mm^2
                % RTfunction=@(t)(1 -1*exp(-0.002*t-0.000001*t.*t-0.000000001*t.*t.*t));
                % 
                % [G,t]= Waveform_ToolBox.Factory(Dn,Gn,RT,dT,Bval,RTfunction);
                %
                
                %%% The user can define it's own Ramp function, otherwise linear %%%
                if nargin<6 || isempty(RTfunction)
                    RTfunction=@(t)t/Rt;
                end                
                
                %%% Build a Structure with all the information we need
                [GM] = Waveform_ToolBox.BuildGStruct(Dn,Gn,Rt);

                %%% Generate a vector of gradient
                [G_tmp] = Waveform_ToolBox.GStruct2GVect(GM,RTfunction);
                G=G_tmp(1:dT:end); % mT/m
                t=1:dT: GM(end).StopTime;
                
                if ~(nargin<5||isempty(Bval))
                    b_fact=Waveform_ToolBox.Bcalc(G,dT);                
                    Ampl=sqrt(Bval./b_fact); % mT/m
                    G=repmat(G,1,length(Ampl)).*Ampl;
                end

          end
          function [G,t]= Factory2(Dn,Gn,Rst,Rt,dT,Bval,RTfunction)
                
              
                % This function build a waveform base on the timing and in mT/m 
                % dT in us
                % This give a b-value in  in s/mm^2
                %
                %
                %
                % Example :
                % 
                % Dn=[7000 273000 7000];        % us
                % Gn=[1  0 -1];                 % Normalized or mT/m
                % RampTime=1520;                % us
                % dT=10;                        % us
                % Bval=[200 350];               % s/mm^2
                % RTfunction=@(t)(1 -1*exp(-0.002*t-0.000001*t.*t-0.000000001*t.*t.*t));
                % 
                % [G,t]= Waveform_ToolBox.Factory(Dn,Gn,RT,dT,Bval,RTfunction);
                %
                
                %%% The user can define it's own Ramp function, otherwise linear %%%
                if nargin<7 || isempty(RTfunction)
                    RTfunction=@(t)t/Rt;
                end                
                
                %%% Build a Structure with all the information we need
                [GM] = Waveform_ToolBox.BuildGStruct(Dn,Gn,Rt);

                %%% Generate a vector of gradient
                [G_tmp] = Waveform_ToolBox.GStruct2GVect(GM,RTfunction);
                G1=G_tmp(1:dT:end); % mT/m
               
                %%% Regenare only the relevant gradients
                Dn2=Dn(find(Rst));
                Gn2=Gn(find(Rst));
                [GM] = Waveform_ToolBox.BuildGStruct(Dn2,Gn2,Rt);
                %%% Generate a vector of gradient
                [G_tmp] = Waveform_ToolBox.GStruct2GVect(GM,RTfunction);
                G=G_tmp(1:dT:end); % mT/m
                t=1:dT: GM(end).StopTime;

                if ~(nargin<6||isempty(Bval))
                    b_fact=Waveform_ToolBox.Bcalc(G1,dT);                
                    Ampl=sqrt(Bval./b_fact); % mT/m
                    G=repmat(G,1,length(Ampl)).*Ampl;
                end

          end
          function [GM] = BuildGStruct(Dn,Gn,Rt)
               tmp_time=0;
               GM=[];
                for cpt=1:1:length(Dn)   
                       if Dn(cpt)>0 && Dn(cpt)>2*Rt    % Trapezoid
                        GM(cpt).StartTime=tmp_time;
                        GM(cpt).RampTime=Rt;
                        GM(cpt).FlapTop=Dn(cpt)-2*GM(cpt).RampTime;
                        GM(cpt).StopTime=GM(cpt).StartTime+GM(cpt).FlapTop+2*GM(cpt).RampTime;
                        GM(cpt).EffectiveTime=GM(cpt).FlapTop+GM(cpt).RampTime;
                        GM(cpt).Ampl=Gn(cpt);
                        tmp_time=GM(cpt).StopTime;
                       elseif Dn(cpt)>0 && Dn(cpt)<=2*Rt    % Triangle
                        GM(cpt).StartTime=tmp_time;
                        GM(cpt).RampTime=Dn(cpt)/2;
                        GM(cpt).FlapTop=0;
                        GM(cpt).StopTime=GM(cpt).StartTime+GM(cpt).FlapTop+2*GM(cpt).RampTime;
                        GM(cpt).EffectiveTime=GM(cpt).FlapTop+GM(cpt).RampTime;
                        GM(cpt).Ampl=Gn(cpt)*GM(cpt).RampTime/Rt;
                        tmp_time=GM(cpt).StopTime; 
                       else                             % Blank
                        GM(cpt).StartTime=tmp_time;
                        GM(cpt).RampTime=0;
                        GM(cpt).FlapTop=0;
                        GM(cpt).StopTime=GM(cpt).StartTime;
                        GM(cpt).EffectiveTime=0;
                        GM(cpt).Ampl=0;
                       end
                end
              
          end
               
          function [G] = GStruct2GVect(GM,RTfunction)
   
                G=zeros(round(GM(end).StopTime),1);
                for cpt=1:1:length(GM)
                    if GM(cpt).EffectiveTime~=0  

                        t1=GM(cpt).StartTime+1;
                        t2=GM(cpt).StartTime+GM(cpt).RampTime;
                        %G(t1:1:t2)=GM(cpt).Ampl*RTfunction(1:1:GM(cpt).RampTime); % Ramp1
                        G(t1:1:t2)=GM(cpt).Ampl*RTfunction(1:1:GM(cpt).RampTime)./RTfunction(GM(cpt).RampTime); % Ramp1
                        
                        t1=t2+1;
                        t2=t2+GM(cpt).FlapTop;
                        G(t1:1:t2)=GM(cpt).Ampl;     % G

                        t1=t2+1;
                        t2=t2+GM(cpt).RampTime;
                        G(t1:1:t2)=GM(cpt).Ampl*RTfunction(GM(cpt).RampTime:-1:1)./RTfunction(GM(cpt).RampTime); % Ram
                        %G(t1:1:t2)=GM(cpt).Ampl*(RTfunction(GM(cpt).RampTime:-1:1)); % Ramp1   
                    end        
                end 
          end
            
          function [GM,M0b,M0a] = Balancer(GM)
                %Todo
             
                M0b=0;
                M0a=0;
                
                for cpt=1:1:length(GM)
                    M0b=M0b+(GM(cpt).FlapTop+GM(cpt).RampTime)*GM(cpt).Ampl;
                end
                M0b
                for cpt=1:1:length(GM)
                   if (GM(cpt).Ampl>0)
                       if M0b>0
                            GM(cpt).Ampl=GM(cpt).Ampl-M0b/(GM(cpt).FlapTop+GM(cpt).RampTime);
                        break;
                       end
                   elseif (GM(cpt).Ampl<0)
                       if M0b<0
                            GM(cpt).Ampl=GM(cpt).Ampl-M0b/(GM(cpt).FlapTop+GM(cpt).RampTime);
                        break;
                       end
                   end
                    
                end
                
                 for cpt=1:1:length(GM)
                    M0a=M0a+(GM(cpt).FlapTop+GM(cpt).RampTime)*GM(cpt).Ampl;
                 end
                M0a
          end
          
          
          function [G,Bmax, TEmin]=WaveformTEMin(Btarget, T90,T180,TEepi,Rt,Gmax,dT,Type)
              
              % Find the TE by bisection, we stop when we have the b-value
             % +/- 10
              eps=10; % bvalue 
              err=999999999;
              
              TE1=30000;
              TE2=250000;
%               TE_prev1=TE1;
%               TE_prev2=TE2;
%               itermax=1000;
%               iter=1;
%               tmpTE=[];
%               tmpb=[];
%                while (err >= eps)&(iter<=itermax)
%                       
%                      
%                       [~,  b1]=  Waveform_ToolBox.WaveformComposer(TE1,T90,T180,TEepi,Rt,Gmax,dT,Type)  ;
%                       [~,  b2]=  Waveform_ToolBox.WaveformComposer(TE2,T90,T180,TEepi,Rt,Gmax,dT,Type)  ;
%                         
%                         if abs(b1-Btarget)>abs(b2-Btarget) % change the low limit
%                             TE1=TE1+abs(TE2-TE1)/4;
%                             err=abs(b2-Btarget);
%                             TEmin=TE2;
%                         else                            % change the high limit
%                             TE2=TE2-abs(TE2-TE1)/4;
%                             err=abs(b1-Btarget);
%                             TEmin=TE1;
%                         end
%                          iter=iter+1;
%                end

               % [G,  Bmax]=  Waveform_ToolBox.WaveformComposer(TEmin,T90,T180,TEepi,Rt,Gmax,dT,Type)  ;
                tmpTE=linspace(0,TE2,100);
                for cpt=1:1:100
                    
                    [G,  Bmax]=  Waveform_ToolBox.WaveformComposer(tmpTE(cpt),T90,T180,TEepi,Rt,Gmax,dT,Type)  ;
                    tmpb(cpt)=Bmax;
                end
                idx=find(tmpb~=0);
                TEmin = interp1(tmpb(idx),tmpTE(idx),Btarget);
                TEmin=ceil(TEmin/10)*10; % round to dT=10;
                [G,  Bmax]=  Waveform_ToolBox.WaveformComposer(TEmin,T90,T180,TEepi,Rt,Gmax,dT,Type)  ;
             %  figure,plot(tmpb,tmpTE);
            end
            
                  
                
               
          function [G,Bmax]=WaveformComposer(TE,T90,T180,TEepi,Rt,Gmax,dT,Type,Bval,RTfunction)     
              lt = ceil( (TE/2 - T180/2 - max(T90, TEepi)) /10 )*10 ;       
              
              FillPre=TE/2- T90 -  T180/2 -lt;
              FillPost=TE/2- TEepi -  T180/2 -lt;
              Delta=lt+FillPre+T180+FillPost;
              
              [FillPre FillPost Delta];
              if Type==1 % AMC
                  [D1 D2]=Waveform_ToolBox.M2Solver(Delta,lt,Rt);
                  Dn=[T90 D1 D2 FillPre T180 FillPost D2 D1];
                  Gn=[0   Gmax -Gmax 0       0    0       Gmax -Gmax ];
              elseif Type==2 %Bip 
                  Dn=[T90 lt/2 lt/2 FillPre T180 FillPost lt/2 lt/2];
                  Gn=[0   Gmax -Gmax 0       0    0       -Gmax Gmax ];
              elseif Type==3 % AMC inv
                  [D1 D2]=Waveform_ToolBox.M2Solver(Delta,lt,Rt);
                  Dn=[T90 D2 D1 FillPre T180 FillPost D1 D2];
                  Gn=[0   Gmax -Gmax 0       0    0       Gmax -Gmax ];
              else %Mono
                  Dn=[T90 lt  FillPre T180 FillPost  lt];
                  Gn=[0   Gmax  0       0    0        -Gmax ];
                 
              end
              
              if nargin<9 || isempty(Bval)
                    Bval=[];
              end   
              
              if nargin<10 || isempty(RTfunction)
                    RTfunction=@(t)t/Rt;
              end     
              
              [G,t]= Waveform_ToolBox.Factory(Dn,Gn,Rt,dT,Bval,RTfunction);
               Bmax =Waveform_ToolBox.Bcalc(G,dT);
          end
          
          function [D1 D2]=M2Solver(Delta,lt,rt)
            dT=10;
            Ratio=999999999;
            Delt2=0;
            check=0;
            for Delt1=0:dT:lt/2
                Delt2=Delt1*(Delta-rt)/(Delta+rt-2*Delt1);
                if(Ratio>abs((Delt1+Delt2)-lt))
                   Ratio= abs((Delt1+Delt2)-lt);
                   D1=Delt1;
                   D2=Delt2;
                   check=1;
                end
            end
            if check
            D2=ceil(D2/dT)*dT;
            else
                D2=0;
                D1=0;
            end
          end

          function [G]=EPIComposer(n,G,D,rt)
              sign =1;
              Dn=[];
              Gn=[];
              for cpt=1:1:n
                  if mod(cpt,2)==0 % Even
                    sign=1;
                  else % Odd
                    sign=-1;
                  end          
                  Dn=[Dn D];
                  Gn=[Gn sign*G];
              end
              %%% Build a Structure with all the information we need
              [GM] = Waveform_ToolBox.BuildGStruct(Dn,Gn,rt);
              RTfunction=@(t)t/rt;
              G=Waveform_ToolBox.GStruct2GVect(GM,RTfunction);
          end
    end   
    
end