function [G_out, M_out]= diff_grad_gen(Dn,Gn,Bval,Rt,dT)

   % %%%% Initialisation %%%% 
    Gamma=       267.513e6;  % rad/ (s . Tesla )
   % Bval=Bval*1e6; % s/m²

    for cpt_b=1:1:length(Bval)
        %%% Create a gradient structure %%%
        tmp_time=0;
        for cpt=1:1:length(Dn)   
               if Dn(cpt)>0 && Dn(cpt)>2*Rt    % Trapezoid
                GM(cpt).StartTime=tmp_time;
                GM(cpt).RampTime=Rt;
                GM(cpt).FlapTop=Dn(cpt)-2*GM(cpt).RampTime;
                GM(cpt).StopTime=GM(cpt).StartTime+GM(cpt).FlapTop+2*GM(cpt).RampTime;
                GM(cpt).EffectiveTime=GM(cpt).FlapTop+GM(cpt).RampTime;
                GM(cpt).Ampl=Gn(cpt);
                tmp_time=GM(cpt).StopTime;
               elseif Dn(cpt)>0 && Dn(cpt)<2*Rt    % Triangle
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

    %     %%% Balance the Moment M0 to be perfectly equal to 0 %%%
        [GM M] =Balancer(GM);

        %%% Generate a vector of gradient
        [G] = Gradient2Vect(GM,dT);

        %%%  Calculate the B-value
        G=G.*1e-3 ;                  % 1e-3 T/m -> T/mm
        n = length(G(1:10:end));
        C=tril(ones(n));
        C2 = C'*C;
        INV = ones(n,1);   
        dT_tmp=10*1e-6;             % s

        b_fact = (Gamma^2)*(G(1:10:end)*dT_tmp)'*(C2*(G(1:10:end)*dT_tmp))*dT_tmp;
        Ampl=sqrt(Bval(cpt_b)/b_fact); % T/mm

        G_out(:,cpt_b)=G.*Ampl;  % T/mm

        b_fact = (Gamma^2)*(G(1:10:end)*dT_tmp)'*(C2*(G(1:10:end)*dT_tmp))*dT_tmp;

        M_out(:,cpt_b)=M.*Ampl; % T/mm * s^n

        Venc= pi /(Gamma*abs(M(2)));  % pi /  (rad/(s . Tesla ) & T/mm * s^2 ) => 1 / (s/mm) => mm/s
    end
end
    
function [GM M] =Balancer(GM)

    [M2 M1 M0] =MomentCalculator(GM);

    tmp_Pos=0;
    tmp_Neg=0;
    for cpt=1:1:length(GM)   
            if GM(cpt).Ampl>0
               tmp_Pos=tmp_Pos+GM(cpt).EffectiveTime*abs(GM(cpt).Ampl);
            elseif GM(cpt).Ampl<0
               tmp_Neg=tmp_Neg+GM(cpt).EffectiveTime*abs(GM(cpt).Ampl);
            end
    end

    tmp_M0=tmp_Pos-tmp_Neg;
    
    for cpt=length(GM):-1:1   
            if tmp_M0>0  % More positive than negative gradients
                if GM(cpt).Ampl>0
                    GM(cpt).Ampl=1-(tmp_M0)/(GM(cpt).EffectiveTime);
                    break;
                end
            else     % More negative than positive gradients
                if GM(cpt).Ampl<0
                    GM(cpt).Ampl=-1-(tmp_M0)/(GM(cpt).EffectiveTime);
                    break;
                end
            end         
    end
    
   [M2 M1 M0] =MomentCalculator(GM);  
    M=[M0 M1 M2];
    

end
function [G] = Gradient2Vect(GM,dT)
   
    G=zeros(round(GM(end).StopTime/dT),1);
    for cpt=1:1:length(GM)
        if GM(cpt).EffectiveTime~=0  
            
            t1=GM(cpt).StartTime+1;
            t2=GM(cpt).StartTime+GM(cpt).RampTime;
            G(t1:dT:t2)=GM(cpt).Ampl*(1:dT:GM(cpt).RampTime)./GM(cpt).RampTime; % Ramp1
            
            t1=t2+1;
            t2=t2+GM(cpt).FlapTop;
            G(t1:dT:t2)=GM(cpt).Ampl;     % G
            
            t1=t2+1;
            t2=t2+GM(cpt).RampTime;
            
            G(t1:dT:t2)=GM(cpt).Ampl*(1-(1:dT:GM(cpt).RampTime)./GM(cpt).RampTime); % Ramp1   
        end        
    end
    
end

function [M2 M1 M0] =MomentCalculator(G)

M0=0;
M1=0;
M2=0;
t0=0;
t1=0;

for cpt=1:1:size(G,2)
    t0=G(cpt).StartTime;
    t1=G(cpt).StartTime+G(cpt).RampTime;
    M0=M0+Moment0(0,G(cpt).Ampl,t0,t1,'RampUp');
    M1=M1+Moment1(0,G(cpt).Ampl,t0,t1,'RampUp');
    M2=M2+Moment2(0,G(cpt).Ampl,t0,t1,'RampUp');

    t0=t1;
    t1=t1+G(cpt).FlapTop;
    M0=M0+Moment0(G(cpt).Ampl,G(cpt).Ampl,t0,t1,'FlapTop');
    M1=M1+Moment1(G(cpt).Ampl,G(cpt).Ampl,t0,t1,'FlapTop');
    M2=M2+Moment2(G(cpt).Ampl,G(cpt).Ampl,t0,t1,'FlapTop');

    t0=t1;
    t1=t1+G(cpt).RampTime;
    M0=M0+Moment0(G(cpt).Ampl,0,t0,t1,'RampDown');
    M1=M1+Moment1(G(cpt).Ampl,0,t0,t1,'RampDown');
    M2=M2+Moment2(G(cpt).Ampl,0,t0,t1,'RampDown');
end



end

function M0=Moment0(G,G2,t0,t1,event)

    t=t1-t0;
    if(strcmp(event,'RampUp'))
        M0=(G2*t/2);
    elseif(strcmp(event,'FlapTop'))
        M0=(G*t);
    elseif(strcmp(event,'RampDown'))
        M0=G*t/2;
    elseif(strcmp(event,'Bridge'))
        M0=t*(G+G2)/2;
    end

end

function M1=Moment1(G,G2,t0,t1,event,scale)

   if(~exist('scale', 'var'))
       scale=1;
   end
    t=t1-t0;
    if(strcmp(event,'RampUp'))
        M1=(G2*t*t/3+scale*t0*Moment0(G,G2,0,t,event));
    elseif(strcmp(event,'FlapTop'))
        M1=(G*t*t/2+scale*t0*Moment0(G,G2,0,t,event));
    elseif(strcmp(event,'RampDown'))
        M1=(G*t*t/6+scale*t0*Moment0(G,G2,0,t,event));
    elseif(strcmp(event,'Bridge'))
        M1=(t*t*(G/6+G2/3)+scale*t0*Moment0(G,G2,0,t,event));
    end

end

function M2=Moment2(G,G2,t0,t1,event,scale)
   
    if(~exist('scale', 'var'))
        scale=1;
    end
    t=t1-t0;
    if(strcmp(event,'RampUp'))
        M2=(G2*t*t*t/4+scale*2*t0*Moment1(G,G2,0,t,event,1/2));
    elseif(strcmp(event,'FlapTop'))
        M2=(G*t*t*t/3+scale*2*t0*Moment1(G,G2,0,t,event,1/2));
    elseif(strcmp(event,'RampDown'))
        M2=(G*t*t*t/12+scale*2*t0*Moment1(G,G2,0,t,event,1/2));
    elseif(strcmp(event,'Bridge'))
        M2=(t*t*t*(G/12+G2/4)+scale*2*t0*Moment1(G,G2,0,t,event,1/2));
    end

end