clear all

load('mature_Dem9_Long.mat')
load('State_Sub9_Prop.mat')

DensityPathogen9_Function(mature_Dem9_Long,State_Sub9)

function DensityPathogen9_Function(mature_Dem9_Long,State_Sub9)

steps = 40;
%Set the number of beds
nBeds = 9;
% Set number of age classes
nAge = 3;
%Set initial population states
Istate = [zeros(2,nBeds);  .1*State_Sub9(3,1) 0 0 0 0 0 0 0 0];
Sstate = State_Sub9 - Istate;

% Set projection matrices
pS_surv = 0.971; %survival prob for sus pop
reduction = .05;
pI_surv = (1-reduction)*pS_surv; %survival prob for inf pop
fecS = 10^6;
fecI = 10^6;
qS = .19*(1/fecS);
qI = (1-reduction)*qS;
MS_r = [pS_surv*(11/12) 0;  pS_surv*1/12 pS_surv*1]; %aging during reproduction for sus pop
MI_r = [pI_surv*(11/12) 0;  pI_surv*1/12 pI_surv*1]; %aging during reproduction for inf pop
MS_a = [0 0 0; qS pS_surv*(11/12) 0; 0 pS_surv*1/12 pS_surv*1]; %aging matrix for sus pop
MI_a = [0 0 0; qI pI_surv*(11/12) 0; 0 pI_surv*1/12 pI_surv*1]; %aging matrix for inf pop


%Bed distribution proportions
b_0 = 0.3;
b_m1 = 0.2;
b_p1 = 0.3;
b_p2 = 0.2;


mature = Sstate(3,:)+Istate(3,:);
mature = [mature ;zeros(steps,nBeds)];

for step_counter = 1:steps
    Istate_copy = Istate;
    Sstate_copy = Sstate;
    if mod(step_counter,4) ==0 
        %Reproduction season

        for bed_counter = 1:nBeds
            numI = sum(Istate_copy(:,bed_counter));
            theta = numI ~= 0;
            Sstate(2:3,bed_counter) = (1-theta*phi(bed_counter))*MS_r*Sstate_copy(2:3,bed_counter); %age
            Istate(2:3,bed_counter) = MI_r*Istate_copy(2:3,bed_counter)+theta*phi(bed_counter)*MS_r*Sstate_copy(2:3,bed_counter);

            Sstate(1,bed_counter) = Sstate(1,bed_counter)+b_0*(1-theta*phi(bed_counter))*fecS*Sstate_copy(3,bed_counter); %same bed
            idx = max(1,bed_counter-1);
            Sstate(1,idx) = Sstate(1,idx) + (bed_counter-1>0)*(b_m1*(1-theta*phi(idx))*fecS*Sstate_copy(3,bed_counter)); %goch sent upstream
            idx = min(nBeds,bed_counter+1);
            Sstate(1,idx) = Sstate(1,idx) + (bed_counter+1<=nBeds)*(b_p1*(1-theta*phi(idx))*fecS*Sstate_copy(3,bed_counter)); %goch 1 downstream
            idx = min(nBeds,bed_counter+2);
            Sstate(1,idx) = Sstate(1,idx) + (bed_counter+2<=nBeds)*(b_p2*(1-theta*phi(idx))*fecS*Sstate_copy(3,bed_counter)); %goch 2 downstream

            Istate(1,bed_counter) = Istate(1,bed_counter)+b_0*theta*phi(bed_counter)*fecS*Sstate_copy(3,bed_counter)+b_0*fecI*Istate_copy(3,bed_counter); %same bed
            idx = max(1,bed_counter-1);
            Istate(1,idx) = Istate(1,idx) + (bed_counter-1>0)*((b_m1*theta*phi(idx)*fecS*Sstate_copy(3,bed_counter))+b_m1*fecI*Istate_copy(3,bed_counter)); %goch sent upstream
            idx = min(nBeds,bed_counter+1);
            Istate(1,idx) = Istate(1,idx) + (bed_counter+1<=nBeds)*((b_p1*theta*phi(idx)*fecS*Sstate_copy(3,bed_counter))+b_p1*fecI*Istate_copy(3,bed_counter)); %goch 1 downstream
            idx = min(nBeds,bed_counter+2);
            Istate(1,idx) = Istate(1,idx) + (bed_counter+2<=nBeds)*((b_p2*theta*phi(idx)*fecS*Sstate_copy(3,bed_counter))+b_p2*fecI*Istate_copy(3,bed_counter)); %goch 2 downstream

        end
    else
        %No reproduction
        for bed_counter = 1:nBeds
        numI = sum(Istate(:,bed_counter));
        theta = numI ~= 0;
        Sstate(:,bed_counter) = (1-theta*phi(bed_counter))*MS_a*Sstate_copy(:,bed_counter);
        Istate(:,bed_counter) = theta*phi(bed_counter)*MS_a*Sstate_copy(:,bed_counter)+MI_a*Istate_copy(:,bed_counter);
        end
    end
    mature(step_counter+1,:) = Sstate(3,:)+Istate(3,:);
end
% env_cont_level
% bed_cont_level
% state
MarkerIndex = 1:4:steps;
MarkerShapes = {'^','o','s','d','v'};
h=figure('Position',.6*get(0,'ScreenSize'));
hold on
for k = 3:nBeds-2
    hold on
    plot(0:steps,mature(:,k),'LineWidth',2,'Marker',MarkerShapes(mod(k+2,5)+1),'MarkerIndices',MarkerIndex,'MarkerSize',12);
end
ylim padded;
legend('3','4','5','6','7')
title("Density Pathogen")
set(gca,'FontSize',32)
xlabel('Time');
ylabel('Mature Pop Size')

h2=figure('Position',.6*get(0,'ScreenSize'));
hold on
for k = 3:nBeds-2
    hold on
    plot(0:steps,abs(mature_Dem9_Long(:,k)-mature(:,k))./mature_Dem9_Long(:,k),'LineWidth',2,'Marker',MarkerShapes(mod(k+2,5)+1),'MarkerSize',12);
end
ylim padded;
legend('3','4','5','6','7')
title("Density Pathogen")
set(gca,'FontSize',32)
xlabel('Time');
ylabel('Percent Change')
    
    function prob_inf = phi(bed)
    alpha = 2;
    ILam = sum(Istate_copy(2:3,bed));
    Lam = sum(Sstate_copy(2:3,bed))+ILam;
    prob_inf = 1-exp(-alpha*ILam/Lam);
    end
end