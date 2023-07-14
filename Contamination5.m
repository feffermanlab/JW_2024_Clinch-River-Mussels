
%Set the number of beds
nBeds = 5;
% Set number of age classes
nAge = 3;
%Set initial population states
%state = [zeros(1,nBeds);100*ones(nAge-2,nBeds);100*ones(1,nBeds)];
state = State_Sub5;
% state = [zeros(1,nBeds);
%     100   30.8097   55.6399   75.9820   69.9554;
%    100  103.3045  186.5601  254.7670  234.5599]; %recovers to steady
%    state

% Set projection matrices
p_surv = 0.971;
fec = 10^6;
M_r = [p_surv*(11/12) 0;  p_surv*1/12 p_surv*1]; %aging during reproduction
M_a = [0 0 0; .21*(1/fec) p_surv*(11/12) 0; 0 p_surv*1/12 p_surv*1]; %aging matrix

absorption = 2;
omega_age = 1;
omega_repro = 0;
cont_mass = 1; %total cont mass
cont_conc = [.01 .05 .03 .02];


%Bed distribution proportions
b_0 = 0.3;
b_m1 = 0.2;
b_p1 = 0.3;
b_p2 = 0.2;

bed_pop = sum(state(2:3,:));

bed_cont_mass = zeros(1,nBeds);

steps = 40;
mature = state(3,:);
mature = [mature; zeros(steps,nBeds)];
mass_list = zeros(steps,nBeds);
env_list = zeros(steps,nBeds+1);
for step_counter = 1:steps
    %Contaminant accumulation
    env_cont_conc = [cont_conc(mod(step_counter-1,4)+1) zeros(1,nBeds)];
    WV = cont_mass/cont_conc(mod(step_counter-1,4)+1);
    old_bed_pop = bed_pop;
    bed_pop = sum(state(2:3,:)); %Mature only
    bed_mass_copy = bed_cont_mass;
    curr_mass = cont_mass;
    %update bed cont levels
    for bed_counter = 1:nBeds
        cont_removed = bed_pop(bed_counter)*absorption*env_cont_conc(bed_counter);
        lambda = old_bed_pop(bed_counter)>bed_pop(bed_counter);
        ratio = (old_bed_pop(bed_counter)-bed_pop(bed_counter))/old_bed_pop(bed_counter);
        bed_cont_mass(bed_counter) = bed_mass_copy(bed_counter)*(1-lambda*ratio)+cont_removed;
        curr_mass = max(curr_mass - cont_removed,0);
        env_cont_conc(bed_counter+1) = curr_mass/WV;
    end
    %update bed population
    if mod(step_counter,4) ==0 
        %Reproduction season
        state_copy = state;
        for bed_counter = 1:nBeds
            state(2:3,bed_counter) = exp(-omega_age*(bed_cont_mass(bed_counter)/sum(state_copy(2:3,bed_counter))))*M_r*state_copy(2:3,bed_counter); %age
            state(1,bed_counter) = state(1,bed_counter)+b_0*exp(-omega_repro*bed_cont_mass(bed_counter))*fec*state_copy(3,bed_counter); %same bed
            state(1,max(1,bed_counter-1)) = state(1,max(1,bed_counter-1)) + (bed_counter-1>0)*(b_m1*exp(-omega_repro*bed_cont_mass(bed_counter))*fec*state_copy(3,bed_counter)); %goch sent upstream
            state(1,min(nBeds,bed_counter+1)) = state(1,min(nBeds,bed_counter+1)) + (bed_counter+1<=nBeds)*(b_p1*exp(-omega_repro*bed_cont_mass(bed_counter))*fec*state_copy(3,bed_counter)); %goch 1 downstream
            state(1,min(nBeds,bed_counter+2)) = state(1,min(nBeds,bed_counter+2)) + (bed_counter+2<=nBeds)*(b_p2*exp(-omega_repro*bed_cont_mass(bed_counter))*fec*state_copy(3,bed_counter)); %goch 2 downstream

        end
    else
        %No reproduction
        state_copy = state;
        for bed_counter = 1:nBeds
        state(:,bed_counter) = exp(-omega_age*(bed_cont_mass(bed_counter)/sum(state_copy(2:3,bed_counter))))*M_a*state_copy(:,bed_counter);
        end
    end
    env_list(step_counter,:) = env_cont_conc;
    mass_list(step_counter,:) = bed_cont_mass;
    mature(step_counter+1,:) = state(3,:);
end
% env_cont_conc
% bed_cont_mass
% state
MarkerIndex = 1:4:steps;
MarkerShapes = {'^','o','s','d','v'};
h1=figure('Position',.6*get(0,'ScreenSize'));
hold on
for k = 1:nBeds
    hold on
    plot(0:steps,mature(:,k),'LineWidth',2,'Marker',MarkerShapes(k),'MarkerSize',12);
end
ylim padded;
legend('1','2','3','4','5')
title("Contaminant")
set(gca,'FontSize',32)
xlabel('Time');
ylabel('Mature Pop Size')

h2=figure('Position',.6*get(0,'ScreenSize'));
hold on
for k = 1:nBeds
    hold on
    plot(0:steps,abs(mature_Dem5_Long(:,k)-mature(:,k))./mature_Dem5_Long(:,k),'LineWidth',2,'Marker',MarkerShapes(k),'MarkerSize',12);
end
ylim padded;
legend('1','2','3','4','5')
title("Contaminant")
set(gca,'FontSize',32)
xlabel('Time');
ylabel('Percent Change')
