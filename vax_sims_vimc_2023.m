%% Simulations for the big burden model

% Pull in burden dynamic link, which includes country data and agedist
%load('./burden_dynamic_link_08Dec2017b.mat') %IHMEs is the same
load('burden_dynamic_link_Yale.mat');
load('model_samples_ALL_Aug2021.mat');
clear R0country R0modelcasescountry repcountry avgage_country

% Import list of VIMC countries
country_list = readtable('typhoid-93-list_GBDregions.csv','ReadVariableNames',true);

vimc_country_names=country_list.Country; %gavi73.country_name; %{'Congo, the Democratic Republic of the','Ethiopia','India','Nigeria','Pakistan'};  
vimc_countries=country_list.ISO3; %gavi73.ISO3; %{'COD','ETH','IND','NGA','PAK'};  
country_regions=country_list.GBD_region;

% Import list of 145 countries
all_countries = readtable('countrycodes_iso.txt','ReadVariableNames',true);

iso_lmic = iso;
iso_lmic(strcmp(iso_lmic.countryiso, 'TWN'),:) = [];

%% Bring in coverage data

%coverage for VIMC countries
vimc_cov = readtable('coverage_202310gavi-4_typhoid-campaign-routine-default.csv','ReadVariableNames',true);
vimc_cov_bluesky = readtable('coverage_202310gavi-4_typhoid-campaign-routine-bluesky.csv','ReadVariableNames',true);

vimc_cov(vimc_cov.year<2000, :) = [];
vimc_cov_bluesky(vimc_cov_bluesky.year<2000, :) = [];


%how many vimc countries are there?
ncountry=length(vimc_countries);

% Bring in population data
input_pop_data = readtable('202310gavi-4_dds-202208_int_pop_both.csv','ReadVariableNames',true);
input_pop_data(:,5:121)=[];
vimc_pop_mat=zeros(101,101,ncountry);

% Age-independent fixed parameters
params.delta = 1/4;
params.alpha = 0.01; % 0.01;
params.omega = 1/104; % -log(1-.25)/52;
% params.epsilon = 2.68*10^(-11);

% Age-specific fixed parameters
agepar.u = [1/(0.75*52); 1/(1.25*52); 1/(3*52); 1/520; 1/520; 0];
agepar.theta = [0.003; 0.003; 0.003; 0.003; 0.021; 0.021];
agepar.theta2 = zeros(6,1);

% Birth rates
mub = [36.6; 23.6; 15.0]; 
% ESA estimate of birth rate in low income: 36.6 for 2010-15
% ESA estimate of birth rate in middle income: 19.7 for 2010-15
% ESA estimate of birth rate in lower-middle income: 23.6 for 2010-15
% ESA estimate of birth rate in upper-middle income: 15.0 for 2010-15

% These are only these ages: 0-4, 5-9, 10-14, 15-20, 20-25. The rest are
% 1-sum(distribution)
low_inc = [0.158914104, 0.14116809, 0.124802931, 0.108397066, 0.091790805];
lmid_inc = [0.108066773, 0.103190497, 0.098970907, 0.094718397, 0.090066389];
umid_inc = [0.072224644, 0.068942809, 0.066816143, 0.068359158, 0.080071929];
% mid_inc = [0.091922178	0.087764088	0.084487406	0.08284539	0.085564619];
% Download from ESA. https://esa.un.org/unpd/wpp/DataQuery/

typesoutput = {'lowinc'; 'lmid'; 'umid'};

% Bring in data on life expectancy at birth (for DALYs calculation)
%life_expectancy = readtable('202110gavi-2_dds-201910_lx0_both.csv','ReadVariableNames',true);
life_expectancy = readtable('202310gavi-4_dds-202208_life_ex_both.csv','ReadVariableNames',true);

% age groups: 0-9m, 9m-2y, 2y-4y, 5y-14y, 15y-25y, 25y+
al = 6;

% Bring in R0, m1, m2, rC, and rep
%load('model_samples.mat')

% New vaccination parameters (same 2000 draws for all countries.)
load('vax_samp2000dec.mat')

nruns = 20; %original for Bilcke paper is 2000, set to 200 for full model runs

%empty vectors to store stochastic parameter estimates (to be used later)
h2o_nruns=zeros(nruns,1);
rep=zeros(nruns,ncountry);

%empty structure to store output
output2023 = struct([]);

% what happens in matlab at values of R0 close to 1?
warning('off', 'MATLAB:warn_r14_stucture_assignment') 

%% Calculate multipliers for deaths and Dalys

% same method as Joke uses in R code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Antimicrobial Resistance %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%proportion of cases AMR (uniform from 0 to 1), replicated ncountry times
propAMR = repmat(rand(nruns,1),1,ncountry);

%relative risk of deaths/dalys (uniform from 1 to 3), replicated ncountry times
rAMR = repmat(1 + 2*rand(nruns,1),1,ncountry);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multipliers for Deaths:%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate inverse logit estimates and then re-logit to get real values
phosp = repmat(normrnd(-3.13,1.67,nruns,1),1,ncountry); %COMMON estimate for all VIMC countries:

% deaths among hospitalized = hospcases * prob death given hospitalization
% pdeathhosp = normrnd(-3.067,0.888,nruns,ncountry); -- old estimate (Pieters data only)
pdeathhosp = repmat(normrnd(-3.73, 0.307*2, nruns, 1),1,ncountry); %new COMMON estimate inculding Marchello studies 
pdeathhosp_R1 = normrnd(-2.52, 0.29,nruns,1); %subsaharan africa region-specific estimate
pdeathhosp_R4 = normrnd(-4.54, 0.48,nruns,1); %south asia region-specific estimate
pdeathhosp_R5 = normrnd(-4.13, 0.29,nruns,1); %southeast asia region-specific estimate


%for parameters with country-specific estimates, replace appropriate row with country-specific values
for c=1:ncountry
    
    country=vimc_countries{c};
    region=country_regions(c);
    
    if strcmp(country,'BGD')
        phosp(:,c) = normrnd(-4.87, 0.69, nruns, 1);
    elseif strcmp(country,'IND')
        phosp(:,c) = normrnd(-2.78,0.37,nruns,1);
    elseif strcmp(country,'KEN')
        phosp(:,c) = normrnd(-2.77,1.09,nruns,1);
    elseif strcmp(country,'PAK')
        phosp(:,c) = normrnd(-4.21,0.58,nruns,1);
    elseif strcmp(country,'IDN')
        phosp(:,c) = normrnd(-4.56, 0.50, nruns, 1);
    elseif strcmp(country, 'VNM')
        phosp(:,c) = normrnd(0.25, 0.36, nruns, 1);
    end
    
    if region==1 %subsaharan Africa
        pdeathhosp(:,c) = pdeathhosp_R1;
    elseif region==4 %south asia
        pdeathhosp(:,c) = pdeathhosp_R4;
    elseif region==5 %southeast asia
        pdeathhosp(:,c) = pdeathhosp_R5;
    end
end

%Probability of hospitalization & Probability of death given hospitalization
phosp=exp(phosp)./(1+exp(phosp));
pdeathhosp=exp(pdeathhosp)./(1+exp(pdeathhosp));

%%%% Proportion of deaths occuring in hospital %%%%
propdeathhosp = repmat(1-.75*rand(nruns,1),1,ncountry); %vector of length nruns, replicated ncountry times

% proportion of cases seeking med care 
pmedcare = repmat(normrnd(.57,.09,nruns,1),1,ncountry);
pnomed = 1-pmedcare;
poutp = pmedcare-phosp; %nrun by ncountry matrix with some country-specific estimates and the rest with value of common estimate

% find any zeros and replace them with average 
pmedcare(poutp<0) =.57;
phosp(poutp<0) = .038;
poutp(poutp<0) = .57-.038; 

%%%%%%%%%%%%%%%%%%%%%%%%%
% Multipliers for Dalys:%
%%%%%%%%%%%%%%%%%%%%%%%%%

% years of life affected by disability
% hosp cases * duration of infection hosp *disability weight
% outpatient cases * duration of infection outpatient *disability weight
% cases not seeking med care * doi no med care * disability weight

% prob seeking healthcare = .57, sd .09
% proportion not seeking care = 1- prop seeking care
% prop outpatient = prop seeking care - phosp

% duration of illness (doi) hosp and outpatient = 16 days +-2
doi_hosp = repmat(normrnd(16/365.25,2/365.25,nruns,1),1,ncountry);
doi_outp = repmat(normrnd(16/365.25,2/365.25,nruns,1),1,ncountry);
% doi not seeking med care = avg 8 days, 1-100% range of 16 days
doi_nomed = repmat((16/365.25)*rand(nruns,1),1,ncountry);
% disability weights: 
% hosp: 0·21 +- 0.04
dwt_hosp = repmat(normrnd(.21,.04,nruns,1),1,ncountry);
% following bilcke at all, half of runs are simulations with outpatient 
% disability weight for infectious disease, acute, moderate = 0·053 +- 0·012 
% and half for severe = 0·21 +- 0.04
dwt_outp = repmat([normrnd(.053,.012,nruns/2,1);normrnd(.21,.04,nruns/2,1)],1,ncountry);
% for those not seeking medical care, half assumed disability weight for
% infectious disease, acute, mild = 0·005 +- 0·002
% and half for moderate 0·053 +- 0·012 
dwt_nomed = repmat([normrnd(.005,.002,nruns/2,1);normrnd(.053,.012,nruns/2,1)],1,ncountry);

%To get deaths, multiply cases by phosp*pdeathhosp/propdeathosp:
deaths_mult = phosp.*pdeathhosp./propdeathhosp;

%To get ylds, multiply cases by phosp*doihosp*dwthosp + poutp*doioutp*dwtoutp + pnomed*doinomed*dwtnomed
yld_mult = phosp.*doi_hosp.*dwt_hosp+poutp.*doi_outp.*dwt_outp+pnomed.*doi_nomed.*dwt_nomed;


for c=1:ncountry %Loop through countries 

country=vimc_countries{c};
%country_name=vimc_country_names{c};

% set z to be the index of whichever country is being modeled
z = find(strcmp(all_countries.ISO3,country));

% age distribution for this particular country
if agedist_country(z)==1
    dist = low_inc; 
elseif agedist_country(z)==2
    dist = lmid_inc;
else
    dist = umid_inc; 
end

% Age groups in ODE: 0-9m, 9m-2y, 2y-4y, 5y-14y, 15y-25y, 25y+
% Ages for population data: 0-4, 5-9, 10-14, 15-20, 20-25. 
agedist = [dist(1)/5; dist(1)/5; dist(1)*3/5; sum(dist(2:3)); sum(dist(4:5)); 1-sum(dist)];
population = 1e5*agedist;

agepar.mub = [-log(1-mub(agedist_country(z))/1000)/52; zeros(5, 1)]; 
agepar.mu = ([agepar.mub(1); agepar.u(1:(end-1), 1)].*[sum(population); population(1:(end-1),1)])./population(:,1); 
agepar.mu = agepar.mu - agepar.u;

% get population in matrix form for multiplying
vimc_pop_mat(:,:,c) = table2array(input_pop_data(strcmp(input_pop_data.country_code,country),5:105));
country_pop=vimc_pop_mat(:,:,c);

% years of life lost = deaths*life expectancy at age of death
% set up vector of years of life lost for each age

le_mat=zeros(20,101);
le_country=life_expectancy(strcmp(life_expectancy.country_code,country),:); %extract values for country
le_country(1:220,:)=[]; %delete first 220 rows, which correspond to 1950-1995
for i=1:height(le_country)
    for y=1:20
        if le_country.year(i)==1995+5*y
            for a=1:101
                if (a-1)>=le_country.age_from(i) && (a-1)<=le_country.age_to(i)
                    le_mat(y,a)=le_country.value(i);
                end
            end
        end
    end
end
le_mat=[le_mat; le_mat(end,:)];

yll_mat=zeros(101,101);
for a=1:101
    yll_mat(:,a)=interp1(1:5:101,le_mat(:,a),1:101,'spline');
end

%% Generate stochastic simulations

tic % start the timer to keep track of how long it takes to run 

for i=1:nruns 
    
estpar.mult1 = double(m1_sample(i,1)); 
estpar.mult2 = double(m2_sample(i,1)); 
estpar.r = 1; 
estpar.rC = double(rC_sample(i,1)); 
coeff = [estpar.mult1,estpar.mult1,estpar.mult2,1,1,1]; 
estpar.R0 = double(R0samples(i,z)); 
estpar.beta = estpar.R0*coeff'.*(agepar.mub(1) + params.delta)./((1+params.delta*(agedist'*agepar.theta)*estpar.rC./agepar.mub(1))*(coeff*agedist));
params.epsilon = 0;
h2o=unifrnd(-0.025,0);  

vacpar.omega_v = -log(1-vax_samp.omega(i))/52; % convert to WEEKS
vacpar.veff = vax_samp.nu(i);

    % Simulate equilibrium
    % Run to equilibrium only once and then evaluate both unc and con and all interventions with this.
    options = odeset('RelTol',1e-6);
    
    % How long to reach equilibrium
    if estpar.beta(4)<0.2
    tspan=7000*52+1;
    elseif estpar.beta(4)<0.3
    tspan=3000*52+1;
    else
    tspan=200*52+1;
    end
    
    pop0 = zeros(al*18,1);
    pop0(1:al) = population; % putting everyone in susceptible 2???
    % throw an infected in one of the four categories to start the simulation
    ageptzero = 3; 
    pop0(ageptzero) = pop0(ageptzero)-100;
    pop0(al+ageptzero) = 100;
    
    rvnone = zeros(tspan, 6);
    massvacnone = zeros(tspan, 6);

    estpar.beta_h2o = repmat(estpar.beta, 1, 52*101+1).*repmat(exp(h2o*(-15.5:(1/52):85.5)), al, 1);
    
    %Burn-in period to reach equilibrium
    [t, pop] = ode45(@(t, pop) fn_SIR_vaccine_vnv_strata(t, pop, estpar.beta_h2o(:,1), agepar.mub, agepar.u, ...
        agepar.mu, params.omega, vacpar.omega_v, estpar.r, estpar.rC, params.alpha, params.delta, ...
        agepar.theta, agepar.theta2, params.epsilon, rvnone, vacpar.veff, massvacnone, al, population), ...
        1:tspan, pop0, options);
    pop0eqm = pop(end, :);
    %save initial population
    pre_pop = pop;
    totpop = sum(pop,2); %get total population
    %get age-specific incidence
    pre_inctime_unvacc = diff(pop(:,(16*al+1):(17*al)));
    
    % For VIMC, not doing different coverage scenarios, just given coverage
    cov_est = vimc_cov;
    rcov = cov_est.coverage(strcmp(cov_est.country_code, iso_lmic.countryiso(z)) & strcmp(cov_est.activity_type, 'routine'));
    ccov = cov_est.coverage(strcmp(cov_est.country_code, iso_lmic.countryiso(z)) & strcmp(cov_est.activity_type, 'campaign'));

    % routine coverage - for 2000 to 2100
    if length(rcov)<101
        add_rows = zeros(101-length(rcov),1);
        routine_coverage = [add_rows; rcov];
    end
    routine_doses=routine_coverage'.*country_pop(1,:);
    
    tmp = repmat(routine_coverage, 1, 52);
    tmp(isnan(tmp)) = 0;
    tmp = tmp';
    tmp = tmp(:);

    if length(tmp)>1
        tmp1 = [0; tmp];
    else
        tmp1 = zeros(5253,1);
    end
    
    v1 = [zeros(length(tmp1),1), tmp1, zeros(length(tmp1),4)];
    % only the second age group gets vaccination, hence 0's for every other
    % age group
    
    % campaign coverage - for 2000 to 2100
    pre_int_length = find(routine_coverage,1)-1; % length of time before vaccine roll-out
    int_length = length(ccov); %number of years for campaign
    
    campaign_coverage=zeros(101,1);
    campaign_coverage(pre_int_length+1:pre_int_length+int_length,1)=ccov;
    
    campaign_doses=campaign_coverage'.*sum(country_pop(2:15,:));
    
    tmp2 = ccov;
    tmp2(isnan(tmp2)) = [];
    
    % this is coded to tolerate the scenario when campaigns are done over
    % the span of several years
    camp_ages = zeros(length(tmp1), 1);
    for k=1:length(tmp2)
        camp_ages((52*(pre_int_length+k-1)+2):(52*(pre_int_length+k-1)+5)) = min((1-(1-tmp2(k)).^.25),.62);
    end
    
    %repeat for IA2030
    cov_est_bluesky = vimc_cov_bluesky;
    
    %routine vaccination
    tmp_old_ia = repmat(cov_est_bluesky.coverage(strcmp(cov_est_bluesky.country_code, iso_lmic.countryiso(z)) & strcmp(cov_est_bluesky.activity_type, 'routine')), 1, 52);
    tmp_old_ia(isnan(tmp_old_ia(:,1)), :) = [];
    if length(tmp_old_ia)<101
        add_rows = zeros(101-length(tmp_old_ia),52);
        tmp_ia = [add_rows;tmp_old_ia];
    else
        tmp_ia = tm_old_ia;
    end
    
    pre_int_length_ia = find(tmp_ia,1); % length of time before vaccine roll-out
    int_length_ia = 101-pre_int_length_ia; %number of years post-intervention
    
    tmp_ia = tmp_ia';
    tmp_ia = tmp_ia(:);

    if length(tmp_ia)>1
        tmp1_ia = [0; tmp_ia];
    else
        tmp1_ia = zeros(5253,1);
    end
    
    v1_ia = [zeros(length(tmp1_ia),1), tmp1_ia, zeros(length(tmp1_ia),4)];
    
    %campaign
    tmp2_ia = cov_est_bluesky.coverage(strcmp(cov_est_bluesky.country_code, iso_lmic.countryiso(z)) & strcmp(cov_est_bluesky.activity_type, 'campaign'));    
    
    tmp2_ia(isnan(tmp2_ia)) = [];
    camp_ages_ia = zeros(length(tmp1_ia), 1);
    
    for k=1:length(tmp2_ia)
        camp_ages_ia((52*(pre_int_length_ia+k-1)+2):(52*(pre_int_length_ia+k-1)+5)) = min((1-(1-tmp2_ia(k)).^.25),.62);
    end

    % age groups: 0-9m, 9m-2y, 2y-4y, 5y-14y, 15y-25y, 25y+
    massvac{1,1} = zeros(length(camp_ages), 6); % No vaccination
    massvac{2,1} = [zeros(length(camp_ages), 1), repmat(camp_ages, 1, 3), zeros(length(camp_ages), 2)]; % Campaign only (default)
    massvac{3,1} = [zeros(length(camp_ages_ia), 1), repmat(camp_ages_ia, 1, 3), zeros(length(camp_ages_ia), 2)]; % Campaign only (bluesky)
    massvac{4,1} = [zeros(length(camp_ages), 1), repmat(camp_ages, 1, 3), zeros(length(camp_ages), 2)]; % Routine + campaign (default)
    massvac{5,1} = [zeros(length(camp_ages_ia), 1), repmat(camp_ages_ia, 1, 3), zeros(length(camp_ages_ia), 2)]; % Routine + campaign (bluesky)
    
    for int=1:5 % cycle through interventions
        if int<4
            rv = zeros(size(v1,1), size(v1,2));
        elseif int==4
            rv=v1;
        else 
            rv=v1_ia;
        end
        
        if int==1 || (int>1 && length(tmp)>1)  % Don't need to simulate the vaccination scenarios for countries with coverage=0 
        
    [t, pop] = ode45(@(t, pop) fn_SIR_vaccine_vnv_strata_h2o(t, pop, estpar.beta_h2o, agepar.mub, agepar.u, ...
        agepar.mu, params.omega, vacpar.omega_v, estpar.r, estpar.rC, params.alpha, params.delta, ...
        agepar.theta, agepar.theta2, params.epsilon, rv, vacpar.veff, massvac{int, 1}, al, population), ...
        1:(length(tmp1)), pop0eqm, options);
    
        end
    
    % chronic = pop(5:end, (5*al+1):(6*al)) + pop(:, (13*al+1):(14*al));    
    inctime_unvacc = diff(pop(:,(16*al+1):(17*al))); %Incidence among unvaccinated
    inctime_vacc = diff(pop(:,(17*al+1):(18*al))); %Incidence among vaccinated
    %get total incidence by age, weekly
    inctime_tot = inctime_unvacc+inctime_vacc;
    %get total incidence by age, year
    inctime_age = squeeze(single(sum(reshape(inctime_tot, 52, 101, al),1))); %pre_int_length+int_length should always be 101??
    %get population size
    tmp_pop=reshape(pop(:,1:al*14),length(pop),al,14); %population by age and disease state
    age_pop = sum(tmp_pop,3); %sum to get total population by age
    age_pop = age_pop(26:52:end,:); %get midpoint population of each year for division
    
    %get age-specific incidence rate per individual
    inc_rate_age = inctime_age./age_pop;
    %get age-specific incidence rate per 100K
    inc_rate_age_100K=inc_rate_age*1e5;
    
    %calculate individual year age-specific incidence
    % distribute out incidence levels to midpoint of age group
    %try endpoint
    %Age groups in ODE: 0-9m, 9m-2y, 2y-4y, 5y-14y, 15y-25y, 25y+
    x_age = [0 .375 1.375 3.5 10 20 62.5 100];
    
    % interpolate individual year age-specific incidence rate
    interp_inc=pchip(x_age,[inc_rate_age(:,1) inc_rate_age inc_rate_age(:,6)],0:100);
    

    rvac = diff(sum(pop(:,(14*al+1):(15*al)),2));
    cvac = diff(sum(pop(:,(15*al+1):(16*al)),2));
    
    % Store: 
    % output2{i,cov}{int,1} = single(1);
    output2023{i,c}{int,1}.interp_inc=interp_inc; %interpolated individual year age-specific incidence rate
    output2023{i,c}{int,1}.inc_rate_age=inc_rate_age; %age-group age-specific incidence rate
    output2023{i,c}{int,1}.inc_rate_age_100K=inc_rate_age_100K; %age-group age-specific incidence rate
    output2023{i,c}{int,1}.inc_rate_100K=sum(inctime_age,2)./sum(age_pop,2)*1e5;

    output2023{i,c}{int,1}.inctime_age=inctime_age; %age-specific incidence by year, vacc+unvacc
    output2023{i,c}{int,1}.age_pop=age_pop;  %age-specific population size over time
%     for age=1:6
%         %get age-specific incidence
%         output2{i,cov}{int,1}.inctime_unvacc(:,age) = single(sum(reshape(inctime_unvacc(:,age), 52, 10), 1));
%         output2{i,cov}{int,1}.inctime_vacc(:,age) = single(sum(reshape(inctime_vacc(:,age), 52, 10),1));
%     
%     end
   
    output2023{i,c}{int,1}.rvac = single(sum(reshape(rvac, 52, 101), 1));
    output2023{i,c}{int,1}.cvac = single(sum(reshape(cvac, 52, 101), 1));
    output2023{i,c}{int,1}.mult1 = single(estpar.mult1); 
    output2023{i,c}{int,1}.mult2 = single(estpar.mult2);
    output2023{i,c}{int,1}.rC = single(estpar.rC);
    output2023{i,c}{int,1}.beta = single(estpar.beta);
    output2023{i,c}{int,1}.R0 = single(estpar.R0);
    output2023{i,c}{int,1}.omega_v = single(vacpar.omega_v);
    output2023{i,c}{int,1}.veff = single(vacpar.veff);
    output2023{i,c}{int,1}.h2o = h2o;
   
    %% Calculate cases, deaths and dalys 
    output2023{i,c}{int,1}.cases = round(repsamples(i,z)*output2023{i,c}{int,1}.interp_inc.*vimc_pop_mat(:,:,c)');
    output2023{i,c}{int,1}.deaths = round((1-propAMR(i)+rAMR(i)*propAMR(i))*deaths_mult(i)*output2023{i,c}{int,1}.cases);
    output2023{i,c}{int,1}.dalys = (1-propAMR(i)+rAMR(i)*propAMR(i))*(yld_mult(i)*output2023{i,c}{int,1}.cases+...
        deaths_mult(i)*output2023{i,c}{int,1}.cases.*yll_mat);
    end
    
    %create vectors of stochastic parameters to be saved in spreadsheet
    estpar.mult1_nruns(i,1) = estpar.mult1;
    estpar.mult2_nruns(i,1) = estpar.mult2; 
    estpar.rC_nruns(i,1) = estpar.rC;
    vacpar.veff_nruns(i,1) = vacpar.veff;
    vacpar.omega_v_nruns(i,1) = vacpar.omega_v;
    estpar.R0_nruns(i,c) = estpar.R0;
    h2o_nruns(i,1) = h2o;
    rep(i,c) = repsamples(i,z);
end

end

toc

save('output_2023full_bothscen.mat','-v7.3'); 




