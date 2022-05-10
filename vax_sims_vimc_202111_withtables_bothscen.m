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

%% Bring in coverage data
% Previous data from Gavi
cov_unc = readtable('Gavi_cov_unc.csv'); %HB: just using this to get the list of gavi countries?
%cov_con = readtable('Gavi_cov_con.csv');

%coverage for VIMC countries
vimc_cov = readtable('coverage_202110gavi-3_typhoid-routine-default.csv','ReadVariableNames',true);
vimc_cov_ia2030 = readtable('coverage_202110gavi-3_typhoid-routine-ia2030_target.csv','ReadVariableNames',true);

%cov_unc(cov_unc.Year<2019, :) = [];
%cov_con(cov_con.Year<2019, :) = [];
vimc_cov(vimc_cov.year<2000, :) = [];
vimc_cov_ia2030(vimc_cov_ia2030.year<2000, :) = [];

% Make list of the Gavi73 countries only
gavi73 = table(unique(cov_unc.ISO3));
gavi73.Properties.VariableNames = {'ISO3'}; %name the column ISO3
gavi73.cn = nan(73,1); %add a column called cn and fill with NaNs


iso_lmic = iso;
iso_lmic(strcmp(iso_lmic.countryiso, 'TWN'),:) = [];
%add country names to gavi73
for i=1:73
   gavi73.cn(i) = find(strcmp(gavi73.ISO3(i), iso_lmic.countryiso));
end
for i=1:73
    x = find(strcmp(gavi73.ISO3{i},country_list.ISO3));
    gavi73.country_name(i) = country_list.Country(x);
end

%how many vimc countries are there?
ncountry=length(vimc_countries);

% Bring in population data
input_pop_data = readtable('202110gavi-2_dds-201910_2_int_pop_both.csv','ReadVariableNames',true);
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
life_expectancy = readtable('202110gavi-2_dds-201910_lx0_both.csv','ReadVariableNames',true);

% age groups: 0-9m, 9m-2y, 2y-4y, 5y-14y, 15y-25y, 25y+
al = 6;

% Bring in R0, m1, m2, rC, and rep
%load('model_samples.mat')

% New vaccination parameters (same 2000 draws for all countries.)
load('vax_samp2000dec.mat')

nruns = 200; %original for Bilcke paper is 2000, set to 200 for full model runs

%empty vectors to store stochastic parameter estimates (to be used later)
h2o_nruns=zeros(nruns,1);
rep=zeros(nruns,ncountry);

%empty structure to store output
output2 = struct([]);

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

% hospitalized cases = cases*prob hosp
% prob hosp 0·038 (0·038) [0·004-0·249] %Common estimate -- old
%logit_mean = log(.038/(1-.038)); % -- old
%logit_se = mean([abs(log(.004/(1-.004))-log(.038/(1-.038))),log(.249/(1-.249))-log(.038/(1-.038))])/1.96; %% -- old

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

%Probability of Hospitalization & Probability of death given
%hospitalization
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

% doi hosp and outpatient = 16 days +-2
doi_hosp = repmat(normrnd(16/365.25,2/365.25,nruns,1),1,ncountry);
doi_outp = repmat(normrnd(16/365.25,2/365.25,nruns,1),1,ncountry);
% doi not seeking med care = avg 8 days, 1-100% range of 16 days
doi_nomed = repmat((16/365.25)*rand(nruns,1),1,ncountry);
% disability weights: 
% hosp: 0·21 + 0·.04
dwt_hosp = repmat(normrnd(.21,.04,nruns,1),1,ncountry);
% following bilcke at all, half of runs are simulations with outpatient 
% disability weight for infectious disease, acute, moderate = 0·053 +- 0·012 
% and half for severe = 0·21 + 0·.04
dwt_outp = repmat([normrnd(.053,.012,nruns/2,1);normrnd(.21,.04,nruns/2,1)],1,ncountry);
% for those not seeking medical care, half assumed disability weight for
% infectious disease, acute, mild = 0·005 + 0·002
% and half for moderate 0·053 +- 0·012 
dwt_nomed = repmat([normrnd(.005,.002,nruns/2,1);normrnd(.053,.012,nruns/2,1)],1,ncountry);

%To get deaths, multiply cases by phosp*pdeathhosp/propdeathosp:
%deaths_mult = phosp.*pdeathhosp./propdeathhosp;
deaths_mult = phosp.*pdeathhosp./propdeathhosp;

%To get ylds, multiply cases by 
% phosp*doihosp*dwthosp + poutp*doioutp*dwtoutp
% +pnomed*doinomed*dwtnomed
%yld_mult = phosp.*doi_hosp.*dwt_hosp+poutp.*doi_outp.*dwt_outp+pnomed.*doi_nomed.*dwt_nomed;
yld_mult = phosp.*doi_hosp.*dwt_hosp+poutp.*doi_outp.*dwt_outp+pnomed.*doi_nomed.*dwt_nomed;


for c=1:ncountry %Loop through countries in test run

country=vimc_countries{c};
%country_name=vimc_country_names{c};

% set z to be the index of whichever country is being modeled
%z = str2num(getenv('OMPI_COMM_WORLD_RANK'))+1;
%z = 48; %index of nigeria in gavi73
%z = find(strcmp(gavi73.ISO3,country));
z = find(strcmp(all_countries.ISO3,country));

%j will map the gavi73 to the lmic 145
%j=gavi73.cn(z); 
j=z; %all_countries.countryn(z); 

tic

% age distribution for this particular country
if agedist_country(j)==1
    dist = low_inc; % change if we change the assumption
elseif agedist_country(j)==2
    dist = lmid_inc; % change if we change the assumption
else
    dist = umid_inc; % change if we change the assumption
end

% Age groups in ODE: 0-9m, 9m-2y, 2y-4y, 5y-14y, 15y-25y, 25y+
% Ages for population data: 0-4, 5-9, 10-14, 15-20, 20-25. 
agedist = [dist(1)/5; dist(1)/5; dist(1)*3/5; sum(dist(2:3)); sum(dist(4:5)); 1-sum(dist)];
population = 1e5*agedist;

agepar.mub = [-log(1-mub(agedist_country(j))/1000)/52; zeros(5, 1)]; 
agepar.mu = ([agepar.mub(1); agepar.u(1:(end-1), 1)].*[sum(population); population(1:(end-1),1)])./population(:,1); 
agepar.mu = agepar.mu - agepar.u;

%get population in matrix form for multiplying
vimc_pop_mat(:,:,c) = table2array(input_pop_data(strcmp(input_pop_data.country_code,country),5:105));
%test_pop_mat(:,:,c) = table2array(test5_pop(strcmp(test5_pop.country_code,country),4:104));

% years of life lost = deaths*life expectancy - age at death
%set up vector of years of life lost for each age
yr_of_birth=zeros(101,101);
for y=1:101
    for a=1:101
        yr_of_birth(y,a)=2000+y-a;
    end
end

leb_mat=zeros(101,101);
leb_country=life_expectancy(strcmp(life_expectancy.country_code,country),:);
for y=1:101
    for a=1:101
        if yr_of_birth(y,a)<=1950 
            leb_mat(y,a)=leb_country.value(1);
        elseif yr_of_birth(y,a)>2099
            leb_mat(y,a)=leb_country.value(end);
        else
            leb_mat(y,a)=leb_country.value(find(leb_country.year==yr_of_birth(y,a)));
        end
    end
end

%yll_vec = repmat([leb-[0:69],zeros(1,31)],101,1);
yll_vec = zeros(101,101);
for a=1:101
    yll_vec(:,a) = max(0,leb_mat(:,a)-a+1);
end

%% Generate stochastic simulations

tic  

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
    ageptzero = 3; %randsample(1:al, 1, true, pop0(1:al)/sum(pop0(1:al)));
    pop0(ageptzero) = pop0(ageptzero)-100;
    pop0(al+ageptzero) = 100;
    
    rvnone = zeros(tspan, 6);
    massvacnone = zeros(tspan, 6);

    %estpar.beta_h2o = repmat(estpar.beta, 1, 1+52*10).*repmat(1:-0.25/length(1:(52*10)):0.75, al, 1);
    %estpar.beta_h2o = repmat(estpar.beta, 1, 1+52*int_length).*repmat(1:-0.25/length(1:(52*int_length)):0.75, al, 1);
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
    
    % coverage - this time it is going to be from the starting point until
    % 30 years after the starting point...
    tmp_old = repmat(cov_est.coverage(strcmp(cov_est.country_code, iso_lmic.countryiso(j)) & strcmp(cov_est.activity_type, 'routine')), 1, 52);
    tmp_old(isnan(tmp_old(:,1)), :) = [];
    if length(tmp_old)<101
        add_rows = zeros(101-length(tmp_old),52);
        tmp = [add_rows;tmp_old];
    else
        tmp = tmp_old;
    end
    
    pre_int_length = find(tmp,1); % length of time before vaccine roll-out
    int_length = 101-pre_int_length; %height(rfpNGAcov); %number of years post-intervention
    
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
    
    tmp2 = cov_est.coverage(strcmp(cov_est.country_code, iso_lmic.countryiso(j)) & strcmp(cov_est.activity_type, 'campaign'));    
    
    tmp2(isnan(tmp2)) = [];
    camp_ages = zeros(length(tmp1), 1);
    
    % this is coded to tolerate the scenario when campaigns are done over
    % the span of several years
    for k=1:length(tmp2)
        camp_ages((52*(pre_int_length+k-1)+1):(52*(pre_int_length+k-1)+4)) = min((1-(1-tmp2(k)).^.25),.62);
    end
    
    %repeat for IA2030
    cov_est_ia2030 = vimc_cov_ia2030;
    
    %routine vaccination
    tmp_old_ia = repmat(cov_est_ia2030.coverage(strcmp(cov_est_ia2030.country_code, iso_lmic.countryiso(j)) & strcmp(cov_est_ia2030.activity_type, 'routine')), 1, 52);
    tmp_old_ia(isnan(tmp_old_ia(:,1)), :) = [];
    if length(tmp_old_ia)<101
        add_rows = zeros(101-length(tmp_old_ia),52);
        tmp_ia = [add_rows;tmp_old_ia];
    else
        tmp_ia = tm_old_ia;
    end
    
    pre_int_length_ia = find(tmp_ia,1); % length of time before vaccine roll-out
    int_length_ia = 101-pre_int_length_ia; %height(rfpNGAcov); %number of years post-intervention
    
    tmp_ia = tmp_ia';
    tmp_ia = tmp_ia(:);

    if length(tmp_ia)>1
        tmp1_ia = [0; tmp_ia];
    else
        tmp1_ia = zeros(5253,1);
    end
    
    v1_ia = [zeros(length(tmp1_ia),1), tmp1_ia, zeros(length(tmp1_ia),4)];
    
%campaign
    tmp2_ia = cov_est_ia2030.coverage(strcmp(cov_est_ia2030.country_code, iso_lmic.countryiso(j)) & strcmp(cov_est_ia2030.activity_type, 'campaign'));    
    
    tmp2_ia(isnan(tmp2_ia)) = [];
    camp_ages_ia = zeros(length(tmp1_ia), 1);
    
    for k=1:length(tmp2_ia)
        camp_ages_ia((52*(pre_int_length_ia+k-1)+1):(52*(pre_int_length_ia+k-1)+4)) = min((1-(1-tmp2_ia(k)).^.25),.62);
    end

    % age groups: 0-9m, 9m-2y, 2y-4y, 5y-14y, 15y-25y, 25y+
    massvac{1,1} = zeros(length(camp_ages), 6); % No vaccination
    massvac{2,1} = [zeros(length(camp_ages), 1), repmat(camp_ages, 1, 3), zeros(length(camp_ages), 2)]; % Campaign only (default)
    massvac{3,1} = [zeros(length(camp_ages_ia), 1), repmat(camp_ages_ia, 1, 3), zeros(length(camp_ages_ia), 2)]; % Campaign only (IA2030)
    massvac{4,1} = [zeros(length(camp_ages), 1), repmat(camp_ages, 1, 3), zeros(length(camp_ages), 2)]; % Routine + campaign (default)
    massvac{5,1} = [zeros(length(camp_ages_ia), 1), repmat(camp_ages_ia, 1, 3), zeros(length(camp_ages_ia), 2)]; % Routine + campaign (IA2030)
    
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
    inctime_age = squeeze(single(sum(reshape(inctime_tot, 52, pre_int_length+int_length, al),1))); %pre_int_length+int_length should always be 101??
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
    output2{i,c}{int,1}.interp_inc=interp_inc; %interpolated individual year age-specific incidence rate
    output2{i,c}{int,1}.inc_rate_age=inc_rate_age; %age-group age-specific incidence rate
    output2{i,c}{int,1}.inc_rate_age_100K=inc_rate_age_100K; %age-group age-specific incidence rate
    output2{i,c}{int,1}.inc_rate_100K=sum(inctime_age,2)./sum(age_pop,2)*1e5;

    output2{i,c}{int,1}.inctime_age=inctime_age; %age-specific incidence by year, vacc+unvacc
    output2{i,c}{int,1}.age_pop=age_pop;  %age-specific population size over time
%     for age=1:6
%         %get age-specific incidence
%         output2{i,cov}{int,1}.inctime_unvacc(:,age) = single(sum(reshape(inctime_unvacc(:,age), 52, 10), 1));
%         output2{i,cov}{int,1}.inctime_vacc(:,age) = single(sum(reshape(inctime_vacc(:,age), 52, 10),1));
%     
%     end
   
    output2{i,c}{int,1}.rvac = single(sum(reshape(rvac, 52, 101), 1));
    output2{i,c}{int,1}.cvac = single(sum(reshape(cvac, 52, 101), 1));
    output2{i,c}{int,1}.mult1 = single(estpar.mult1); 
    output2{i,c}{int,1}.mult2 = single(estpar.mult2);
    output2{i,c}{int,1}.rC = single(estpar.rC);
    output2{i,c}{int,1}.beta = single(estpar.beta);
    output2{i,c}{int,1}.R0 = single(estpar.R0);
    output2{i,c}{int,1}.omega_v = single(vacpar.omega_v);
    output2{i,c}{int,1}.veff = single(vacpar.veff);
    output2{i,c}{int,1}.h2o = h2o;
   
    %% Calculate cases, deaths and dalys 
    output2{i,c}{int,1}.cases = round(repsamples(i,z)*output2{i,c}{int,1}.interp_inc.*vimc_pop_mat(:,:,c)');
    output2{i,c}{int,1}.deaths = round((1-propAMR(i)+rAMR(i)*propAMR(i))*deaths_mult(i)*output2{i,c}{int,1}.cases);
    output2{i,c}{int,1}.dalys = (1-propAMR(i)+rAMR(i)*propAMR(i))*(yld_mult(i)*output2{i,c}{int,1}.cases+...
        deaths_mult(i)*output2{i,c}{int,1}.cases.*yll_vec);
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

%filename = ['./output2_testcountry', sprintf('%02.0f', j),'_h2o.mat'];
save('output2_202110full_bothscen_HB.mat','-v7.3'); 
%end


%% Generate output tables formatted for VIMC

% For stochastic outputs, one for baseline (no vacc), one for camp 15 (default):
%for each run id, need year in order, then age group, cohort size, cases, dalys
%stochastic_output_novacc = table; 

% For central output: same as stochastic, but mean values
central_output_novacc = table; 

% read in disease, years, age, country, cohort size, and then set up default tables as copies
% stochastic output tables
stochastic_output_novacc.disease(1:nruns*ncountry*101*101,1) = {'Typhoid'};

stochastic_output_campaign = stochastic_output_novacc;
stochastic_output_camproutine = stochastic_output_novacc;
central_output_campaign = central_output_novacc;
central_output_camproutine = central_output_novacc;

for c=1:ncountry
    
country=vimc_countries{c}; %test_countries{c};
country_name=vimc_country_names{c};   
    
stochastic_output_novacc.run_id = repmat(reshape(repmat(1:nruns,101*101,1),nruns*101*101,1),ncountry,1);
stochastic_output_novacc.year =reshape(repmat(2000:2100,1,nruns*ncountry*101),1,nruns*ncountry*101*101)';
stochastic_output_novacc.age = reshape(repmat(0:100,101,nruns*ncountry),nruns*ncountry*101*101,1);
stochastic_output_novacc.country(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = {country};
stochastic_output_novacc.country_name(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = {country_name};
stochastic_output_novacc.cohort_size(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = repmat(reshape(vimc_pop_mat(:,:,c)',101*101,1),nruns,1);

stochastic_output_campaign.run_id = repmat(reshape(repmat(1:nruns,101*101,1),nruns*101*101,1),ncountry,1);
stochastic_output_campaign.year =reshape(repmat(2000:2100,1,nruns*ncountry*101),1,nruns*ncountry*101*101)';
stochastic_output_campaign.age = reshape(repmat(0:100,101,nruns*ncountry),nruns*ncountry*101*101,1);
stochastic_output_campaign.country(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = {country};
stochastic_output_campaign.country_name(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = {country_name};
stochastic_output_campaign.cohort_size(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = repmat(reshape(vimc_pop_mat(:,:,c)',101*101,1),nruns,1);

stochastic_output_campaign_ia.run_id = repmat(reshape(repmat(1:nruns,101*101,1),nruns*101*101,1),ncountry,1);
stochastic_output_campaign_ia.year =reshape(repmat(2000:2100,1,nruns*ncountry*101),1,nruns*ncountry*101*101)';
stochastic_output_campaign_ia.age = reshape(repmat(0:100,101,nruns*ncountry),nruns*ncountry*101*101,1);
stochastic_output_campaign_ia.country(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = {country};
stochastic_output_campaign_ia.country_name(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = {country_name};
stochastic_output_campaign_ia.cohort_size(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = repmat(reshape(vimc_pop_mat(:,:,c)',101*101,1),nruns,1);

stochastic_output_camproutine.run_id = repmat(reshape(repmat(1:nruns,101*101,1),nruns*101*101,1),ncountry,1);
stochastic_output_camproutine.year =reshape(repmat(2000:2100,1,nruns*ncountry*101),1,nruns*ncountry*101*101)';
stochastic_output_camproutine.age = reshape(repmat(0:100,101,nruns*ncountry),nruns*ncountry*101*101,1);
stochastic_output_camproutine.country(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = {country};
stochastic_output_camproutine.country_name(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = {country_name};
stochastic_output_camproutine.cohort_size(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = repmat(reshape(vimc_pop_mat(:,:,c)',101*101,1),nruns,1);

stochastic_output_camproutine_ia.run_id = repmat(reshape(repmat(1:nruns,101*101,1),nruns*101*101,1),ncountry,1);
stochastic_output_camproutine_ia.year =reshape(repmat(2000:2100,1,nruns*ncountry*101),1,nruns*ncountry*101*101)';
stochastic_output_camproutine_ia.age = reshape(repmat(0:100,101,nruns*ncountry),nruns*ncountry*101*101,1);
stochastic_output_camproutine_ia.country(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = {country};
stochastic_output_camproutine_ia.country_name(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = {country_name};
stochastic_output_camproutine_ia.cohort_size(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = repmat(reshape(vimc_pop_mat(:,:,c)',101*101,1),nruns,1);

%central output tables
central_output_novacc.disease(1:ncountry*101*101,1) = {'Typhoid'};
central_output_novacc.year = reshape(repmat(2000:2100,1,ncountry*101),ncountry*101*101,1);
central_output_novacc.age = reshape(repmat(0:100,101,ncountry),ncountry*101*101,1);
central_output_novacc.country(101*101*(c-1)+1:101*101*c,1) = {country};
central_output_novacc.country_name(101*101*(c-1)+1:101*101*c,1) = {country_name};
central_output_novacc.cohort_size(101*101*(c-1)+1:101*101*c,1) = reshape(vimc_pop_mat(:,:,c)',101*101,1);

central_output_campaign.disease(1:ncountry*101*101,1) = {'Typhoid'};
central_output_campaign.year = reshape(repmat(2000:2100,1,ncountry*101),ncountry*101*101,1);
central_output_campaign.age = reshape(repmat(0:100,101,ncountry),ncountry*101*101,1);
central_output_campaign.country(101*101*(c-1)+1:101*101*c,1) = {country};
central_output_campaign.country_name(101*101*(c-1)+1:101*101*c,1) = {country_name};
central_output_campaign.cohort_size(101*101*(c-1)+1:101*101*c,1) = reshape(vimc_pop_mat(:,:,c)',101*101,1);

central_output_campaign_ia.disease(1:ncountry*101*101,1) = {'Typhoid'};
central_output_campaign_ia.year = reshape(repmat(2000:2100,1,ncountry*101),ncountry*101*101,1);
central_output_campaign_ia.age = reshape(repmat(0:100,101,ncountry),ncountry*101*101,1);
central_output_campaign_ia.country(101*101*(c-1)+1:101*101*c,1) = {country};
central_output_campaign_ia.country_name(101*101*(c-1)+1:101*101*c,1) = {country_name};
central_output_campaign_ia.cohort_size(101*101*(c-1)+1:101*101*c,1) = reshape(vimc_pop_mat(:,:,c)',101*101,1);

central_output_camproutine.disease(1:ncountry*101*101,1) = {'Typhoid'};
central_output_camproutine.year = reshape(repmat(2000:2100,1,ncountry*101),ncountry*101*101,1);
central_output_camproutine.age = reshape(repmat(0:100,101,ncountry),ncountry*101*101,1);
central_output_camproutine.country(101*101*(c-1)+1:101*101*c,1) = {country};
central_output_camproutine.country_name(101*101*(c-1)+1:101*101*c,1) = {country_name};
central_output_camproutine.cohort_size(101*101*(c-1)+1:101*101*c,1) = reshape(vimc_pop_mat(:,:,c)',101*101,1);

central_output_camproutine_ia.disease(1:ncountry*101*101,1) = {'Typhoid'};
central_output_camproutine_ia.year = reshape(repmat(2000:2100,1,ncountry*101),ncountry*101*101,1);
central_output_camproutine_ia.age = reshape(repmat(0:100,101,ncountry),ncountry*101*101,1);
central_output_camproutine_ia.country(101*101*(c-1)+1:101*101*c,1) = {country};
central_output_camproutine_ia.country_name(101*101*(c-1)+1:101*101*c,1) = {country_name};
central_output_camproutine_ia.cohort_size(101*101*(c-1)+1:101*101*c,1) = reshape(vimc_pop_mat(:,:,c)',101*101,1);


%for c=1:ncountry
for i =1:size(output2,1)  
    % Read in cases, deaths, dalys for run_id ==1
    % no vacc
    %stochastic_output_novacc.cases(logical((stochastic_output_novacc.run_id==i).*strcmp(stochastic_output_novacc.country,country)),1) = reshape(output2{i,c}{1,1}.cases,101*101,1);
    %stochastic_output_novacc.deaths(logical((stochastic_output_novacc.run_id==i).*strcmp(stochastic_output_novacc.country,country)),1) = reshape(output2{i,c}{1,1}.deaths,101*101,1);
    %stochastic_output_novacc.dalys(logical((stochastic_output_novacc.run_id==i).*strcmp(stochastic_output_novacc.country,country)),1) = reshape(output2{i,c}{1,1}.dalys,101*101,1);
    stochastic_output_novacc.cases(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2{i,c}{1,1}.cases,101*101,1);
    stochastic_output_novacc.deaths(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2{i,c}{1,1}.deaths,101*101,1);
    stochastic_output_novacc.dalys(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2{i,c}{1,1}.dalys,101*101,1);

    % Campaign only
    stochastic_output_campaign.cases(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2{i,c}{2,1}.cases,101*101,1);
    stochastic_output_campaign.deaths(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2{i,c}{2,1}.deaths,101*101,1);
    stochastic_output_campaign.dalys(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2{i,c}{2,1}.dalys,101*101,1);
    
    % Campagin only -- IA2030
    stochastic_output_campaign_ia.cases(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2{i,c}{3,1}.cases,101*101,1);
    stochastic_output_campaign_ia.deaths(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2{i,c}{3,1}.deaths,101*101,1);
    stochastic_output_campaign_ia.dalys(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2{i,c}{3,1}.dalys,101*101,1);

    % Campaign & Routine
    stochastic_output_camproutine.cases(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2{i,c}{4,1}.cases,101*101,1);
    stochastic_output_camproutine.deaths(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2{i,c}{4,1}.deaths,101*101,1);
    stochastic_output_camproutine.dalys(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2{i,c}{4,1}.dalys,101*101,1);
    
    % Campaign & Routine -- IA2030
    stochastic_output_camproutine_ia.cases(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2{i,c}{5,1}.cases,101*101,1);
    stochastic_output_camproutine_ia.deaths(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2{i,c}{5,1}.deaths,101*101,1);
    stochastic_output_camproutine_ia.dalys(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2{i,c}{5,1}.dalys,101*101,1);

end

% generate means for central output
central_output_novacc.cases(101*101*(c-1)+1:101*101*c,1) = reshape(mean(reshape(stochastic_output_novacc.cases(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1);
central_output_novacc.deaths(101*101*(c-1)+1:101*101*c,1) = reshape(mean(reshape(stochastic_output_novacc.deaths(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1);
central_output_novacc.dalys(101*101*(c-1)+1:101*101*c,1) = reshape(mean(reshape(stochastic_output_novacc.dalys(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1);

central_output_campaign.cases(101*101*(c-1)+1:101*101*c,1) = reshape(mean(reshape(stochastic_output_campaign.cases(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1);
central_output_campaign.deaths(101*101*(c-1)+1:101*101*c,1) = reshape(mean(reshape(stochastic_output_campaign.deaths(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1);
central_output_campaign.dalys(101*101*(c-1)+1:101*101*c,1) = reshape(mean(reshape(stochastic_output_campaign.dalys(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1);

central_output_campaign_ia.cases(101*101*(c-1)+1:101*101*c,1) = reshape(mean(reshape(stochastic_output_campaign_ia.cases(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1);
central_output_campaign_ia.deaths(101*101*(c-1)+1:101*101*c,1) = reshape(mean(reshape(stochastic_output_campaign_ia.deaths(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1);
central_output_campaign_ia.dalys(101*101*(c-1)+1:101*101*c,1) = reshape(mean(reshape(stochastic_output_campaign_ia.dalys(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1);

central_output_camproutine.cases(101*101*(c-1)+1:101*101*c,1) = reshape(mean(reshape(stochastic_output_camproutine.cases(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1);
central_output_camproutine.deaths(101*101*(c-1)+1:101*101*c,1) = reshape(mean(reshape(stochastic_output_camproutine.deaths(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1);
central_output_camproutine.dalys(101*101*(c-1)+1:101*101*c,1) = reshape(mean(reshape(stochastic_output_camproutine.dalys(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1);

central_output_camproutine_ia.cases(101*101*(c-1)+1:101*101*c,1) = reshape(mean(reshape(stochastic_output_camproutine_ia.cases(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1);
central_output_camproutine_ia.deaths(101*101*(c-1)+1:101*101*c,1) = reshape(mean(reshape(stochastic_output_camproutine_ia.deaths(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1);
central_output_camproutine_ia.dalys(101*101*(c-1)+1:101*101*c,1) = reshape(mean(reshape(stochastic_output_camproutine_ia.dalys(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1);

end

%convert stochastic outputs to tables
stochastic_output_novacc_table = struct2table(stochastic_output_novacc);
stochastic_output_campaign_table = struct2table(stochastic_output_campaign);
stochastic_output_campaign_table_ia = struct2table(stochastic_output_campaign);
stochastic_output_camproutine_table = struct2table(stochastic_output_camproutine);
stochastic_output_camproutine_table_ia = struct2table(stochastic_output_camproutine);

%also central outputs for IA2030
central_output_campaign_table_ia = struct2table(central_output_campaign_ia);
central_output_camproutine_table_ia = struct2table(central_output_camproutine_ia);

%mean_total_cases_novacc = median(squeeze(sum(reshape(stochastic_output_novacc.cases,101,101,nruns),2)),2);
%mean_total_deaths_novacc = median(squeeze(sum(reshape(stochastic_output_novacc.deaths,101,101,nruns),2)),2);
%mean_total_dalys_novacc = median(squeeze(sum(reshape(stochastic_output_novacc.dalys,101,101,nruns),2)),2);

%mean_total_cases_default = median(squeeze(sum(reshape(stochastic_output_default.cases,101,101,nruns),2)),2);
%mean_total_deaths_default = median(squeeze(sum(reshape(stochastic_output_default.deaths,101,101,nruns),2)),2);
%mean_total_dalys_default = median(squeeze(sum(reshape(stochastic_output_default.dalys,101,101,nruns),2)),2);

%total_output_novacc = table([2000:2100]',mean_total_cases_novacc, mean_total_deaths_novacc, mean_total_dalys_novacc,...
%    'VariableName',{'Year','MeanTotalCases','MeanTotalDeaths','MeanTotalDalys'});
%total_output_default = table([2000:2100]',mean_total_cases_default, mean_total_deaths_default, mean_total_dalys_default,...
%    'VariableName',{'Year','MeanTotalCases','MeanTotalDeaths','MeanTotalDalys'});

%NO VAX
    writetable(stochastic_output_novacc_table,'stochastic_output_novacc_202110_final.csv')
    
    writetable(central_output_novacc,'central_output_novacc_202110_final.csv')

%DEFAULT
    writetable(stochastic_output_campaign_table,'stochastic_output_campaign-DEFAULT_202110_final.csv')
    writetable(stochastic_output_camproutine_table,'stochastic_output_camproutine-DEFAULT_202110_final.csv')

    writetable(central_output_campaign,'central_output_campaign-DEFAULT_202110_final.csv')
    writetable(central_output_camproutine,'central_output_camproutine-DEFAULT_202110_final.csv') %forgot to change this one to "test"

%IA2030
    %writetable(stochastic_output_novacc_table,'stochastic_output_novacc_202110_IA2030_test.csv')
    writetable(stochastic_output_campaign_table_ia,'stochastic_output_campaign_202110_IA2030_final.csv')
    writetable(stochastic_output_camproutine_table_ia,'stochastic_output_camproutine_202110_IA2030_final.csv')

    writetable(central_output_campaign_table_ia,'central_output_campaign_202110_IA2030_final.csv')
    writetable(central_output_camproutine_table_ia,'central_output_camproutine_202110_IA2030_final.csv')



%writetable(total_output_novacc,'total_output_novacc.csv')
%writetable(total_output_default,'total_output_default.csv')

%% Compile parameter file to export

%pull out common estimates for phosp and pdeathhosp
phosp_common=phosp(:,1);
%pdeathhosp_common=pdeathhosp(:,1);

%make vectors of countries with specific phosp and pdeathhosp estimates
phosp_countries=["BGD", "IND", "KEN", "IDN", "VNM", "PAK"]'; %countries with specific phosp esimates

param_table = table; %'Size',[nruns,27],...

param_table.run_id = (1:nruns)';
for i = 1:nruns
    for c=1:ncountry
        z = find(strcmp(all_countries.ISO3,vimc_countries{c}));
        eval(strcat('param_table.R0_',vimc_countries{c},'(i)=','R0samples(i,z);'));
        eval(strcat('param_table.rep_',vimc_countries{c},'(i)=','repsamples(i,z);'));
        eval(strcat('param_table.pdeathhosp_',vimc_countries{c},'(i)=','pdeathhosp(i,c);'));
       
    end
    for p=1:length(phosp_countries)
        phosp_indx = find(strcmp(vimc_countries,phosp_countries(p)));
        eval(strcat('param_table.phosp_',phosp_countries(p),'(i)=','phosp(i,phosp_indx);'));
    end
    %OPTION 1 -- report pdeathhosp by region
    %for h=1:length(pdeathhosp_regions)
       % pdeathhosp_indx = find(strcmp(country_regions,pdeathhosp_regions(h)));
       % eval(strcat('param_table.pdeathhosp_',pdeathhosp_regions(h),'(i)=','pdeathhosp(i,pdeathhosp_indx);'));
    %end
    %OPTION 2 -- report pdeathhosp by country (for countries that have a
    %regin-specific estimate)
    %for d=1:length(vimc_countries_subset) %STILL NEED TO CREATE THIS SUBSETTED DATA
       % pdeathhosp_indx = find(strcmp(vimc_countries,vimc_countries_subset.GBD_region(d)));
       % eval(strcat('param_table.pdeathhosp_',vimc_countries(d),'(i)=','pdeathhosp(i,pdeathhosp_indx);'));
   % end
    
        
    %read in parameters that are in output2 structure
    %param_table.R0(i) = output2{i,1}{1,1}.R0;
    param_table.phosp_common(i)=phosp_common(i);
    %param_table.pdeathhosp_common(i)=pdeathhosp_common(i);
    %param.table.pdeathhosp_GBDRegion1=pdeathhosp_R1;
    %param.table.pdeathhosp_GBDRegion4=pdeathhosp_R4;
    %param.table.pdeathhosp_GBDRegion5=pdeathhosp_R5;
    param_table.mult1(i) = output2{i,c}{1,1}.mult1;
    param_table.mult2(i) = output2{i,c}{1,1}.mult2;
    param_table.rC(i) = output2{i,c}{1,1}.rC;
    param_table.h2o(i) = output2{i,c}{1,1}.h2o;
    param_table.veff(i) = output2{i,c}{1,1}.veff;
    param_table.omega_v(i) = vax_samp.omega(i);
  
    
    %param_table.delta(i) = params.delta;
    %param_table.alpha(i) = params.alpha;
    %param_table.omega(i) = params.omega;
    %param_table.mub(i) = mub(agedist_country(gavi73.cn(z))); %birth rate
    %param_table.mub(i) = mub(agedist_country(all_countries.countryn(z))); %birth rate
    %param_table.theta_y(i) = agepar.theta(1);
    %param_table.theta_o(i) = agepar.theta(end);
    %param_table.r(i) = 1;
    %param_table.leb(i)= leb;
    for str2 = {'propAMR','rAMR','pmedcare','propdeathhosp',...
            'doi_hosp','doi_outp','doi_nomed','dwt_hosp','dwt_outp','dwt_nomed'}
        eval(strcat('param_table.',str2{:},'(i)=',str2{:},'(i);'));
    end
    
end

writetable(param_table,'stoch_param_table_final.csv');


