%% Generate output tables formatted for VIMC

% For stochastic outputs, one for baseline (no vacc), one for camp 15 (default):
%for each run id, need year in order, then age group, cohort size, cases, dalys

% stochastic output tables
stochastic_output_novacc = table; 
stochastic_output_campaign = table;
stochastic_output_camproutine = table;
stochastic_output_campaign_bs = table;
stochastic_output_camproutine_bs = table;

% For central output: same as stochastic, but mean values
central_output_novacc = table; 
central_output_campaign = table;
central_output_camproutine = table;
central_output_campaign_bs = table;
central_output_camproutine_bs = table;

for c=1:ncountry
    
country=vimc_countries{c}; 
country_name=vimc_country_names{c};   
    
% read in disease, years, age, country, cohort size, and then set up default tables as copies
stochastic_output_novacc.disease(1:nruns*ncountry*101*101,1) = {'Typhoid'};
stochastic_output_novacc.run_id = repmat(reshape(repmat(1:nruns,101*101,1),nruns*101*101,1),ncountry,1);
stochastic_output_novacc.year =reshape(repmat(2000:2100,1,nruns*ncountry*101),1,nruns*ncountry*101*101)';
stochastic_output_novacc.age = reshape(repmat(0:100,101,nruns*ncountry),nruns*ncountry*101*101,1);
stochastic_output_novacc.country(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = {country};
stochastic_output_novacc.country_name(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = {country_name};
stochastic_output_novacc.cohort_size(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = repmat(reshape(vimc_pop_mat(:,:,c)',101*101,1),nruns,1);

stochastic_output_campaign.disease(1:nruns*ncountry*101*101,1) = {'Typhoid'};
stochastic_output_campaign.run_id = repmat(reshape(repmat(1:nruns,101*101,1),nruns*101*101,1),ncountry,1);
stochastic_output_campaign.year =reshape(repmat(2000:2100,1,nruns*ncountry*101),1,nruns*ncountry*101*101)';
stochastic_output_campaign.age = reshape(repmat(0:100,101,nruns*ncountry),nruns*ncountry*101*101,1);
stochastic_output_campaign.country(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = {country};
stochastic_output_campaign.country_name(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = {country_name};
stochastic_output_campaign.cohort_size(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = repmat(reshape(vimc_pop_mat(:,:,c)',101*101,1),nruns,1);

stochastic_output_campaign_bs.disease(1:nruns*ncountry*101*101,1) = {'Typhoid'};
stochastic_output_campaign_bs.run_id = repmat(reshape(repmat(1:nruns,101*101,1),nruns*101*101,1),ncountry,1);
stochastic_output_campaign_bs.year =reshape(repmat(2000:2100,1,nruns*ncountry*101),1,nruns*ncountry*101*101)';
stochastic_output_campaign_bs.age = reshape(repmat(0:100,101,nruns*ncountry),nruns*ncountry*101*101,1);
stochastic_output_campaign_bs.country(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = {country};
stochastic_output_campaign_bs.country_name(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = {country_name};
stochastic_output_campaign_bs.cohort_size(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = repmat(reshape(vimc_pop_mat(:,:,c)',101*101,1),nruns,1);

stochastic_output_camproutine.disease(1:nruns*ncountry*101*101,1) = {'Typhoid'};
stochastic_output_camproutine.run_id = repmat(reshape(repmat(1:nruns,101*101,1),nruns*101*101,1),ncountry,1);
stochastic_output_camproutine.year =reshape(repmat(2000:2100,1,nruns*ncountry*101),1,nruns*ncountry*101*101)';
stochastic_output_camproutine.age = reshape(repmat(0:100,101,nruns*ncountry),nruns*ncountry*101*101,1);
stochastic_output_camproutine.country(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = {country};
stochastic_output_camproutine.country_name(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = {country_name};
stochastic_output_camproutine.cohort_size(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = repmat(reshape(vimc_pop_mat(:,:,c)',101*101,1),nruns,1);

stochastic_output_camproutine_bs.disease(1:nruns*ncountry*101*101,1) = {'Typhoid'};
stochastic_output_camproutine_bs.run_id = repmat(reshape(repmat(1:nruns,101*101,1),nruns*101*101,1),ncountry,1);
stochastic_output_camproutine_bs.year =reshape(repmat(2000:2100,1,nruns*ncountry*101),1,nruns*ncountry*101*101)';
stochastic_output_camproutine_bs.age = reshape(repmat(0:100,101,nruns*ncountry),nruns*ncountry*101*101,1);
stochastic_output_camproutine_bs.country(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = {country};
stochastic_output_camproutine_bs.country_name(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = {country_name};
stochastic_output_camproutine_bs.cohort_size(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = repmat(reshape(vimc_pop_mat(:,:,c)',101*101,1),nruns,1);

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

central_output_campaign_bs.disease(1:ncountry*101*101,1) = {'Typhoid'};
central_output_campaign_bs.year = reshape(repmat(2000:2100,1,ncountry*101),ncountry*101*101,1);
central_output_campaign_bs.age = reshape(repmat(0:100,101,ncountry),ncountry*101*101,1);
central_output_campaign_bs.country(101*101*(c-1)+1:101*101*c,1) = {country};
central_output_campaign_bs.country_name(101*101*(c-1)+1:101*101*c,1) = {country_name};
central_output_campaign_bs.cohort_size(101*101*(c-1)+1:101*101*c,1) = reshape(vimc_pop_mat(:,:,c)',101*101,1);

central_output_camproutine.disease(1:ncountry*101*101,1) = {'Typhoid'};
central_output_camproutine.year = reshape(repmat(2000:2100,1,ncountry*101),ncountry*101*101,1);
central_output_camproutine.age = reshape(repmat(0:100,101,ncountry),ncountry*101*101,1);
central_output_camproutine.country(101*101*(c-1)+1:101*101*c,1) = {country};
central_output_camproutine.country_name(101*101*(c-1)+1:101*101*c,1) = {country_name};
central_output_camproutine.cohort_size(101*101*(c-1)+1:101*101*c,1) = reshape(vimc_pop_mat(:,:,c)',101*101,1);

central_output_camproutine_bs.disease(1:ncountry*101*101,1) = {'Typhoid'};
central_output_camproutine_bs.year = reshape(repmat(2000:2100,1,ncountry*101),ncountry*101*101,1);
central_output_camproutine_bs.age = reshape(repmat(0:100,101,ncountry),ncountry*101*101,1);
central_output_camproutine_bs.country(101*101*(c-1)+1:101*101*c,1) = {country};
central_output_camproutine_bs.country_name(101*101*(c-1)+1:101*101*c,1) = {country_name};
central_output_camproutine_bs.cohort_size(101*101*(c-1)+1:101*101*c,1) = reshape(vimc_pop_mat(:,:,c)',101*101,1);


for i =1:size(output2023,1)  
    % Read in cases, deaths, dalys for run_id ==1
    % no vacc
    %stochastic_output_novacc.cases(logical((stochastic_output_novacc.run_id==i).*strcmp(stochastic_output_novacc.country,country)),1) = reshape(output2{i,c}{1,1}.cases,101*101,1);
    %stochastic_output_novacc.deaths(logical((stochastic_output_novacc.run_id==i).*strcmp(stochastic_output_novacc.country,country)),1) = reshape(output2{i,c}{1,1}.deaths,101*101,1);
    %stochastic_output_novacc.dalys(logical((stochastic_output_novacc.run_id==i).*strcmp(stochastic_output_novacc.country,country)),1) = reshape(output2{i,c}{1,1}.dalys,101*101,1);
    stochastic_output_novacc.cases(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2023{i,c}{1,1}.cases,101*101,1);
    stochastic_output_novacc.deaths(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2023{i,c}{1,1}.deaths,101*101,1);
    stochastic_output_novacc.dalys(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2023{i,c}{1,1}.dalys,101*101,1);

    % Campaign only
    stochastic_output_campaign.cases(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2023{i,c}{2,1}.cases,101*101,1);
    stochastic_output_campaign.deaths(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2023{i,c}{2,1}.deaths,101*101,1);
    stochastic_output_campaign.dalys(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2023{i,c}{2,1}.dalys,101*101,1);
    
    % Campaign only -- Blue-sky
    stochastic_output_campaign_bs.cases(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2023{i,c}{3,1}.cases,101*101,1);
    stochastic_output_campaign_bs.deaths(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2023{i,c}{3,1}.deaths,101*101,1);
    stochastic_output_campaign_bs.dalys(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2023{i,c}{3,1}.dalys,101*101,1);

    % Campaign & Routine
    stochastic_output_camproutine.cases(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2023{i,c}{4,1}.cases,101*101,1);
    stochastic_output_camproutine.deaths(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2023{i,c}{4,1}.deaths,101*101,1);
    stochastic_output_camproutine.dalys(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2023{i,c}{4,1}.dalys,101*101,1);
    
    % Campaign & Routine -- Blue-sky
    stochastic_output_camproutine_bs.cases(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2023{i,c}{5,1}.cases,101*101,1);
    stochastic_output_camproutine_bs.deaths(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2023{i,c}{5,1}.deaths,101*101,1);
    stochastic_output_camproutine_bs.dalys(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2023{i,c}{5,1}.dalys,101*101,1);

end

% generate means for central output
central_output_novacc.cases(101*101*(c-1)+1:101*101*c,1) = round(reshape(mean(reshape(stochastic_output_novacc.cases(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1),2);
central_output_novacc.deaths(101*101*(c-1)+1:101*101*c,1) = round(reshape(mean(reshape(stochastic_output_novacc.deaths(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1),2);
central_output_novacc.dalys(101*101*(c-1)+1:101*101*c,1) = round(reshape(mean(reshape(stochastic_output_novacc.dalys(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1),2);

central_output_campaign.cases(101*101*(c-1)+1:101*101*c,1) = round(reshape(mean(reshape(stochastic_output_campaign.cases(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1),2);
central_output_campaign.deaths(101*101*(c-1)+1:101*101*c,1) = round(reshape(mean(reshape(stochastic_output_campaign.deaths(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1),2);
central_output_campaign.dalys(101*101*(c-1)+1:101*101*c,1) = round(reshape(mean(reshape(stochastic_output_campaign.dalys(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1),2);

central_output_campaign_bs.cases(101*101*(c-1)+1:101*101*c,1) = round(reshape(mean(reshape(stochastic_output_campaign_bs.cases(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1),2);
central_output_campaign_bs.deaths(101*101*(c-1)+1:101*101*c,1) = round(reshape(mean(reshape(stochastic_output_campaign_bs.deaths(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1),2);
central_output_campaign_bs.dalys(101*101*(c-1)+1:101*101*c,1) = round(reshape(mean(reshape(stochastic_output_campaign_bs.dalys(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1),2);

central_output_camproutine.cases(101*101*(c-1)+1:101*101*c,1) = round(reshape(mean(reshape(stochastic_output_camproutine.cases(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1),2);
central_output_camproutine.deaths(101*101*(c-1)+1:101*101*c,1) = round(reshape(mean(reshape(stochastic_output_camproutine.deaths(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1),2);
central_output_camproutine.dalys(101*101*(c-1)+1:101*101*c,1) = round(reshape(mean(reshape(stochastic_output_camproutine.dalys(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1),2);

central_output_camproutine_bs.cases(101*101*(c-1)+1:101*101*c,1) = round(reshape(mean(reshape(stochastic_output_camproutine_bs.cases(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1),2);
central_output_camproutine_bs.deaths(101*101*(c-1)+1:101*101*c,1) = round(reshape(mean(reshape(stochastic_output_camproutine_bs.deaths(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1),2);
central_output_camproutine_bs.dalys(101*101*(c-1)+1:101*101*c,1) = round(reshape(mean(reshape(stochastic_output_camproutine_bs.dalys(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1),2);

end

%convert stochastic outputs to tables
%stochastic_output_novacc_table = round(struct2table(stochastic_output_novacc));
%stochastic_output_campaign_table = round(struct2table(stochastic_output_campaign));
%stochastic_output_campaign_table_bs = round(struct2table(stochastic_output_campaign));
%stochastic_output_camproutine_table = round(struct2table(stochastic_output_camproutine));
%stochastic_output_camproutine_table_bs = round(struct2table(stochastic_output_camproutine));

%also central outputs for Blue-sky
%central_output_campaign_table_bs = round(struct2table(central_output_campaign_bs),2);
%central_output_camproutine_table_bs = round(struct2table(central_output_camproutine_bs),2);

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

%% Write csv files

%NO VAX
writetable(stochastic_output_novacc,'typhoid-no-vaccination_stochastic_output2023_Yale.csv')
writetable(central_output_novacc,'typhoid-no-vaccination_central_output2023_Yale.csv')

%DEFAULT
writetable(stochastic_output_campaign,'typhoid-campaign-default_stochastic_output2023_Yale.csv')
writetable(stochastic_output_camproutine,'typhoid-campaign-routine-default_stochastic_output2023_Yale.csv')

writetable(central_output_campaign,'typhoid-campaign-default_central_output2023_Yale.csv')
writetable(central_output_camproutine,'typhoid-campaign-routine-default_central_output2023_Yale.csv') 

%Blue-sky
writetable(stochastic_output_campaign_bs,'typhoid-campaign-bluesky_stochastic_output2023_Yale.csv')
writetable(stochastic_output_camproutine_bs,'typhoid-campaign-routine-bluesky_stochastic_output2023_Yale.csv')

writetable(central_output_campaign_bs,'typhoid-campaign-bluesky_central_output2023_Yale.csv')
writetable(central_output_camproutine_bs,'typhoid-campaign-routine-bluesky_central_output2023_Yale.csv')



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
    param_table.mult1(i) = output2023{i,c}{1,1}.mult1;
    param_table.mult2(i) = output2023{i,c}{1,1}.mult2;
    param_table.rC(i) = output2023{i,c}{1,1}.rC;
    param_table.h2o(i) = output2023{i,c}{1,1}.h2o;
    param_table.veff(i) = output2023{i,c}{1,1}.veff;
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

writetable(param_table,'stoch_param_table_2023_Yale.csv');

