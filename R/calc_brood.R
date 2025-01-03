#==================================================

library(dplyr)


#explore RDS and map objects to the Cohort Analysis rates tab on the HAR_CalclationofStockandRecruitment_t02017.xlsx



har_ream<-readRDS("data/har2024msf.RDS") # need msf version because of refactoring



harCY<-har_ream$allcm1stock_data$CY_era
names(harCY)



by_seq<-harCY$FirstBY:harCY$LastBY

harCY$GearHR.BY[[which(by_seq==1981)]]
harCY$mortality_Ocean_byAge
harCY$AdltEqv #column AEQ rate

harCY$Flg1
#this is how to calculate the OcnHR_er
TotOcnCatch <- stock_data$mortality_byAge[[MortType]][[Size]][BY, j] - stock_data$mortality_Mature_byAge[[MortType]][[Size]][BY, j]
TotOcnCatch <- TotOcnCatch * stock_data$AdltEqv[BY, j]
stock_data$OcnHR_er[ZZ, j] <- stock_data$mortality_Ocean_byAge[ZZ, j] / stock_data$TotRunAtAge[j, SubScrpt]



harCY$TempHR[, 1, 2,2:5 ]

harCY$TempHR[1, 3, 2, 2:5]

harCY$mortality_Ocean_byAge

apply(harCY$HR$kept,2,sum)

har_ream$allcm1stock_data$CY_era$OcnHR_er


#HRJ numbers
harCY$terminal
harCY$escStray$Canada$reported

harCY$AdltEqv_Mortality_Kept

harCY$AdltEqv_Mortality_Total

harCY$ReducedCohort #[BY, i, Age]
#cohort after natural mortality
harCY$ReducedCohort[1,,2:5][1,]
#terminal run
harCY$ReducedCohort[1,,2:5][7,]

harCY$AdltEqv
names(harCY)

names(harCY)[grep("Mat", names(harCY),perl = TRUE)]

harCY$HR_BY[["kept"]][[1]][1:79,][!harCY$terminal]

#PretermER
preterminal1981<-replace(harCY$HR_BY[["kept"]][[1]][1:79,],harCY$terminal,0)+replace(harCY$HR_BY[["IM"]][[1]][1:79,],harCY$terminal,0)
apply(preterminal1981,2,sum, na.rm=T)
apply(preterminal1981[apply(preterminal1981,1,sum, na.rm=T)>0&preterminal1981[,2]<0.001&preterminal1981[,3]<0.001,],2,sum)

#TotTermHR
terminal1981<-replace(harCY$HR_BY[["kept"]][[1]][1:79,],!harCY$terminal,0)+replace(harCY$HR_BY[["IM"]][[1]][1:79,],!harCY$terminal,0)
apply(terminal1981,2,sum, na.rm=T)



pretermnet<-which(harCY$fishName%in%harCY$fishName[grep("^(?!T).*(N)$", harCY$fishName,perl = TRUE)])
matnetages<-harCY$TermNetSwitchAge:harCY$MaxStkAge

mature_ocean_net <- matrix(0, ncol=harCY$MaxStkAge, nrow=length(by_seq))

preterminal_er <- matrix(0, ncol=harCY$MaxStkAge, nrow=length(by_seq))

#Ocean ER on immature




for(i in seq_along(by_seq)){
  mature_ocean_net[i,matnetages] <- apply( harCY$HR_BY[["kept"]][[i]][pretermnet,matnetages]+
                                harCY$HR_BY[["IM"]][[i]][pretermnet,matnetages],2,sum)

  preterminal_er[i,] <- apply(replace(harCY$HR_BY[["kept"]][[i]][1:79,],harCY$terminal,0)+
                    replace(harCY$HR_BY[["IM"]][[i]][1:79,],harCY$terminal,0),2,sum)[1:harCY$MaxStkAge]
}

maturation_rate_long <- reshape2::melt(harCY$MatRte)|>
                        rename( "brood"=Var1, "age"=Var2, "maturation_rate"=value )|>
                        filter(age>1)
                        

mature_ocean_net_long<-reshape2::melt(mature_ocean_net)|>
                        rename( "brood"=Var1, "age"=Var2, "mature_ocean_net_hr"=value )|>
                        mutate(brood = brood+harCY$FirstBY-1)|>
                        filter(age>1)


preterminal_er_long <- reshape2::melt(preterminal_er)|>
                        rename( "brood"=Var1, "age"=Var2, "preterminal_er"=value )|>
                        mutate(brood = brood+harCY$FirstBY-1)|>
                        filter(age>1)


aeq_long<-reshape2::melt(harCY$AdltEqv)|>
                        rename( "brood"=Var1, "age"=Var2, "aeq"=value )|>
                        filter(age>1&age<6)

preterminal_er_long$natural_survival<-harCY$SurvRte[preterminal_er_long$age-1]


head(maturation_rate_long)
#extract REAM objects

names(harCY)[grep("Age", names(harCY))]
1-harCY$SurvRte

harCY$StartAge
MatureOceanNetER1981<-harCY$HR_BY[["kept"]][[1]][pretermnet,]+harCY$HR_BY[["IM"]][[1]][pretermnet,]
apply(MatureOceanNetER1981,2,sum)

MatureOceanNetER1982<-harCY$HR_BY[["kept"]][[2]][pretermnet,]+harCY$HR_BY[["IM"]][[2]][pretermnet,]
apply(MatureOceanNetER1982,2,sum)


str_detect(harCY$fishName,)
harCY$fishName[c(6,18,50,51)]
#need to compute MatureOceanNetER

dim(harCY$HR_BY[["kept"]][[1]][1:79,])

###These are the ERs on the file!!!!
apply(harCY$HR_BY[["kept"]][[1]],2,sum)+apply(harCY$HR_BY[["IM"]][[1]],2,sum)
#need to split into:
#pre terminal ER
#mature_ocean_net_er
#true_term_er
#tot_term_hr =true_term_er+mature_ocean_net_er

#read isn river spawners at age -- these numbers are originaly from HAR calculation spreadsheer tab River Spawners
#It needs updating
river_spawners<-read.csv("data/river_spawners_1981_2015.csv")
head(river_spawners)
dim(river_spawners)

#brood removals
brood_removals<-read.csv("data/brood_removals_1981_2015.csv")
head(brood_removals)

#terminal catch
terminal_catch<-read.csv("data/terminal_catch_1982_2015.csv")
head(terminal_catch)

#hatchery_origin_esc
hatchery_origin_esc<-read.csv("data/hatchery_origin_esc_1981_2015.csv")

female_spawners <-read.csv("data/female_spawners.csv")
female_spawners$brood <- female_spawners$year-female_spawners$age
female_spawners$lookup <- as.numeric(paste0(female_spawners$brood,female_spawners$age))


hatchery_origin_esc$brood <- hatchery_origin_esc$year-hatchery_origin_esc$age
hatchery_origin_esc$lookup <- as.numeric(paste0(hatchery_origin_esc$brood,hatchery_origin_esc$age))
hatchery_origin_esc$cwt_exp_total <- rowSums(hatchery_origin_esc[,c("cwt_exp_rs",
                                       "cwt_exp_br", 
                                       "cwt_exp_additional")]
                                       ,na.rm=T)




all_escapement_data <- full_join(river_spawners,brood_removals) |>
                      full_join(terminal_catch)|>
                      full_join(hatchery_origin_esc)|>
                      full_join(female_spawners)|>
                      full_join(mature_ocean_net_long)|>
                      full_join(maturation_rate_long) |>
                      full_join(preterminal_er_long)|>
                      full_join( aeq_long)


#Calculated quantities
all_escapement_data$return_river <- rowSums(all_escapement_data[,c("river_spawners",
                                       "brood_removals")]
                                       ,na.rm=T)

all_escapement_data$return_river_natural <- ifelse(all_escapement_data$return_river-all_escapement_data$cwt_exp_total>1,
                                                   all_escapement_data$return_river-all_escapement_data$cwt_exp_total,
                                                   100)


#calculate proportion wild
all_escapement_data$proportion_wild <- ifelse(all_escapement_data$return_river>0,
                                             1-all_escapement_data$cwt_exp_total/all_escapement_data$return_river,
                                             1)


all_escapement_data$terminal_catch_natural <- ifelse(all_escapement_data$terminal_catch*all_escapement_data$proportion_wild>0,
                                                     all_escapement_data$terminal_catch*all_escapement_data$proportion_wild,
                                                     0)
                                                       
all_escapement_data$terminal_run_natural <- rowSums(all_escapement_data[,c("return_river_natural",
                                       "terminal_catch_natural")]
                                       ,na.rm=T)


all_escapement_data$total_mature_run <- ifelse(1-all_escapement_data$mature_ocean_net_hr>0,
                                 all_escapement_data$terminal_run_natural/(1-all_escapement_data$mature_ocean_net_hr),
                                 0)

all_escapement_data$ocean_post_fishery_abundance <- all_escapement_data$total_mature_run/all_escapement_data$maturation_rate

all_escapement_data$surviving_cohort <- all_escapement_data$ocean_post_fishery_abundance-all_escapement_data$total_mature_run


all_escapement_data$ocean_pre_fishery_abundance <- ifelse((1-all_escapement_data$preterminal_er)>0,
                                                   all_escapement_data$ocean_post_fishery_abundance/(1-all_escapement_data$preterminal_er),
                                                   0)

all_escapement_data$cohort_abundance <-all_escapement_data$ocean_pre_fishery_abundance/(all_escapement_data$natural_survival)


all_escapement_data$ocean_mortality_immature<-all_escapement_data$ocean_pre_fishery_abundance-
                                           all_escapement_data$ocean_post_fishery_abundance


all_escapement_data$ocean_catch_aeq<-all_escapement_data$ocean_mortality_immature*all_escapement_data$aeq

all_escapement_data$adult_recruits<-all_escapement_data$ocean_catch_aeq+all_escapement_data$total_mature_run

all_escapement_data$maturity_rate_natural_population <- ifelse(all_escapement_data$ocean_post_fishery_abundance>0,
                                                           all_escapement_data$total_mature_run/all_escapement_data$ocean_post_fishery_abundance,
                                                           1)


#recruits 
recruits<-reshape2::dcast(all_escapement_data, brood ~ age, value.var = 'adult_recruits')

female_spawner<-reshape2::dcast(all_escapement_data, year ~ age, value.var = 'female_spawner',fun.aggregate =sum)
female_spawner$female_spawners_3_5<-apply(female_spawner[,3:5],1,sum)

all_spawners<-reshape2::dcast(river_spawners, year ~ age, value.var = 'river_spawners',fun.aggregate =sum)
all_spawners$all_spawners_3_5<-apply(all_spawners[,3:5],1,sum)

spawners_prop<-left_join(female_spawner[,c("year","female_spawners_3_5")],all_spawners[,c("year","all_spawners_3_5")])

spawners_prop$female_rate<-spawners_prop$female_spawners_3_5/spawners_prop$all_spawners_3_5

all_escapement_data$adjustment_escapement <- all_escapement_data$river_spawners/mean(spawners_prop$female_rate,na.rm=T)

