#==================================================
#load required libraries
library(dplyr)
library(reshape2)

#explore RDS and map objects to the Cohort Analysis rates tab on the HAR_CalclationofStockandRecruitment_t02017.xlsx


#read in required data




#ream output need msf version because of refactoring
ream<-readRDS("data/har2024msf.RDS") # 


#read in escapement data
#It needs updating
river_spawners<-read.csv("data/river_spawners_1981_2015.csv")

#brood removals
brood_removals<-read.csv("data/brood_removals_1981_2015.csv")

#terminal catch
terminal_catch<-read.csv("data/terminal_catch_1982_2015.csv")

#hatchery_origin_esc
hatchery_origin_esc<-read.csv("data/hatchery_origin_esc_1981_2015.csv")

#female spawners at age
female_spawners <-read.csv("data/female_spawners.csv")



#Compiling required data from the ERA analysis

outCY<-ream$allcm1stock_data$CY_era
by_seq<-outCY$FirstBY:outCY$LastBY


#HRJ numbers
#PretermER
#identify pre terminla fisheries 
pretermnet<-which(outCY$fishName%in%outCY$fishName[grep("^(?!T).*(N)$", outCY$fishName,perl = TRUE)])
matnetages<-outCY$TermNetSwitchAge:outCY$MaxStkAge

mature_ocean_net <- matrix(0, ncol=outCY$MaxStkAge, nrow=length(by_seq))

preterminal_er <- matrix(0, ncol=outCY$MaxStkAge, nrow=length(by_seq))

#Ocean ER on immature
for(i in seq_along(by_seq)){
  mature_ocean_net[i,matnetages] <- apply( outCY$HR_BY[["kept"]][[i]][pretermnet,matnetages]+
                                outCY$HR_BY[["IM"]][[i]][pretermnet,matnetages],2,sum)

  preterminal_er[i,] <- apply(replace(outCY$HR_BY[["kept"]][[i]][1:79,],outCY$terminal,0)+
                    replace(outCY$HR_BY[["IM"]][[i]][1:79,],outCY$terminal,0),2,sum)[1:outCY$MaxStkAge]
}

maturation_rate_long <- reshape2::melt(outCY$MatRte)|>
                        rename( "brood"=Var1, "age"=Var2, "maturation_rate"=value )|>
                        filter(age>1)
                        

mature_ocean_net_long<-reshape2::melt(mature_ocean_net)|>
                        rename( "brood"=Var1, "age"=Var2, "mature_ocean_net_hr"=value )|>
                        mutate(brood = brood+outCY$FirstBY-1)|>
                        filter(age>1)


preterminal_er_long <- reshape2::melt(preterminal_er)|>
                        rename( "brood"=Var1, "age"=Var2, "preterminal_er"=value )|>
                        mutate(brood = brood+outCY$FirstBY-1)|>
                        filter(age>1)


aeq_long <- reshape2::melt(outCY$AdltEqv)|>
                        rename( "brood"=Var1, "age"=Var2, "aeq"=value )|>
                        filter(age>1&age<6)

preterminal_er_long$natural_survival<-outCY$SurvRte[preterminal_er_long$age-1]


#process the escapement data

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


#ps. in the original spreadsheed, ER of 1 was replaced with an average of previous and current years, i.e. it was assumed that such high 
#preterminal ERs were not possible, which can be reasonable as there were escapement observations


all_escapement_data$preterminal_er_smooth<-all_escapement_data$preterminal_er
#Alternative implementation of ocean pre fishery with smoothing of Er rates
pretermerone<-which(1-round(all_escapement_data$preterminal_er,5)<=0)
for(i in pretermerone){
  b<-all_escapement_data$brood[i]
  a<-all_escapement_data$age[i]

   all_escapement_data$preterminal_er_smooth[i]<- mean(all_escapement_data$preterminal_er[
                                                                (all_escapement_data$brood==(b-1)|
                                                                all_escapement_data$brood==(b+1))&
                                                                all_escapement_data$age==a])

}


all_escapement_data$preterminal_er_smooth[pretermerone]
all_escapement_data$preterminal_er[pretermerone]


all_escapement_data$ocean_pre_fishery_abundance <- ifelse((1-round(all_escapement_data$preterminal_er,5))>0,
                                                   all_escapement_data$ocean_post_fishery_abundance/(1-all_escapement_data$preterminal_er),
                                                   0)




all_escapement_data$ocean_pre_fishery_abundance <- ifelse((1-round(all_escapement_data$preterminal_er,5))>0,
                                                   all_escapement_data$ocean_post_fishery_abundance/(1-all_escapement_data$preterminal_er),
                                                   0)


#(all_escapement_data$preterminal_er)[which(all_escapement_data$ocean_catch_aeq>1000000)]
#all_escapement_data[which(all_escapement_data$ocean_catch_aeq>1000000),]


all_escapement_data$cohort_abundance <- all_escapement_data$ocean_pre_fishery_abundance/
                                        (all_escapement_data$natural_survival)


all_escapement_data$ocean_mortality_immature<-all_escapement_data$ocean_pre_fishery_abundance-
                                           all_escapement_data$ocean_post_fishery_abundance


all_escapement_data$ocean_catch_aeq<-all_escapement_data$ocean_mortality_immature*all_escapement_data$aeq

all_escapement_data$adult_recruits<-all_escapement_data$ocean_catch_aeq+all_escapement_data$total_mature_run

all_escapement_data$maturity_rate_natural_population <- ifelse(all_escapement_data$ocean_post_fishery_abundance>0,
                                                           all_escapement_data$total_mature_run/all_escapement_data$ocean_post_fishery_abundance,
                                                           1)

#recruits 
recruits<-reshape2::dcast(all_escapement_data, brood ~ age, value.var = 'adult_recruits')
recruits$rec_total<-apply(recruits[,3:5],1,sum)
recruits$year<-recruits$brood

female_spawner<-reshape2::dcast(all_escapement_data, year ~ age, value.var = 'female_spawner',fun.aggregate =sum)
female_spawner$female_spawners_3_5<-apply(female_spawner[,3:5],1,sum)

all_spawners<-reshape2::dcast(river_spawners, year ~ age, value.var = 'river_spawners',fun.aggregate =sum)
all_spawners$all_spawners_3_5<-apply(all_spawners[,3:5],1,sum)

spawners_prop<-left_join(female_spawner[,c("year","female_spawners_3_5")],all_spawners[,c("year","all_spawners_3_5")])

spawners_prop$female_rate<-spawners_prop$female_spawners_3_5/spawners_prop$all_spawners_3_5

all_escapement_data$adjusted_escapement <- all_escapement_data$river_spawners/mean(spawners_prop$female_rate,na.rm=T)


#Compute the SR series
adjescapement<-reshape2::dcast(all_escapement_data, year ~ age, value.var = 'adjusted_escapement',
                          fun.aggregate =sum)

adjescapement$esc_total<- apply(adjescapement[,3:5],1,sum)

sr_data<-left_join(adjescapement|>select(year,esc_total),
          recruits|>select(year,rec_total))

library(ggplot2)
ggplot(sr_data)+
geom_point(aes(x=esc_total, y=rec_total, color=as.factor(year)))+
theme_bw()


left_join

escapemeyt
