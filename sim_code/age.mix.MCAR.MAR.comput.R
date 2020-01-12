# Master model for simulation and inferrence of age-mixing patterns in transmission clusters
# The output is a vector of values of

# 1. Epidemic and sexual behaviour statistics

# "R.prev.15.25.w", "R.prev.15.25.m", "R.prev.25.40.w", "R.prev.25.40.m", "R.prev.40.50.w", "R.prev.40.50.m",
# "R.AAD.male", "R.SDAD.male", "R.slope.male", "R.WSD.male", "R.BSD.male", "R.intercept.male"
# "R.p.prev.6months.m", # "R.p.prev.6months.f",
# "R.inc.15.25.m", "R.inc.15.25.w", "R.inc.25.40.m", "R.inc.25.40.w", "R.inc.40.50.m", "R.inc.40.50.w"

# 2. What is produced by age.mixing.MCAR.fun at 100% of coverage

# 3. With MCAR and MAR (35:95, by 5) scenarios (13 * 4), each scenario returns measurements which are describbed
# in age.mixing.MCAR.fun and age.mixing.MCAR.fun scripts):

# Missing Completly at Random has 13 scenarios
# Missing At Random has 39 scenarios, with 13 when we assume we have more women in the sample (seq.gender.ratio = 70%)
# the second we have fewer women (seq.gender.ratio = 30%), and the third we have same amount of men and women (seq.gender.ratio = 50%)


# Make sure you have seq-gen, FastTree, and comandline ClusterPicker_1.2.3 in you working directory

# Define directory

# work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop

# work.dir <- "/home/dniyukuri/lustre/agemix.25.10.2018.2" # on CHPC




age.mix.MCAR.MAR.comput <- function(inputvector=inputvector){
  
  
  
  
  # source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/age_mixing_large_AD/needed.functions.RSimpactHelp.R")
  # source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/age_mixing_large_AD/advanced.transmission.network.builder.R")
  # source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/age_mixing_large_AD/age.mixing.MCAR.fun.R")
  # source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/age_mixing_large_AD/age.mixing.MAR.fun.R")
  
  source("/home/dniyukuri/lustre/age_mixing_large_AD/needed.functions.RSimpactHelp.R")
  source("/home/dniyukuri/lustre/age_mixing_large_AD/advanced.transmission.network.builder.R")
  source("/home/dniyukuri/lustre/age_mixing_large_AD/age.mixing.MCAR.fun.R")
  source("/home/dniyukuri/lustre/age_mixing_large_AD/age.mixing.MAR.fun.R")
  
  # work.dir <- "/home/niyukuri/Dropbox/25.10.2018.age.mix2/age_mixing_large_AD/" # on PC
  
  work.dir <- "/home/dniyukuri/lustre/age_mixing_large_AD" # on PCHPC
  
  
  # work.dir <- "/home/david/age_mixing_AD_clusters" # on PCHPC
  
  
  setwd(paste0(work.dir))
  
  
  library(RSimpactCyan)
  library(RSimpactHelper)
  library(Rcpp)
  library(ape)
  library(expoTree)
  library(data.table)
  library(readr)
  library(phangorn)
  library(lme4)
  library(nlme)
  library(dplyr)
  library(adephylo)
  library(treedater)
  library(geiger)
  library(picante)
  library(igraph)
  library(phyloTop)
  library(phytools)
  library(Rsamtools)
  library(robustbase)
  library(intergraph)
  library(lubridate)
  library(tidyr)
  library(data.table)
  
  
  ###########################################
  # Step 1: Setup and running simpact      #
  ###########################################
  

  
  ## Run Simpact for specific parameter combination
  
  age.distr <- agedistr.creator(shape = 5, scale = 65)
  #
  cfg.list <- input.params.creator(population.eyecap.fraction = 0.2,
                                   population.simtime = 40, 
                                   population.nummen = 10000, 
                                   population.numwomen = 10000,
                                   hivseed.time = 10, 
                                   hivseed.type = "amount",
                                   hivseed.amount = 10, 
                                   hivseed.age.min = 20,
                                   hivseed.age.max = 50,
                                   formation.hazard.agegapry.meanage = -0.025,
                                   debut.debutage = 15
  )
  
  # # Assumption of nature of sexual network
  # #########################################
  #
  cfg.list["population.msm"] = "no"
  
  
  # # Sexual behaviour
  # ###################
  #
  seedid <- inputvector[1]
  
  cfg.list["dissolution.alpha_0"] <- inputvector[2] # [1] # -0.52 c("unif", -1, 0)
  cfg.list["dissolution.alpha_4"] <- inputvector [3] # [2] # -0.05 c("unif", -0.5, 0)
  cfg.list["formation.hazard.agegapry.baseline"] <- inputvector[4] # [3] # 2 c("unif", 1, 3)
  cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[5] # [4] # 0 c("unif", -0.5, 0.5)
  cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[5] # [4] # 0
  cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[6] # [5] # 3 c("unif", 2, 4)
  cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[6] # [5] # 3 
  cfg.list["formation.hazard.agegapry.gap_agescale_man"] <- inputvector[7] # [6] # 0.25 c("unif", 0, 1)
  cfg.list["formation.hazard.agegapry.gap_agescale_woman"] <- inputvector[7] # [6] # 0.25
  cfg.list["formation.hazard.agegapry.numrel_man"] <- inputvector[8] # [7] # -0.3 c("unif", -1, 0)
  cfg.list["formation.hazard.agegapry.numrel_woman"] <- inputvector[8] # [7] # -0.3
  cfg.list["formation.hazard.agegapry.numrel_diff"] <- inputvector[9] # [8] # -0.1 c("unif", -0.9, 0)
  
  
  # # HIV transmission
  # ###################
  #
  
  cfg.list["hivtransmission.param.a"] <- inputvector[10] # [10] # -1 c("unif", -2, 0)
  cfg.list["hivtransmission.param.b"] <- inputvector[11] # [11] # -90 c("unif", -100, -80)
  cfg.list["hivtransmission.param.c"] <- inputvector[12] # [12] # 0.5 c("unif", 0, 1)
  cfg.list["hivtransmission.param.f1"] <- inputvector[13] # [13] # 0.04879016 c("unif", 0, 0.5)
  cfg.list["hivtransmission.param.f2"] <- inputvector[14] # [14] # -0.1386294 c("unif", -0.5, 0)
  
  # Disease progression > may be remove in parameter to estimates
  
  cfg.list["person.vsp.toacute.x"] <- inputvector[15] # [15] # 5 c("unif", 3, 7)
  cfg.list["person.vsp.toaids.x"] <- inputvector[16] # [16] # 7 c("unif", 5, 9)
  cfg.list["person.vsp.tofinalaids.x"] <- inputvector[17] # [17] # 12 c("unif", 10, 14)
  
  
  #
  # # Demographic
  # ##############
  #
  
  cfg.list["conception.alpha_base"] <- inputvector[18] # [18] # -2.7 c("unif", -3.5, -1.7)
  
  
  # # Assumptions to avoid negative branch lengths
  # ###############################################
  # # + sampling == start ART
  # # when someone start ART, he/she is sampled and becomes non-infectious
  
  cfg.list["monitoring.fraction.log_viralload"] <- 0
  
  
  #
  # ## Add-ons
  #
  ### BEGIN Add-on
  cfg.list["mortality.aids.survtime.C"] <- 65
  cfg.list["mortality.aids.survtime.k"] <- -0.2
  cfg.list["monitoring.fraction.log_viralload"] <- 0 #0.3
  cfg.list["dropout.interval.dist.type"] <- "uniform"
  cfg.list["dropout.interval.dist.uniform.min"] <- 1000
  cfg.list["dropout.interval.dist.uniform.max"] <- 2000
  
  cfg.list["person.survtime.logoffset.dist.type"] <- "normal"
  cfg.list["person.survtime.logoffset.dist.normal.mu"] <- 0
  cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1
  
  cfg.list["person.agegap.man.dist.type"] <- "normal" #fixed
  #cfg.list["person.agegap.man.dist.fixed.value"] <- -6
  cfg.list["person.agegap.woman.dist.type"] <- "normal" #"fixed"
  #cfg.list["person.agegap.woman.dist.fixed.value"] <- -6
  
  cfg.list["mortality.aids.survtime.C"] <- 65
  cfg.list["mortality.aids.survtime.k"] <- -0.2
  cfg.list["monitoring.cd4.threshold"] <- 0 # 0 means nobody qualifies for ART
  cfg.list["diagnosis.baseline"] <- -2
  
  
  cfg.list["person.eagerness.man.dist.gamma.a"] <- 0.23 # 0.23
  cfg.list["person.eagerness.woman.dist.gamma.a"] <- 0.23 # 0.23
  cfg.list["person.eagerness.man.dist.gamma.b"] <- 45 # 45
  cfg.list["person.eagerness.woman.dist.gamma.b"] <- 45 # 45
  
  #### END Add-ons
  
  
  # # ART intervention
  # ###################
  #
  # # ART acceptability paramter and the ART  interventions
  
  cfg.list["person.art.accept.threshold.dist.fixed.value"] <- 0.6
  
  # Let's introduce ART, and evaluate whether the HIV prevalence drops less  rapidly
  art.intro <- list()
  art.intro["time"] <- 20
  art.intro["diagnosis.baseline"] <- -2 # 0#100
  art.intro["monitoring.cd4.threshold"] <- 100 # 1200
  
  ### add something about diagnosis
  art.intro["diagnosis.agefactor"] <- 0
  art.intro["diagnosis.genderfactor"] <- 0
  art.intro["diagnosis.diagpartnersfactor"] <- 0
  art.intro["diagnosis.isdiagnosedfactor"] <- 0
  ### end of add-on about diagnosis
  #art.intro["monitoring.interval.piecewise.cd4s"] <- "0,1300"
  # Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2013:500
  art.intro1 <- list()
  art.intro1["time"] <- 22
  art.intro1["diagnosis.baseline"] <- -2 # 0#100
  art.intro1["monitoring.cd4.threshold"] <- 150 # 1200
  
  art.intro2 <- list()
  art.intro2["time"] <- 25 # inputvector[5] ######### 30
  art.intro2["monitoring.cd4.threshold"] <- 200
  
  art.intro3 <- list()
  art.intro3["time"] <- 30 # inputvector[4] + inputvector[5] + inputvector[6] ########### 33
  art.intro3["monitoring.cd4.threshold"] <- 350
  
  art.intro4 <- list()
  art.intro4["time"] <- 33 # inputvector[4] + inputvector[5] + inputvector[6] + inputvector[7] ########### 36
  art.intro4["monitoring.cd4.threshold"] <- 500
  
  art.intro5 <- list()
  art.intro5["time"] <- 36
  art.intro5["monitoring.cd4.threshold"] <- 700 # This is equivalent to immediate access
  
  # tasp.indicator <- inputvector[9] # 1 if the scenario is TasP, 0 if the scenario is current status
  interventionlist <- list(art.intro, art.intro1, art.intro2, art.intro3, art.intro4, art.intro5)
  
  intervention <- interventionlist
  
  # Events
  cfg.list["population.maxevents"] <- as.numeric(cfg.list["population.simtime"][1]) * as.numeric(cfg.list["population.nummen"][1]) * 3
  
  # Avoid overlaping in same directory
  
  #creating subfolder with unique name for each simulation
  generate.filename <- function(how.long){
    
    rn <- sample(1:100,1)
    t <- as.numeric(Sys.time())
    set.seed((t - floor(t)) * 1e8)
    chars <- c(letters, LETTERS)
    sub.dir.sim.id <-  paste0(sample(chars,how.long), collapse = "")
    
    noise.sample1 <- sample(8:15,1, replace = TRUE)
    sub.dir.sim.id.ext <- paste0(sample(chars,noise.sample1), collapse = "")
    noise.sample <- sample(1:1000,1)
    noise.sample2 <- sample(8:17,1, replace = TRUE)
    sub.dir.sim.id <- paste0(sub.dir.sim.id.ext,
                             paste0(sample(chars,noise.sample2), collapse = ""),noise.sample, rn)
    
    return(sub.dir.sim.id)
  }
  
  
  
  sub.dir.rename <- paste0(work.dir,"/temp/",generate.filename(10))
  
  
  
  
  # Running Simpact 
  #################
  
  results <- tryCatch(simpact.run(configParams = cfg.list,
                                  destDir = sub.dir.rename,
                                  agedist = age.distr,
                                  seed = seedid,
                                  intervention = intervention),
                      error = simpact.errFunction)
  
  
  if (length(results) == 0){
    results.mcar <- rep(NA, 5908) # 37 + 82 + 33 + 1 + 1 + 53 = 207
  }else{
    if (as.numeric(results["eventsexecuted"]) >= (as.numeric(cfg.list["population.maxevents"]) - 1)){
      results.mcar <- rep(NA, 5908)
    }else{
      
      
      DataListALL <- readthedata(results)
      
      
      datalist.agemix <- DataListALL
      
      
      
      # datalist.agemix <- get(load("datalist.agemix.RData"))
      
      
      
      
      ###########################################
      # Step 2: Construct transmission networks #
      ###########################################
      
      
      simpact.trans.net.adv <- advanced.transmission.network.builder(datalist = datalist.agemix, endpoint = 40)
      
      
      net.size.vector <- vector() # i_th seed in the list of seeds
      
      for(i in 1:length(simpact.trans.net.adv)){
        
        tree.n <- simpact.trans.net.adv[[i]] # transmission network for i^th seed
        
        net.size.vector <- c(net.size.vector, nrow(as.data.frame(tree.n)))
        
      }
      
      big.index <- which(net.size.vector>=100)
      
      if(length(big.index) >= 1){ 
        
        # simpact.trans.net <- transmission.network.builder(datalist = datalist.agemix, endpoint = 40)
        
        
        
        # simpact.trans.net.projection <- transmission.network.builder(datalist = datalist.agemix, endpoint = 45)
        
        
        
        ###############################
        # Step 3: Sequence simulation #
        ###############################
        
        
        trans.net <- simpact.trans.net.adv # simpact.trans.net # all transmission networks
        
        
        dirseqgen <- work.dir
        
        seeds.num <- inputvector[1]
        
        # Sequence simulation is done for at least a transmission network with 6 individuals
        # This means that limitTransmEvents equal at least 7
        
        sequence.simulation.seqgen.par(dir.seq = dirseqgen,
                                       sub.dir.rename = sub.dir.rename,
                                       simpact.trans.net = simpact.trans.net.adv, # simpact.trans.net,
                                       seq.gen.tool = "seq-gen",
                                       seeds.num = seeds.num,
                                       endpoint = 40,
                                       limitTransmEvents = 7, # no less than 7
                                       hiv.seq.file = "hiv.seq.C.pol.j.fasta",
                                       clust = TRUE) # hiv.seq.file lodged in work.dir
        
        # Transform the sequence format to be handled by ClusterPicker
        sequ.dna <- read.dna(file = paste0(sub.dir.rename,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta"), format = "interleaved")
        write.dna(sequ.dna, file = paste0(sub.dir.rename,"/C.Epidemic.fas") , format = "fasta")
        
        
        
        ###################################################################
        # Step 4: Epidemic statistics and sexual behaviour: full data set #
        ###################################################################
        
        
        
        # (i) Age mixing in relationships
        ##################################
        
        # 
        
        agemix.rels.df <- agemix.df.maker(datalist.agemix)
        
        # 
        agemix.model <- pattern.modeller(dataframe = agemix.rels.df,
                                         agegroup = c(15, 50),
                                         timepoint = 40, # datalist.agemix$itable$population.simtime[1],
                                         timewindow = 5)#1)#3)
        # 
        # # men.lme <- tryCatch(agemixing.lme.fitter(data = dplyr::filter(agemix.model[[1]], Gender =="male")),
        # #                     error = agemixing.lme.errFunction) # Returns an empty list if the lme model can't be fitted
        #
        # men.lmer <- ampmodel(data = dplyr::filter(agemix.model[[1]], Gender =="male"))
        
        data = dplyr::filter(agemix.model[[1]], Gender =="male")
        
        if( nrow(data) > length(unique(data$ID)) & length(unique(data$ID)) > 1 ){
          
          men.lmer <- lmer(pagerelform ~ agerelform0 + (1 | ID),
                           data = dplyr::filter(agemix.model[[1]], Gender =="male"),
                           REML = TRUE,
                           control=lmerControl(check.nobs.vs.nlev = "ignore",
                                               check.nobs.vs.rankZ = "ignore",
                                               check.nobs.vs.nRE="ignore"))
          
          bignumber <- NA # let's try if NA works (instead of 9999 for example)
          AAD.male <- ifelse(length(men.lmer) > 0, mean(dplyr::filter(agemix.model[[1]], Gender =="male")$AgeGap), bignumber)
          SDAD.male <- ifelse(length(men.lmer) > 0, sd(dplyr::filter(agemix.model[[1]], Gender =="male")$AgeGap), bignumber)
          #powerm <- ifelse(length(men.lme) > 0, as.numeric(attributes(men.lme$apVar)$Pars["varStruct.power"]), bignumber)
          slope.male <- ifelse(length(men.lmer) > 0, summary(men.lmer)$coefficients[2, 1], bignumber) #summary(men.lmer)$tTable[2, 1], bignumber)
          WSD.male <- ifelse(length(men.lmer) > 0, summary(men.lmer)$sigma, bignumber) #WVAD.base <- ifelse(length(men.lme) > 0, men.lme$sigma^2, bignumber)
          
          BSD.male <- ifelse(length(men.lmer) > 0, bvar(men.lmer), bignumber) # Bad name for the function because it actually extracts between subject standard deviation # BVAD <- ifelse(length(men.lmer) > 0, getVarCov(men.lme)[1,1], bignumber)
          
          intercept.male <- ifelse(length(men.lmer) > 0, summary(men.lmer)$coefficients[1,1] - 15, bignumber)
          
          # c(AAD.male, SDAD.male, slope.male, WSD.male, BSD.male, intercept.male)
          
          ## AAD: average age difference across all relationship
          ## VAD: variance of these age differences
          ## SDAD: standard deviation of age differences
          ## BSD: between-subject standard deviation of age differences
          
          mix.rels.dat <- c(AAD.male, SDAD.male, slope.male, WSD.male, BSD.male, intercept.male)
          
          names(mix.rels.dat) <-  c("R.AAD.male", "R.SDAD.male", "R.slope.male", "R.WSD.male", "R.BSD.male", "R.intercept.male")
          
        }else{
          
          mix.rels.dat <- rep(NA, 6)
          
          names(mix.rels.dat) <-  c("R.AAD.male", "R.SDAD.male", "R.slope.male", "R.WSD.male", "R.BSD.male", "R.intercept.male")
          
        }
        
        # age.scatter.df <- agemix.model[[1]]
        
        #  (ii) Point 	prevalence of concurrency in the adult population
        ##################################################################
        
        
        # Concurrency point prevalence 6 months before a survey, among men
        
        pp.cp.6months.male.rels <- tryCatch(concurr.pointprev.calculator(datalist = datalist.agemix,
                                                                         timepoint = 40 - 0.5),
                                            error=function(e) return(rep(NA, 1)))
        
        
        # pp.cp.6months.male.rels <- tryCatch(concurr.pointprev.calculator(datalist = datalist.agemix,
        #                                                         timepoint = 40 - 0.5) %>%
        #   dplyr::select(concurr.pointprev) %>%
        #   dplyr::slice(1) %>%
        #   as.numeric(),
        #   error=function(e) return(rep(NA, 1)))
        
        # 
        # pp.cp.6months.female.rels <- concurr.pointprev.calculator(datalist = datalist.agemix,
        #                                                           timepoint = 40 - 0.5) %>%
        #   dplyr::select(concurr.pointprev) %>%
        #   dplyr::slice(2) %>%
        #   as.numeric()
        # 
        # 
        
        
        # (iii) Prevalence
        ##################
        
        
        # hiv.prev.lt25.women <- prevalence.calculator(datalist = datalist.agemix,
        #                                              agegroup = c(15, 25),
        #                                              timepoint = 40) %>%
        #   dplyr::select(pointprevalence) %>%
        #   dplyr::slice(2) %>%
        #   as.numeric()
        # 
        # hiv.prev.lt25.men <- prevalence.calculator(datalist = datalist.agemix,
        #                                            agegroup = c(15, 25),
        #                                            timepoint = 40) %>%
        #   dplyr::select(pointprevalence) %>%
        #   dplyr::slice(1) %>%
        #   as.numeric()
        # 
        # hiv.prev.25.40.women <- prevalence.calculator(datalist = datalist.agemix,
        #                                               agegroup = c(25, 40),
        #                                               timepoint = 40) %>%
        #   dplyr::select(pointprevalence) %>%
        #   dplyr::slice(2) %>%
        #   as.numeric()
        # 
        # hiv.prev.25.40.men <- prevalence.calculator(datalist = datalist.agemix,
        #                                             agegroup = c(25, 40),
        #                                             timepoint = 40) %>%
        #   dplyr::select(pointprevalence) %>%
        #   dplyr::slice(1) %>%
        #   as.numeric()
        # 
        # hiv.prev.40.50.women <- prevalence.calculator(datalist = datalist.agemix,
        #                                               agegroup = c(40, 50),
        #                                               timepoint = 40) %>%
        #   dplyr::select(pointprevalence) %>%
        #   dplyr::slice(2) %>%
        #   as.numeric()
        # 
        # hiv.prev.40.50.men <- prevalence.calculator(datalist = datalist.agemix,
        #                                             agegroup = c(40, 50),
        #                                             timepoint = 40) %>%
        #   dplyr::select(pointprevalence) %>%
        #   dplyr::slice(1) %>%
        #   as.numeric()
        
        # Women
        
        hiv.prev.lt25.women <- prevalence.calculator(datalist = datalist.agemix,
                                                     agegroup = c(15, 25),
                                                     timepoint = 40) %>%
          dplyr::select(pointprevalence) %>%
          dplyr::slice(2) %>%
          as.numeric()
        
        hiv.prev.25.30.women <- prevalence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(25, 30),
                                                      timepoint = 40) %>%
          dplyr::select(pointprevalence) %>%
          dplyr::slice(2) %>%
          as.numeric()
        
        hiv.prev.30.35.women <- prevalence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(30, 35),
                                                      timepoint = 40) %>%
          dplyr::select(pointprevalence) %>%
          dplyr::slice(2) %>%
          as.numeric()
        
        hiv.prev.35.40.women <- prevalence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(35, 40),
                                                      timepoint = 40) %>%
          dplyr::select(pointprevalence) %>%
          dplyr::slice(2) %>%
          as.numeric()
        
        hiv.prev.40.45.women <- prevalence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(40, 45),
                                                      timepoint = 40) %>%
          dplyr::select(pointprevalence) %>%
          dplyr::slice(2) %>%
          as.numeric()
        
        hiv.prev.45.50.women <- prevalence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(45, 50),
                                                      timepoint = 40) %>%
          dplyr::select(pointprevalence) %>%
          dplyr::slice(2) %>%
          as.numeric()
        
        
        
        # Men
        
        hiv.prev.lt25.men <- prevalence.calculator(datalist = datalist.agemix,
                                                   agegroup = c(15, 25),
                                                   timepoint = 40) %>%
          dplyr::select(pointprevalence) %>%
          dplyr::slice(1) %>%
          as.numeric()
        
        hiv.prev.25.30.men <- prevalence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(25, 30),
                                                    timepoint = 40) %>%
          dplyr::select(pointprevalence) %>%
          dplyr::slice(1) %>%
          as.numeric()
        
        hiv.prev.30.35.men <- prevalence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(30, 35),
                                                    timepoint = 40) %>%
          dplyr::select(pointprevalence) %>%
          dplyr::slice(1) %>%
          as.numeric()
        
        hiv.prev.35.40.men <- prevalence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(35, 40),
                                                    timepoint = 40) %>%
          dplyr::select(pointprevalence) %>%
          dplyr::slice(1) %>%
          as.numeric()
        
        hiv.prev.40.45.men <- prevalence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(40, 45),
                                                    timepoint = 40) %>%
          dplyr::select(pointprevalence) %>%
          dplyr::slice(1) %>%
          as.numeric()
        
        hiv.prev.45.50.men <- prevalence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(45, 50),
                                                    timepoint = 40) %>%
          dplyr::select(pointprevalence) %>%
          dplyr::slice(1) %>%
          as.numeric()
        
        
        
        
        
        # agegroup <- c("15-24", "25-29", "30-34", "35-39", "40-44", "45-49")
        
        men.prev <- c(hiv.prev.lt25.men, hiv.prev.25.30.men, hiv.prev.30.35.men, hiv.prev.35.40.men, hiv.prev.40.45.men, hiv.prev.45.50.men)
        names(men.prev) <- c("prev.m.15.24", "prev.m.25.29", "prev.m.30.34", "prev.m.35.39", "prev.m.40.44", "prev.m.45.49")
        women.prev <- c(hiv.prev.lt25.women, hiv.prev.25.30.women, hiv.prev.30.35.women, hiv.prev.35.40.women, hiv.prev.40.45.women, hiv.prev.45.50.women)
        names(women.prev) <- c("prev.w.15.24", "prev.w.25.29", "prev.w.30.34", "prev.w.35.39", "prev.w.40.44", "prev.w.45.49")        
        
        
        
        
        # (iv) Incidence
        #################
        
        
        
        epi.rels.incidence.df.15.24.men <- incidence.calculator(datalist = datalist.agemix,
                                                                agegroup = c(15, 25),
                                                                timewindow = c(35, 40),
                                                                only.active = "No") %>%
          dplyr::select(incidence) %>%
          dplyr::slice(1) %>%
          as.numeric()
        
        epi.rels.incidence.df.15.24.women <- incidence.calculator(datalist = datalist.agemix,
                                                                  agegroup = c(15, 25),
                                                                  timewindow = c(35, 40),
                                                                  only.active = "No") %>%
          dplyr::select(incidence) %>%
          dplyr::slice(2) %>%
          as.numeric()
        
        
        
        epi.rels.incidence.df.25.39.men <- incidence.calculator(datalist = datalist.agemix,
                                                                agegroup = c(25, 40),
                                                                timewindow = c(35, 40),
                                                                only.active = "No") %>%
          dplyr::select(incidence) %>%
          dplyr::slice(1) %>%
          as.numeric()
        
        epi.rels.incidence.df.25.39.women <- incidence.calculator(datalist = datalist.agemix,
                                                                  agegroup = c(25, 40),
                                                                  timewindow = c(35, 40),
                                                                  only.active = "No") %>%
          dplyr::select(incidence) %>%
          dplyr::slice(2) %>%
          as.numeric()
        
        
        
        epi.rels.incidence.df.40.49.men <- incidence.calculator(datalist = datalist.agemix,
                                                                agegroup = c(40, 50),
                                                                timewindow = c(35, 40),
                                                                only.active = "No") %>%
          dplyr::select(incidence) %>%
          dplyr::slice(1) %>%
          as.numeric()
        
        epi.rels.incidence.df.40.49.women <- incidence.calculator(datalist = datalist.agemix,
                                                                  agegroup = c(40, 50),
                                                                  timewindow = c(35, 40),
                                                                  only.active = "No") %>%
          dplyr::select(incidence) %>%
          dplyr::slice(2) %>%
          as.numeric()
        
        
        
        summary.epidemic.rels.df <- c(men.prev, women.prev,
                                      epi.rels.incidence.df.15.24.men, epi.rels.incidence.df.15.24.women, 
                                      epi.rels.incidence.df.25.39.men, epi.rels.incidence.df.25.39.women,
                                      epi.rels.incidence.df.40.49.men, epi.rels.incidence.df.40.49.women,
                                      pp.cp.6months.male.rels, # , # pp.cp.6months.female.rels, 
                                      mix.rels.dat) 
        
        names(summary.epidemic.rels.df) <- c(names(men.prev), names(women.prev),
                                             "R.inc.15.25.m", "R.inc.15.25.w", "R.inc.25.40.m", 
                                             "R.inc.25.40.w", "R.inc.40.50.m", "R.inc.40.50.w",
                                             "R.p.prev.6months.m",
                                             names(mix.rels.dat)) 
        
        
        
        # 
        # # True age structure in transmission transmission network for selected individuals #
        # #####################################################################################
        # 
        # age.groups.filtered.transmission.net.fun <- function(table.transmission.net.cov = table.simpact.trans.net.cov,
        #                                                      age.group.15.25 = c(15,25),
        #                                                      age.group.25.40 = c(25,40),
        #                                                      age.group.40.50 = c(40,50)){
        #   
        #   num.women.15.25 <- dplyr::filter(table.transmission.net.cov, 
        #                                    table.transmission.net.cov$GenderRec=="1" & table.transmission.net.cov$ageSampTimeRec >= age.group.15.25[1] & table.transmission.net.cov$ageSampTimeRec < age.group.15.25[2])
        #   
        #   num.men.15.25 <- dplyr::filter(table.transmission.net.cov, 
        #                                  table.transmission.net.cov$GenderRec=="0" & table.transmission.net.cov$ageSampTimeRec >= age.group.15.25[1] & table.transmission.net.cov$ageSampTimeRec < age.group.15.25[2])
        #   
        #   
        #   num.women.25.40 <- dplyr::filter(table.transmission.net.cov, 
        #                                    table.transmission.net.cov$GenderRec=="1" & table.transmission.net.cov$ageSampTimeRec >= age.group.25.40[1] & table.transmission.net.cov$ageSampTimeRec < age.group.25.40[2])
        #   
        #   
        #   num.men.25.40 <- dplyr::filter(table.transmission.net.cov, 
        #                                  table.transmission.net.cov$GenderRec=="0" & table.transmission.net.cov$ageSampTimeRec >= age.group.25.40[1] & table.transmission.net.cov$ageSampTimeRec < age.group.25.40[2])
        #   
        #   
        #   
        #   num.women.40.50 <- dplyr::filter(table.transmission.net.cov, 
        #                                    table.transmission.net.cov$GenderRec=="1" & table.transmission.net.cov$ageSampTimeRec >= age.group.40.50[1] & table.transmission.net.cov$ageSampTimeRec < age.group.40.50[2])
        #   
        #   
        #   
        #   num.men.40.50 <- dplyr::filter(table.transmission.net.cov, 
        #                                  table.transmission.net.cov$GenderRec=="0" & table.transmission.net.cov$ageSampTimeRec >= age.group.40.50[1] & table.transmission.net.cov$ageSampTimeRec < age.group.40.50[2])
        #   
        #   
        #   numbers.indiv.women.15.25 <- nrow(num.women.15.25)
        #   numbers.indiv.men.15.25 <- nrow(num.men.15.25)
        #   numbers.indiv.women.25.40 <- nrow(num.women.25.40)
        #   numbers.indiv.men.25.40 <- nrow(num.men.25.40)
        #   numbers.indiv.women.40.50 <- nrow(num.women.40.50)
        #   numbers.indiv.men.40.50 <- nrow(num.men.40.50)
        #   
        #   numbers.individuals.age.groups <- c(numbers.indiv.women.15.25, numbers.indiv.men.15.25, 
        #                                       numbers.indiv.women.25.40, numbers.indiv.men.25.40,
        #                                       numbers.indiv.women.40.50, numbers.indiv.men.40.50)
        #   
        #   names(numbers.individuals.age.groups) <- c("num.women.true.net.15.25", "num.men.true.net.15.25", 
        #                                              "num.women.true.net.25.40", "num.men.true.net.25.40",
        #                                              "num.women.true.net.40.50", "num.men.true.net.40.50")
        #   
        #   
        #   
        #   num.women.15.25$ageSampTimeDon <- num.women.15.25$SampTime - num.women.15.25$TOBDon
        #   num.men.15.25$ageSampTimeDon <- num.men.15.25$SampTime - num.men.15.25$TOBDon
        #   
        #   num.women.25.40$ageSampTimeDon <- num.women.25.40$SampTime - num.women.25.40$TOBDon
        #   num.men.25.40$ageSampTimeDon <- num.men.25.40$SampTime - num.men.25.40$TOBDon
        #   
        #   num.women.40.50$ageSampTimeDon <- num.women.40.50$SampTime - num.women.40.50$TOBDon
        #   num.men.40.50$ageSampTimeDon <- num.men.40.50$SampTime - num.men.40.50$TOBDon
        #   
        #   
        #   # Age differences
        #   
        #   AD.num.women.15.25 <- abs(num.women.15.25$ageSampTimeDon - num.women.15.25$ageSampTimeRec)
        #   AD.num.men.15.25 <- abs(num.men.15.25$ageSampTimeDon - num.men.15.25$ageSampTimeRec)
        #   
        #   AD.num.women.25.40 <- abs(num.women.25.40$ageSampTimeDon - num.women.25.40$ageSampTimeRec)
        #   AD.num.men.25.40 <- abs(num.men.25.40$ageSampTimeDon - num.men.25.40$ageSampTimeRec)
        #   
        #   AD.num.women.40.50 <- abs(num.women.40.50$ageSampTimeDon - num.women.40.50$ageSampTimeRec)
        #   AD.num.men.40.50 <- abs(num.men.40.50$ageSampTimeDon - num.men.40.50$ageSampTimeRec)
        #   
        #   
        #   mean.AD.num.women.15.25 <- mean(AD.num.women.15.25)
        #   med.AD.num.women.15.25 <- median(AD.num.women.15.25)
        #   sd.AD.num.women.15.25 <- sd(AD.num.women.15.25)
        #   
        #   mean.AD.num.men.15.25 <- mean(AD.num.men.15.25)
        #   med.AD.num.men.15.25 <- median(AD.num.men.15.25)
        #   sd.AD.num.men.15.25 <- sd(AD.num.men.15.25)
        #   
        #   
        #   mean.AD.num.women.25.40 <- mean(AD.num.women.25.40)
        #   med.AD.num.women.25.40 <- median(AD.num.women.25.40)
        #   sd.AD.num.women.25.40 <- sd(AD.num.women.25.40)
        #   
        #   mean.AD.num.men.25.40 <- mean(AD.num.men.25.40)
        #   med.AD.num.men.25.40 <- median(AD.num.men.25.40)
        #   sd.AD.num.men.25.40 <- sd(AD.num.men.25.40)
        #   
        #   
        #   mean.AD.num.women.40.50 <- mean(AD.num.women.40.50)
        #   med.AD.num.women.40.50 <- median(AD.num.women.40.50)
        #   sd.AD.num.women.40.50 <- sd(AD.num.women.40.50)
        #   
        #   mean.AD.num.men.40.50 <- mean(AD.num.men.40.50)
        #   med.AD.num.men.40.50 <- median(AD.num.men.40.50)
        #   sd.AD.num.men.40.50 <- sd(AD.num.men.40.50)
        #   
        #   
        #   mean.AD.age.groups <- c(mean.AD.num.women.15.25, mean.AD.num.men.15.25,
        #                           mean.AD.num.women.25.40, mean.AD.num.men.25.40,
        #                           mean.AD.num.women.40.50, mean.AD.num.men.40.50)
        #   
        #   names(mean.AD.age.groups) <- c("mean.AD.num.women.true.net.15.25", "mean.AD.num.men.true.net.15.25",
        #                                  "mean.AD.num.women.true.net.25.40", "mean.AD.num.men.true.net.25.40",
        #                                  "mean.AD.num.women.true.net.40.50", "mean.AD.num.men.true.net.40.50")
        #   
        #   
        #   
        #   med.AD.age.groups <- c(med.AD.num.women.15.25, med.AD.num.men.15.25,
        #                          med.AD.num.women.25.40, med.AD.num.men.25.40,
        #                          med.AD.num.women.40.50, med.AD.num.men.40.50)
        #   
        #   names(med.AD.age.groups) <- c("med.AD.num.women.true.net.15.25", "med.AD.num.men.true.net.15.25",
        #                                 "med.AD.num.women.true.net.25.40", "med.AD.num.men.true.net.25.40",
        #                                 "med.AD.num.women.true.net.40.50", "med.AD.num.men.true.net.40.50")
        #   
        #   
        #   
        #   sd.AD.age.groups <- c(sd.AD.num.women.15.25, sd.AD.num.men.15.25,
        #                         sd.AD.num.women.25.40, sd.AD.num.men.25.40,
        #                         sd.AD.num.women.40.50, sd.AD.num.men.40.50)
        #   
        #   names(sd.AD.age.groups) <- c("sd.AD.num.women.true.net.15.25", "sd.AD.num.men.true.net.15.25",
        #                                "sd.AD.num.women.true.net.25.40", "sd.AD.num.men.true.net.25.40",
        #                                "sd.AD.num.women.true.net.40.50", "sd.AD.num.men.true.net.40.50")
        #   
        #   
        #   table.transmission.net.cov$ageSampTimeDon <- table.transmission.net.cov$SampTime - table.transmission.net.cov$TOBDon
        #   
        #   men.df <- dplyr::filter(table.transmission.net.cov, table.transmission.net.cov$GenderDon=="0")
        #   women.df <- dplyr::filter(table.transmission.net.cov, table.transmission.net.cov$GenderDon=="1")
        #   
        #   Age.groups.table <- NULL
        #   
        #   filter.dat <- function(table.dat.fr = table.dat.fr){
        #     v1.dat <- vector()
        #     v2.dat <- vector()
        #     age1.dat <- vector()
        #     age2.dat <- vector()
        #     gender1.dat <- vector()
        #     gender2.dat <- vector()
        #     
        #     for(i in 1:nrow(table.dat.fr)){
        #       
        #       v1 <- table.dat.fr$RecId[i]
        #       v2 <- table.dat.fr$DonId[i]
        #       
        #       index.v1 <- which(table.dat.fr$RecId == v1)
        #       
        #       age1 <- table.dat.fr$ageSampTimeRec[index.v1]
        #       age2 <- table.dat.fr$ageSampTimeDon[index.v1]
        #       
        #       gender1 <- table.dat.fr$GenderRec[index.v1]
        #       gender2 <- table.dat.fr$GenderDon[index.v1]
        #       
        #       v1.dat <- c(v1.dat, v1)
        #       v2.dat <- c(v2.dat, v2)
        #       age1.dat <- c(age1.dat, age1)
        #       age2.dat <- c(age2.dat, age2)
        #       gender1.dat <- c(gender1.dat, gender1)
        #       gender2.dat <- c(gender2.dat, gender2)
        #       
        #     }
        #     
        #     age.table <- data.frame(v1.dat, gender1.dat, age1.dat, v2.dat, gender2.dat, age2.dat)
        #     
        #     return(age.table)
        #     
        #   }
        #   
        #   age.table <-   filter.dat(table.dat.fr = table.transmission.net.cov)
        #   
        #   # men as donors
        #   men.age.table.1 <-  filter.dat(table.dat.fr = men.df) # dplyr::filter(age.table, age.table$gender1.dat==0)
        #   
        #   # women as donors
        #   women.age.table.1 <-   filter.dat(table.dat.fr = women.df) # dplyr::filter(age.table, age.table$gender1.dat==1)
        #   
        #   
        #   # men 15.25 and women
        #   
        #   men.15.25.women.15.25.1 <- vector()
        #   men.15.25.women.25.40.1 <- vector()
        #   men.15.25.women.40.50.1 <- vector()
        #   
        #   if(nrow(men.age.table.1) >1 ){
        #     
        #     for (j in 1:nrow(men.age.table.1)) {
        #       
        #       
        #       if(men.age.table.1$age1.dat[j] >= age.group.15.25[1] & men.age.table.1$age1.dat[j] < age.group.15.25[2]){
        #         
        #         if(men.age.table.1$age2.dat[j] >= age.group.15.25[1] & men.age.table.1$age2.dat[j] < age.group.15.25[2]){
        #           
        #           men.15.25.women.15.25.1 <- c(men.15.25.women.15.25.1, men.age.table.1$age2.dat[j])
        #           
        #         }else if(men.age.table.1$age2.dat[j] >= age.group.25.40[1] & men.age.table.1$age2.dat[j] < age.group.25.40[2]){
        #           
        #           men.15.25.women.25.40.1 <- c(men.15.25.women.25.40.1, men.age.table.1$age2.dat[j])
        #           
        #         }else if (men.age.table.1$age2.dat[j] >= age.group.40.50[1] & men.age.table.1$age2.dat[j] < age.group.40.50[2]){
        #           
        #           men.15.25.women.40.50.1 <- c(men.15.25.women.40.50.1, men.age.table.1$age2.dat[j])
        #         }
        #         
        #       }
        #       
        #       
        #     }
        #     
        #   }
        #   
        #   
        #   # women 15.25 and men
        #   
        #   women.15.25.men.15.25.2 <- vector()
        #   women.15.25.men.25.40.2 <- vector()
        #   women.15.25.men.40.50.2 <- vector()
        #   
        #   if(nrow(women.age.table.1) >1 ){
        #     
        #     for (j in 1:nrow(women.age.table.1)) {
        #       
        #       
        #       if(women.age.table.1$age1.dat[j] >= age.group.15.25[1] & women.age.table.1$age1.dat[j] < age.group.15.25[2]){
        #         
        #         if(women.age.table.1$age2.dat[j] >= age.group.15.25[1] & women.age.table.1$age2.dat[j] < age.group.15.25[2]){
        #           
        #           women.15.25.men.15.25.2 <- c(women.15.25.men.15.25.2, women.age.table.1$age2.dat[j])
        #           
        #         }else if(women.age.table.1$age2.dat[j] >= age.group.25.40[1] & women.age.table.1$age2.dat[j] < age.group.25.40[2]){
        #           
        #           women.15.25.men.25.40.2 <- c(women.15.25.men.25.40.2, women.age.table.1$age2.dat[j])
        #           
        #         }else if (women.age.table.1$age2.dat[j] >= age.group.40.50[1] & women.age.table.1$age2.dat[j] < age.group.40.50[2]){
        #           
        #           women.15.25.men.40.50.2 <- c(women.15.25.men.40.50.2, women.age.table.1$age2.dat[j])
        #         }
        #         
        #       }
        #       
        #       
        #     }
        #     
        #   }
        #   
        #   
        #   
        #   # men 25.40 and women
        #   
        #   men.25.40.women.15.25.1 <- vector()
        #   men.25.40.women.25.40.1 <- vector()
        #   men.25.40.women.40.50.1 <- vector()
        #   
        #   if(nrow(men.age.table.1) > 1 ){
        #     
        #     for (j in 1:nrow(men.age.table.1)) {
        #       
        #       
        #       if(men.age.table.1$age1.dat[j] >= age.group.25.40[1] & men.age.table.1$age1.dat[j] < age.group.25.40[2]){
        #         
        #         if(men.age.table.1$age2.dat[j] >= age.group.15.25[1] & men.age.table.1$age2.dat[j] < age.group.15.25[2]){
        #           
        #           men.25.40.women.15.25.1 <- c(men.25.40.women.15.25.1, men.age.table.1$age2.dat[j])
        #           
        #         }else if(men.age.table.1$age2.dat[j] >= age.group.25.40[1] & men.age.table.1$age2.dat[j] < age.group.25.40[2]){
        #           
        #           men.25.40.women.25.40.1 <- c(men.25.40.women.25.40.1, men.age.table.1$age2.dat[j])
        #           
        #         }else if (men.age.table.1$age2.dat[j] >= age.group.40.50[1] & men.age.table.1$age2.dat[j] < age.group.40.50[2]){
        #           
        #           men.25.40.women.40.50.1 <- c(men.25.40.women.40.50.1, men.age.table.1$age2.dat[j])
        #         }
        #         
        #       }
        #       
        #     }
        #     
        #   }
        #   
        #   
        #   
        #   
        #   
        #   # women 25.40 and men
        #   
        #   women.25.40.men.15.25.2 <- vector()
        #   women.25.40.men.25.40.2 <- vector()
        #   women.25.40.men.40.50.2 <- vector()
        #   
        #   if(nrow(women.age.table.1) >1 ){
        #     
        #     for (j in 1:nrow(women.age.table.1)) {
        #       
        #       
        #       if(women.age.table.1$age1.dat[j] >= age.group.25.40[1] & women.age.table.1$age1.dat[j] < age.group.25.40[2]){
        #         
        #         if(women.age.table.1$age2.dat[j] >= age.group.15.25[1] & women.age.table.1$age2.dat[j] < age.group.15.25[2]){
        #           
        #           women.25.40.men.15.25.2 <- c(women.25.40.men.15.25.2, women.age.table.1$age2.dat[j])
        #           
        #         }else if(women.age.table.1$age2.dat[j] >= age.group.25.40[1] & women.age.table.1$age2.dat[j] < age.group.25.40[2]){
        #           
        #           women.25.40.men.25.40.2 <- c(women.25.40.men.25.40.2, women.age.table.1$age2.dat[j])
        #           
        #         }else if (women.age.table.1$age2.dat[j] >= age.group.40.50[1] & women.age.table.1$age2.dat[j] < age.group.40.50[2]){
        #           
        #           women.25.40.men.40.50.2 <- c(women.25.40.men.40.50.2, women.age.table.1$age2.dat[j])
        #         }
        #         
        #       }
        #       
        #       
        #     }
        #     
        #   }
        #   
        #   
        #   
        #   
        #   # men 40.50 and women
        #   
        #   men.40.50.women.15.25.1 <- vector()
        #   men.40.50.women.25.40.1 <- vector()
        #   men.40.50.women.40.50.1 <- vector()
        #   
        #   
        #   if(nrow(men.age.table.1) >1 ){
        #     
        #     for (j in 1:nrow(men.age.table.1)) {
        #       
        #       
        #       if(men.age.table.1$age1.dat[j] >= age.group.40.50[1] & men.age.table.1$age1.dat[j] < age.group.40.50[2]){
        #         
        #         if(men.age.table.1$age2.dat[j] >= age.group.15.25[1] & men.age.table.1$age2.dat[j] < age.group.15.25[2]){
        #           
        #           men.40.50.women.15.25.1 <- c(men.40.50.women.15.25.1, men.age.table.1$age2.dat[j])
        #           
        #         }else if(men.age.table.1$age2.dat[j] >= age.group.25.40[1] & men.age.table.1$age2.dat[j] < age.group.25.40[2]){
        #           
        #           men.40.50.women.25.40.1 <- c(men.40.50.women.25.40.1, men.age.table.1$age2.dat[j])
        #           
        #         }else if (men.age.table.1$age2.dat[j] >= age.group.40.50[1] & men.age.table.1$age2.dat[j] < age.group.40.50[2]){
        #           
        #           men.40.50.women.40.50.1 <- c(men.40.50.women.40.50.1, men.age.table.1$age2.dat[j])
        #         }
        #         
        #       }
        #       
        #       
        #     }
        #   }
        #   
        #   
        #   
        #   
        #   
        #   # women 40.50 and men
        #   
        #   women.40.50.men.15.25.2 <- vector()
        #   women.40.50.men.25.40.2 <- vector()
        #   women.40.50.men.40.50.2 <- vector()
        #   
        #   
        #   if(nrow(women.age.table.1) > 1 ){
        #     
        #     for (j in 1:nrow(women.age.table.1)) {
        #       
        #       
        #       if(women.age.table.1$age1.dat[j] >= age.group.40.50[1] & women.age.table.1$age1.dat[j] < age.group.40.50[2]){
        #         
        #         if(women.age.table.1$age2.dat[j] >= age.group.15.25[1] & women.age.table.1$age2.dat[j] < age.group.15.25[2]){
        #           
        #           women.40.50.men.15.25.2 <- c(women.40.50.men.15.25.2, women.age.table.1$age2.dat[j])
        #           
        #         }else if(women.age.table.1$age2.dat[j] >= age.group.25.40[1] & women.age.table.1$age2.dat[j] < age.group.25.40[2]){
        #           
        #           women.40.50.men.25.40.2 <- c(women.40.50.men.25.40.2, women.age.table.1$age2.dat[j])
        #           
        #         }else if (women.age.table.1$age2.dat[j] >= age.group.40.50[1] & women.age.table.1$age2.dat[j] < age.group.40.50[2]){
        #           
        #           women.40.50.men.40.50.2 <- c(women.40.50.men.40.50.2, women.age.table.1$age2.dat[j])
        #         }
        #         
        #       }
        #       
        #       
        #     }
        #     
        #   }
        #   
        #   
        #   
        #   
        #   
        #   men.15.25.women.15.25 <- c(men.15.25.women.15.25.1, women.15.25.men.15.25.2)
        #   
        #   men.15.25.women.25.40 <- c(men.15.25.women.25.40.1, women.25.40.men.15.25.2)
        #   
        #   men.15.25.women.40.50 <- c(men.15.25.women.40.50.1, women.40.50.men.15.25.2)
        #   
        #   men.25.40.women.15.25 <- c(men.25.40.women.15.25.1, women.15.25.men.25.40.2)
        #   
        #   men.25.40.women.25.40 <- c(men.25.40.women.25.40.1, women.25.40.men.25.40.2)
        #   
        #   men.25.40.women.40.50 <- c(men.25.40.women.40.50.1, women.40.50.men.25.40.2)
        #   
        #   men.40.50.women.15.25 <- c(men.40.50.women.15.25.1, women.15.25.men.40.50.2)
        #   
        #   men.40.50.women.25.40 <- c(men.40.50.women.25.40.1, women.25.40.men.40.50.2)
        #   
        #   men.40.50.women.40.50 <- c(men.40.50.women.40.50.1, women.40.50.men.40.50.2)
        #   
        #   
        #   Age.groups.table <- matrix(c(length(men.15.25.women.15.25), length(men.15.25.women.25.40), length(men.15.25.women.40.50),
        #                                length(men.25.40.women.15.25), length(men.25.40.women.25.40), length(men.25.40.women.40.50),
        #                                length(men.40.50.women.15.25), length(men.40.50.women.25.40), length(men.40.50.women.40.50)),
        #                              ncol = 3,
        #                              byrow = TRUE)
        #   
        #   colnames(Age.groups.table) <- c("Female.15.25", "Female.25.40", "Female.40.50")
        #   rownames(Age.groups.table) <- c("Male.15.25", "Male.25.40", "Male.40.50")
        #   
        #   Age.groups.table <- as.table(Age.groups.table)
        #   
        #   
        #   men.15.25.T <- sum(length(men.15.25.women.15.25), length(men.15.25.women.25.40), length(men.15.25.women.40.50))
        #   men.25.40.T <- sum(length(men.25.40.women.15.25), length(men.25.40.women.25.40), length(men.25.40.women.40.50))
        #   men.40.50.T <- sum(length(men.40.50.women.15.25), length(men.40.50.women.25.40), length(men.40.50.women.40.50))
        #   
        #   prop.men.age.groups.table <- matrix(c(length(men.15.25.women.15.25)/men.15.25.T, length(men.15.25.women.25.40)/men.15.25.T, length(men.15.25.women.40.50)/men.15.25.T,
        #                                         length(men.25.40.women.15.25)/men.25.40.T, length(men.25.40.women.25.40)/men.25.40.T, length(men.25.40.women.40.50)/men.25.40.T,
        #                                         length(men.40.50.women.15.25)/men.40.50.T, length(men.40.50.women.25.40)/men.40.50.T, length(men.40.50.women.40.50)/men.40.50.T),
        #                                       ncol = 3,
        #                                       byrow = TRUE)
        #   
        #   colnames(prop.men.age.groups.table) <- c("Female.15.25", "Female.25.40", "Female.40.50")
        #   rownames(prop.men.age.groups.table) <- c("prop.Male.15.25", "prop.Male.25.40", "prop.Male.40.50")
        #   
        #   
        #   
        #   
        #   women.15.25.T <- sum(length(men.15.25.women.15.25), length(men.25.40.women.15.25), length(men.40.50.women.15.25))
        #   women.25.40.T <- sum(length(men.15.25.women.25.40), length(men.25.40.women.25.40), length(men.40.50.women.25.40))
        #   women.40.50.T <- sum(length(men.15.25.women.40.50), length(men.25.40.women.40.50), length(men.40.50.women.40.50))
        #   
        #   prop.women.age.groups.table <- matrix(c(length(men.15.25.women.15.25)/women.15.25.T, length(men.25.40.women.15.25)/women.15.25.T, length(men.40.50.women.15.25)/women.15.25.T,
        #                                           length(men.15.25.women.25.40)/women.25.40.T, length(men.25.40.women.25.40)/women.25.40.T, length(men.40.50.women.25.40)/women.25.40.T,
        #                                           length(men.15.25.women.40.50)/women.40.50.T, length(men.25.40.women.40.50)/women.40.50.T, length(men.40.50.women.40.50)/women.40.50.T),
        #                                         ncol = 3,
        #                                         byrow = TRUE)
        #   
        #   colnames(prop.women.age.groups.table) <- c("Male.15.25", "Male.25.40", "Male.40.50")
        #   rownames(prop.women.age.groups.table) <- c("prop.Female.15.25", "prop.Female.25.40", "prop.Female.40.50")
        #   
        #   
        #   # Directionality
        #   
        #   men.15.25.women.15.25.MtoW <- c(men.15.25.women.15.25.1)
        #   
        #   men.15.25.women.25.40.MtoW <- c(men.15.25.women.25.40.1)
        #   
        #   men.15.25.women.40.50.MtoW <- c(men.15.25.women.40.50.1)
        #   
        #   men.25.40.women.15.25.MtoW <- c(men.25.40.women.15.25.1)
        #   
        #   men.25.40.women.25.40.MtoW <- c(men.25.40.women.25.40.1)
        #   
        #   men.25.40.women.40.50.MtoW <- c(men.25.40.women.40.50.1)
        #   
        #   men.40.50.women.15.25.MtoW <- c(men.40.50.women.15.25.1)
        #   
        #   men.40.50.women.25.40.MtoW <- c(men.40.50.women.25.40.1)
        #   
        #   men.40.50.women.40.50.MtoW <- c(men.40.50.women.40.50.1)
        #   
        #   
        #   Age.groups.table.MtoW <- matrix(c(length(men.15.25.women.15.25.MtoW), length(men.15.25.women.25.40.MtoW), length(men.15.25.women.40.50.MtoW),
        #                                     length(men.25.40.women.15.25.MtoW), length(men.25.40.women.25.40.MtoW), length(men.25.40.women.40.50.MtoW),
        #                                     length(men.40.50.women.15.25.MtoW), length(men.40.50.women.25.40.MtoW), length(men.40.50.women.40.50.MtoW)),
        #                                   ncol = 3,
        #                                   byrow = TRUE)
        #   
        #   colnames(Age.groups.table.MtoW) <- c("Female.15.25.MtoW", "Female.25.40.MtoW", "Female.40.50.MtoW")
        #   rownames(Age.groups.table.MtoW) <- c("Male.15.25.MtoW", "Male.25.40.MtoW", "Male.40.50.MtoW")
        #   
        #   Age.groups.table.MtoW <- as.table(Age.groups.table.MtoW)
        #   
        #   
        #   men.15.25.T.MtoW <- sum(length(men.15.25.women.15.25.MtoW), length(men.15.25.women.25.40.MtoW), length(men.15.25.women.40.50.MtoW))
        #   men.25.40.T.MtoW <- sum(length(men.25.40.women.15.25.MtoW), length(men.25.40.women.25.40.MtoW), length(men.25.40.women.40.50.MtoW))
        #   men.40.50.T.MtoW <- sum(length(men.40.50.women.15.25.MtoW), length(men.40.50.women.25.40.MtoW), length(men.40.50.women.40.50.MtoW))
        #   
        #   prop.men.age.groups.table.MtoW <- matrix(c(length(men.15.25.women.15.25.MtoW)/men.15.25.T.MtoW, length(men.15.25.women.25.40.MtoW)/men.15.25.T.MtoW, length(men.15.25.women.40.50.MtoW)/men.15.25.T.MtoW,
        #                                              length(men.25.40.women.15.25.MtoW)/men.25.40.T.MtoW, length(men.25.40.women.25.40.MtoW)/men.25.40.T.MtoW, length(men.25.40.women.40.50.MtoW)/men.25.40.T.MtoW,
        #                                              length(men.40.50.women.15.25.MtoW)/men.40.50.T.MtoW, length(men.40.50.women.25.40.MtoW)/men.40.50.T.MtoW, length(men.40.50.women.40.50.MtoW)/men.40.50.T.MtoW),
        #                                            ncol = 3,
        #                                            byrow = TRUE)
        #   
        #   colnames(prop.men.age.groups.table.MtoW) <- c("Female.15.25.MtoW", "Female.25.40.MtoW", "Female.40.50.MtoW")
        #   rownames(prop.men.age.groups.table.MtoW) <- c("prop.Male.15.25.MtoW", "prop.Male.25.40.MtoW", "prop.Male.40.50.MtoW")
        #   
        #   
        #   
        #   men.15.25.women.15.25.WtoM <- c(women.15.25.men.15.25.2)
        #   
        #   men.15.25.women.25.40.WtoM <- c(women.25.40.men.15.25.2)
        #   
        #   men.15.25.women.40.50.WtoM <- c(women.40.50.men.15.25.2)
        #   
        #   men.25.40.women.15.25.WtoM <- c( women.15.25.men.25.40.2)
        #   
        #   men.25.40.women.25.40.WtoM <- c(women.25.40.men.25.40.2)
        #   
        #   men.25.40.women.40.50.WtoM <- c(women.40.50.men.25.40.2)
        #   
        #   men.40.50.women.15.25.WtoM <- c(women.15.25.men.40.50.2)
        #   
        #   men.40.50.women.25.40.WtoM <- c(women.25.40.men.40.50.2)
        #   
        #   men.40.50.women.40.50.WtoM <- c(women.40.50.men.40.50.2)
        #   
        #   Age.groups.table.WtoM <- matrix(c(length(men.15.25.women.15.25.WtoM), length(men.15.25.women.25.40.WtoM), length(men.15.25.women.40.50.WtoM),
        #                                     length(men.25.40.women.15.25.WtoM), length(men.25.40.women.25.40.WtoM), length(men.25.40.women.40.50.WtoM),
        #                                     length(men.40.50.women.15.25.WtoM), length(men.40.50.women.25.40.WtoM), length(men.40.50.women.40.50.WtoM)),
        #                                   ncol = 3,
        #                                   byrow = TRUE)
        #   
        #   colnames(Age.groups.table.WtoM) <- c("Female.15.25.WtoM", "Female.25.40.WtoM", "Female.40.50.WtoM")
        #   rownames(Age.groups.table.WtoM) <- c("Male.15.25.WtoM", "Male.25.40.WtoM", "Male.40.50.WtoM")
        #   
        #   Age.groups.table.WtoM <- as.table(Age.groups.table.WtoM)
        #   
        #   
        #   men.15.25.T.WtoM <- sum(length(men.15.25.women.15.25.WtoM), length(men.15.25.women.25.40.WtoM), length(men.15.25.women.40.50.WtoM))
        #   men.25.40.T.WtoM <- sum(length(men.25.40.women.15.25.WtoM), length(men.25.40.women.25.40.WtoM), length(men.25.40.women.40.50.WtoM))
        #   men.40.50.T.WtoM <- sum(length(men.40.50.women.15.25.WtoM), length(men.40.50.women.25.40.WtoM), length(men.40.50.women.40.50.WtoM))
        #   
        #   prop.men.age.groups.table.WtoM <- matrix(c(length(men.15.25.women.15.25.WtoM)/men.15.25.T.WtoM, length(men.15.25.women.25.40.WtoM)/men.15.25.T.WtoM, length(men.15.25.women.40.50.WtoM)/men.15.25.T.WtoM,
        #                                              length(men.25.40.women.15.25.WtoM)/men.25.40.T.WtoM, length(men.25.40.women.25.40.WtoM)/men.25.40.T.WtoM, length(men.25.40.women.40.50.WtoM)/men.25.40.T.WtoM,
        #                                              length(men.40.50.women.15.25.WtoM)/men.40.50.T.WtoM, length(men.40.50.women.25.40.WtoM)/men.40.50.T.WtoM, length(men.40.50.women.40.50.WtoM)/men.40.50.T.WtoM),
        #                                            ncol = 3,
        #                                            byrow = TRUE)
        #   
        #   colnames(prop.men.age.groups.table.WtoM) <- c("Female.15.25.WtoM", "Female.25.40.WtoM", "Female.40.50.WtoM")
        #   rownames(prop.men.age.groups.table.WtoM) <- c("prop.Male.15.25.WtoM", "prop.Male.25.40.WtoM", "prop.Male.40.50.WtoM")
        #   
        #   
        #   
        #   outputlist <- NULL
        #   outputlist$Age.groups.table <- Age.groups.table
        #   outputlist$prop.men.age.groups.table <- prop.men.age.groups.table
        #   outputlist$prop.women.age.groups.table <- prop.women.age.groups.table
        #   
        #   outputlist$Age.groups.table.MtoW <- Age.groups.table.MtoW
        #   outputlist$prop.men.age.groups.table.MtoW <- prop.men.age.groups.table.MtoW
        #   outputlist$Age.groups.table.WtoM <- Age.groups.table.WtoM
        #   outputlist$prop.men.age.groups.table.WtoM <- prop.men.age.groups.table.WtoM
        #   
        #   
        #   outputlist$numbers.individuals.age.groups <- numbers.individuals.age.groups
        #   outputlist$mean.AD.age.groups <- mean.AD.age.groups
        #   outputlist$med.AD.age.groups <- med.AD.age.groups
        #   outputlist$sd.AD.age.groups <- sd.AD.age.groups
        #   
        #   
        #   return(outputlist)
        #   
        # }
        # 
        # 
        # 
        # 
        # 
        # 
        # 
        # # Results
        # ###########
        # 
        # # Function to handle NAs
        # 
        # 
        # NA.handle.fun <- function(input=input){
        #   
        #   v.names <- names(input)
        #   
        #   v <- as.numeric(input)
        #   
        #   v.vec <- vector()
        #   
        #   for(i in 1:length(v)){
        #     
        #     v.i <- v[i]
        #     
        #     if(is.na(v.i)==TRUE){
        #       v.j <- 0
        #     }else{
        #       v.j <- v.i
        #     }
        #     v.vec <- c(v.vec, v.j)
        #   }
        #   
        #   names(v.vec) <- v.names
        #   return(v.vec)
        # }
        # 
        # age.group.40.50 <- c(40, 50)
        # 
        # mCAr.IDs.100 <- IDs.Seq.Random(simpact.trans.net = simpact.trans.net.adv, # simpact.trans.net 
        #                                limitTransmEvents = 7,
        #                                timewindow = c(35, 40), 
        #                                seq.cov = 100, 
        #                                age.limit = age.group.40.50[2])
        # 
        # 
        # 
        # # table.simpact.trans.net.cov.100
        # 
        # infectionTable <- vector("list", length(simpact.trans.net.adv))
        # 
        # for(j in 1:length(simpact.trans.net.adv)){
        #   
        #   p <- j
        #   
        #   trans.network.i <- as.data.frame(simpact.trans.net.adv[[p]])
        #   
        #   # trans.network.i <- trans.network.i[-1,]
        #   
        #   id.lab <- paste0(p,".",trans.network.i$id,".C")
        #   
        #   trans.network.i$id.lab <- id.lab
        #   trans.network.i$ageSampTimeRec <- trans.network.i$SampTime - trans.network.i$TOBRec
        #   
        #   infectionTable[[p]] <- trans.network.i
        #   
        #   
        # }
        # 
        # 
        # infecttable <- rbindlist(infectionTable)
        # 
        # 
        # table.simpact.trans.net.adv <- infecttable # rbindlist(simpact.trans.net.adv)
        # 
        # table.simpact.trans.net.cov.100 <- dplyr::filter(table.simpact.trans.net.adv, table.simpact.trans.net.adv$id.lab%in%mCAr.IDs.100)
        # 
        # 
        # 
        # age.structure.transm.net.true.100.List <- age.groups.filtered.transmission.net.fun(table.transmission.net.cov = table.simpact.trans.net.cov.100,
        #                                                                                    age.group.15.25 = c(15,25),
        #                                                                                    age.group.25.40 = c(25,40),
        #                                                                                    age.group.40.50 = c(40,50))
        # 
        # 
        # age.struc.trans.net.true.cov.100 <- age.structure.transm.net.true.100.List$Age.groups.table
        # 
        # cov.100.age.str.M.15.25.F.15.25 <- age.struc.trans.net.true.cov.100[1,][1]
        # cov.100.age.str.M.25.40.F.15.25 <- age.struc.trans.net.true.cov.100[2,][1]
        # cov.100.age.str.M.40.50.F.15.25 <- age.struc.trans.net.true.cov.100[3,][1]
        # 
        # cov.100.age.str.M.15.25.F.25.40 <- age.struc.trans.net.true.cov.100[1,][2]
        # cov.100.age.str.M.25.40.F.25.40 <- age.struc.trans.net.true.cov.100[2,][2]
        # cov.100.age.str.M.40.50.F.25.40 <- age.struc.trans.net.true.cov.100[3,][2]
        # 
        # cov.100.age.str.M.15.25.F.40.50 <- age.struc.trans.net.true.cov.100[1,][3]
        # cov.100.age.str.M.25.40.F.40.50 <- age.struc.trans.net.true.cov.100[2,][3]
        # cov.100.age.str.M.40.50.F.40.50 <- age.struc.trans.net.true.cov.100[3,][3]
        # 
        # table.cov.100.age.str <- c(cov.100.age.str.M.15.25.F.15.25, cov.100.age.str.M.25.40.F.15.25, cov.100.age.str.M.40.50.F.15.25,
        #                            cov.100.age.str.M.15.25.F.25.40, cov.100.age.str.M.25.40.F.25.40, cov.100.age.str.M.40.50.F.25.40,
        #                            cov.100.age.str.M.15.25.F.40.50, cov.100.age.str.M.25.40.F.40.50, cov.100.age.str.M.40.50.F.40.50)
        # 
        # 
        # names(table.cov.100.age.str) <- c("cov.100.M.15.25.F.15.25", "cov.100.M.25.40.F.15.25", "cov.100.M.40.50.F.15.25",
        #                                   "cov.100.M.15.25.F.25.40", "cov.100.M.25.40.F.25.40", "cov.100.M.40.50.F.25.40",
        #                                   "cov.100.M.15.25.F.40.50", "cov.100.M.25.40.F.40.50", "cov.100.M.40.50.F.40.50")
        # 
        # 
        # # Men prop
        # 
        # age.struc.trans.net.true.cov.100.prop.men <- age.structure.transm.net.true.100.List$prop.men.age.groups.table
        # 
        # cov.100.true.age.str.prop.men.15.25.F.15.25 <- age.struc.trans.net.true.cov.100.prop.men[1,][1]
        # cov.100.true.age.str.prop.men.25.40.F.15.25 <- age.struc.trans.net.true.cov.100.prop.men[2,][1]
        # cov.100.true.age.str.prop.men.40.50.F.15.25 <- age.struc.trans.net.true.cov.100.prop.men[3,][1]
        # 
        # cov.100.true.age.str.prop.men.15.25.F.25.40 <- age.struc.trans.net.true.cov.100.prop.men[1,][2]
        # cov.100.true.age.str.prop.men.25.40.F.25.40 <- age.struc.trans.net.true.cov.100.prop.men[2,][2]
        # cov.100.true.age.str.prop.men.40.50.F.25.40 <- age.struc.trans.net.true.cov.100.prop.men[3,][2]
        # 
        # cov.100.true.age.str.prop.men.15.25.F.40.50 <- age.struc.trans.net.true.cov.100.prop.men[1,][3]
        # cov.100.true.age.str.prop.men.25.40.F.40.50 <- age.struc.trans.net.true.cov.100.prop.men[2,][3]
        # cov.100.true.age.str.prop.men.40.50.F.40.50 <- age.struc.trans.net.true.cov.100.prop.men[3,][3]
        # 
        # table.cov.100.true.age.str.prop.men <- c(cov.100.true.age.str.prop.men.15.25.F.15.25, cov.100.true.age.str.prop.men.25.40.F.15.25, cov.100.true.age.str.prop.men.40.50.F.15.25,
        #                                          cov.100.true.age.str.prop.men.15.25.F.25.40, cov.100.true.age.str.prop.men.25.40.F.25.40, cov.100.true.age.str.prop.men.40.50.F.25.40,
        #                                          cov.100.true.age.str.prop.men.15.25.F.40.50, cov.100.true.age.str.prop.men.25.40.F.40.50, cov.100.true.age.str.prop.men.40.50.F.40.50)
        # 
        # names(table.cov.100.true.age.str.prop.men) <- c("cov.100.true.prop.men15.25.F.15.25", "cov.100.true.prop.men25.40.F.15.25", "cov.100.true.prop.men40.50.F.15.25",
        #                                                 "cov.100.true.prop.men15.25.F.25.40", "cov.100.true.prop.men25.40.F.25.40", "cov.100.true.prop.men40.50.F.25.40",
        #                                                 "cov.100.true.prop.men15.25.F.40.50", "cov.100.true.prop.men25.40.F.40.50", "cov.100.true.prop.men40.50.F.40.50")
        # 
        # 
        # table.cov.100.true.age.str.prop.men <- NA.handle.fun(input = table.cov.100.true.age.str.prop.men)
        # 
        # 
        # 
        # # Women prop
        # 
        # age.structure.transm.clust.true.prop.women <- age.structure.transm.net.true.100.List$prop.women.age.groups.table
        # 
        # cov.100.true.age.str.prop.women.15.25.M.15.25 <- age.structure.transm.clust.true.prop.women[1,][1]
        # cov.100.true.age.str.prop.women.25.40.M.15.25 <- age.structure.transm.clust.true.prop.women[2,][1]
        # cov.100.true.age.str.prop.women.40.50.M.15.25 <- age.structure.transm.clust.true.prop.women[3,][1]
        # 
        # cov.100.true.age.str.prop.women.15.25.M.25.40 <- age.structure.transm.clust.true.prop.women[1,][2]
        # cov.100.true.age.str.prop.women.25.40.M.25.40 <- age.structure.transm.clust.true.prop.women[2,][2]
        # cov.100.true.age.str.prop.women.40.50.M.25.40 <- age.structure.transm.clust.true.prop.women[3,][2]
        # 
        # cov.100.true.age.str.prop.women.15.25.M.40.50 <- age.structure.transm.clust.true.prop.women[1,][3]
        # cov.100.true.age.str.prop.women.25.40.M.40.50 <- age.structure.transm.clust.true.prop.women[2,][3]
        # cov.100.true.age.str.prop.women.40.50.M.40.50 <- age.structure.transm.clust.true.prop.women[3,][3]
        # 
        # table.cov.100.true.age.str.prop.women <- c(cov.100.true.age.str.prop.women.15.25.M.15.25, cov.100.true.age.str.prop.women.25.40.M.15.25, cov.100.true.age.str.prop.women.40.50.M.15.25,
        #                                            cov.100.true.age.str.prop.women.15.25.M.25.40, cov.100.true.age.str.prop.women.25.40.M.25.40, cov.100.true.age.str.prop.women.40.50.M.25.40,
        #                                            cov.100.true.age.str.prop.women.15.25.M.40.50, cov.100.true.age.str.prop.women.25.40.M.40.50, cov.100.true.age.str.prop.women.40.50.M.40.50)
        # 
        # names(table.cov.100.true.age.str.prop.women) <- c("cov.100.true.prop.women15.25.M.15.25", "cov.100.true.prop.women25.40.M.15.25", "cov.100.true.prop.women40.50.M.15.25",
        #                                                   "cov.100.true.prop.women15.25.M.25.40", "cov.100.true.prop.women25.40.M.25.40", "cov.100.true.prop.women40.50.M.25.40",
        #                                                   "cov.100.true.prop.women15.25.M.40.50", "cov.100.true.prop.women25.40.M.40.50", "cov.100.true.prop.women40.50.M.40.50")
        # 
        # 
        # 
        # table.cov.100.true.age.str.prop.women <- NA.handle.fun(input = table.cov.100.true.age.str.prop.women)
        # 
        # 
        # 
        # numbers.individuals.age.groups.cov.100 <- age.structure.transm.net.true.100.List$numbers.individuals.age.groups
        # mean.AD.age.groups.cov.100 <- age.structure.transm.net.true.100.List$mean.AD.age.groups 
        # med.AD.age.groups.cov.100 <- age.structure.transm.net.true.100.List$med.AD.age.groups 
        # sd.AD.age.groups.cov.100 <- age.structure.transm.net.true.100.List$sd.AD.age.groups 
        # 
        # names(numbers.individuals.age.groups.cov.100) <- paste0("cov.100.", names(numbers.individuals.age.groups.cov.100))
        # names(mean.AD.age.groups.cov.100) <- paste0("cov.100.", names(mean.AD.age.groups.cov.100))
        # names(med.AD.age.groups.cov.100) <- paste0("cov.100.", names(med.AD.age.groups.cov.100))
        # names(sd.AD.age.groups.cov.100) <- paste0("cov.100.", names(sd.AD.age.groups.cov.100))
        # 
        # 
        # 
        # 
        # cov.100.stats <- c(table.cov.100.age.str, table.cov.100.true.age.str.prop.men, table.cov.100.true.age.str.prop.women,
        #                    numbers.individuals.age.groups.cov.100, 
        #                    mean.AD.age.groups.cov.100, 
        #                    med.AD.age.groups.cov.100,
        #                    sd.AD.age.groups.cov.100)
        # 
        # 
        
        # MCAR
        
        

        
        
        
        res.clust.MCAR.cov.35 <- tryCatch(age.mixing.MCAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 35,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                          error=function(e) return(rep(NA, 111)))
        
        res.clust.MCAR.cov.40 <- tryCatch(age.mixing.MCAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 40,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                          error=function(e) return(rep(NA, 111)))
        
        res.clust.MCAR.cov.45 <- tryCatch(age.mixing.MCAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 45,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                          error=function(e) return(rep(NA, 111)))
        
        res.clust.MCAR.cov.50 <- tryCatch(age.mixing.MCAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 50,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                          error=function(e) return(rep(NA, 111)))
        
        res.clust.MCAR.cov.55 <- tryCatch(age.mixing.MCAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 55,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                          error=function(e) return(rep(NA, 111)))
        
        
        res.clust.MCAR.cov.60 <- tryCatch(age.mixing.MCAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 60,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                          error=function(e) return(rep(NA, 111)))
        
        
        res.clust.MCAR.cov.65 <- tryCatch(age.mixing.MCAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 65,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                          error=function(e) return(rep(NA, 111)))
        
        
        
        
        res.clust.MCAR.cov.70 <- tryCatch(age.mixing.MCAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 70,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                          error=function(e) return(rep(NA, 111)))
        
        
        
        res.clust.MCAR.cov.75 <- tryCatch(age.mixing.MCAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 75,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                          error=function(e) return(rep(NA, 111)))
        
        
        
        
        res.clust.MCAR.cov.80 <- tryCatch(age.mixing.MCAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 80,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                          error=function(e) return(rep(NA, 111)))
        
        res.clust.MCAR.cov.85 <- tryCatch(age.mixing.MCAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 85,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                          error=function(e) return(rep(NA, 111)))
        
        
        
        
        res.clust.MCAR.cov.90 <- tryCatch(age.mixing.MCAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 90,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                          error=function(e) return(rep(NA, 111)))
        
        res.clust.MCAR.cov.95 <- tryCatch(age.mixing.MCAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 95,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                          error=function(e) return(rep(NA, 111)))
        
        # 100 Coverage
        
        res.clust.MCAR.cov.100 <- tryCatch(age.mixing.MCAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                               datalist.agemix = datalist.agemix, 
                                                               work.dir = work.dir,  
                                                               dirfasttree = work.dir, 
                                                               sub.dir.rename = sub.dir.rename,
                                                               limitTransmEvents = 7,
                                                               timewindow = c(35, 40),
                                                               seq.cov = 100,
                                                               age.group.15.25 = c(15,25),
                                                               age.group.25.40 = c(25,40),
                                                               age.group.40.50 = c(40,50),
                                                               cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        
        names.columns <- names(res.clust.MCAR.cov.100)
        
        
        results.mcar <- as.numeric(c(res.clust.MCAR.cov.35, res.clust.MCAR.cov.40, 
                                     res.clust.MCAR.cov.45, res.clust.MCAR.cov.50, 
                                     res.clust.MCAR.cov.55, res.clust.MCAR.cov.60, 
                                     res.clust.MCAR.cov.65, res.clust.MCAR.cov.70, 
                                     res.clust.MCAR.cov.75, res.clust.MCAR.cov.80, 
                                     res.clust.MCAR.cov.85, res.clust.MCAR.cov.90, 
                                     res.clust.MCAR.cov.95, res.clust.MCAR.cov.100))
        
        
        names(results.mcar) <- c(paste0("cov.MCAR.",35,".",paste0(names.columns)), paste0("cov.MCAR.",40,".",paste0(names.columns)),
                                 paste0("cov.MCAR.",45,".",paste0(names.columns)), paste0("cov.MCAR.",50,".",paste0(names.columns)),
                                 paste0("cov.MCAR.",55,".",paste0(names.columns)), paste0("cov.MCAR.",60,".",paste0(names.columns)),
                                 paste0("cov.MCAR.",65,".",paste0(names.columns)), paste0("cov.MCAR.",70,".",paste0(names.columns)),
                                 paste0("cov.MCAR.",75,".",paste0(names.columns)), paste0("cov.MCAR.",80,".",paste0(names.columns)),
                                 paste0("cov.MCAR.",85,".",paste0(names.columns)), paste0("cov.MCAR.",90,".",paste0(names.columns)),
                                 paste0("cov.MCAR.",95,".",paste0(names.columns)), paste0("cov.MCAR.",100,".",paste0(names.columns)))
        
        
        # MAR
        # (a) 0.7
        
        
        
        res.clust.MAR.a.cov.35 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 35,
                                                              seq.gender.ratio = 0.7,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        res.clust.MAR.a.cov.40 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 40,
                                                              seq.gender.ratio = 0.7,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        res.clust.MAR.a.cov.45 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 45,
                                                              seq.gender.ratio = 0.7,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        res.clust.MAR.a.cov.50 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 50,
                                                              seq.gender.ratio = 0.7,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        res.clust.MAR.a.cov.55 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 55,
                                                              seq.gender.ratio = 0.7,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        
        res.clust.MAR.a.cov.60 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 60,
                                                              seq.gender.ratio = 0.7,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        
        res.clust.MAR.a.cov.65 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 65,
                                                              seq.gender.ratio = 0.7,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        
        
        
        res.clust.MAR.a.cov.70 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 70,
                                                              seq.gender.ratio = 0.7,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        
        
        res.clust.MAR.a.cov.75 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 75,
                                                              seq.gender.ratio = 0.7,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        
        
        
        res.clust.MAR.a.cov.80 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 80,
                                                              seq.gender.ratio = 0.7,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        res.clust.MAR.a.cov.85 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 85,
                                                              seq.gender.ratio = 0.7,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        
        
        
        res.clust.MAR.a.cov.90 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 90,
                                                              seq.gender.ratio = 0.7,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        res.clust.MAR.a.cov.95 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 95,
                                                              seq.gender.ratio = 0.7,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        
        
        results.mar.a <- as.numeric(c(res.clust.MAR.a.cov.35, res.clust.MAR.a.cov.40, 
                                      res.clust.MAR.a.cov.45, res.clust.MAR.a.cov.50, 
                                      res.clust.MAR.a.cov.55, res.clust.MAR.a.cov.60, 
                                      res.clust.MAR.a.cov.65, res.clust.MAR.a.cov.70, 
                                      res.clust.MAR.a.cov.75, res.clust.MAR.a.cov.80, 
                                      res.clust.MAR.a.cov.85,res.clust.MAR.a.cov.90, 
                                      res.clust.MAR.a.cov.95))
        

        names(results.mar.a) <- c(paste0("cov.MAR.a.",35,".",paste0(names.columns)), paste0("cov.MAR.a.",40,".",paste0(names.columns)),
                                 paste0("cov.MAR.a.",45,".",paste0(names.columns)), paste0("cov.MAR.a.",50,".",paste0(names.columns)),
                                 paste0("cov.MAR.a.",55,".",paste0(names.columns)), paste0("cov.MAR.a.",60,".",paste0(names.columns)),
                                 paste0("cov.MAR.a.",65,".",paste0(names.columns)), paste0("cov.MAR.a.",70,".",paste0(names.columns)),
                                 paste0("cov.MAR.a.",75,".",paste0(names.columns)), paste0("cov.MAR.a.",80,".",paste0(names.columns)),
                                 paste0("cov.MAR.a.",85,".",paste0(names.columns)), paste0("cov.MAR.a.",90,".",paste0(names.columns)),
                                 paste0("cov.MAR.a.",95,".",paste0(names.columns)))
        
        
        
        # MAR
        # (b) 0.3
        
        
        res.clust.MAR.b.cov.35 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 35,
                                                              seq.gender.ratio = 0.3,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        res.clust.MAR.b.cov.40 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 40,
                                                              seq.gender.ratio = 0.3,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        res.clust.MAR.b.cov.45 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 45,
                                                              seq.gender.ratio = 0.3,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        res.clust.MAR.b.cov.50 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 50,
                                                              seq.gender.ratio = 0.3,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        res.clust.MAR.b.cov.55 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 55,
                                                              seq.gender.ratio = 0.3,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        
        res.clust.MAR.b.cov.60 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 60,
                                                              seq.gender.ratio = 0.3,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        
        res.clust.MAR.b.cov.65 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 65,
                                                              seq.gender.ratio = 0.3,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        
        
        
        res.clust.MAR.b.cov.70 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 70,
                                                              seq.gender.ratio = 0.3,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        
        
        res.clust.MAR.b.cov.75 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 75,
                                                              seq.gender.ratio = 0.3,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        
        
        
        res.clust.MAR.b.cov.80 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 80,
                                                              seq.gender.ratio = 0.3,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        res.clust.MAR.b.cov.85 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 85,
                                                              seq.gender.ratio = 0.3,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        
        
        
        res.clust.MAR.b.cov.90 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 90,
                                                              seq.gender.ratio = 0.3,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        res.clust.MAR.b.cov.95 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 95,
                                                              seq.gender.ratio = 0.3,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        
        
        results.mar.b <- as.numeric(c(res.clust.MAR.b.cov.35, res.clust.MAR.b.cov.40, 
                                      res.clust.MAR.b.cov.45, res.clust.MAR.b.cov.50, 
                                      res.clust.MAR.b.cov.55, res.clust.MAR.b.cov.60, 
                                      res.clust.MAR.b.cov.65, res.clust.MAR.b.cov.70, 
                                      res.clust.MAR.b.cov.75, res.clust.MAR.b.cov.80, 
                                      res.clust.MAR.b.cov.85,res.clust.MAR.b.cov.90, 
                                      res.clust.MAR.b.cov.95))
        
        
        names(results.mar.b) <- c(paste0("cov.MAR.b.",35,".",paste0(names.columns)), paste0("cov.MAR.b.",40,".",paste0(names.columns)),
                                  paste0("cov.MAR.b.",45,".",paste0(names.columns)), paste0("cov.MAR.b.",50,".",paste0(names.columns)),
                                  paste0("cov.MAR.b.",55,".",paste0(names.columns)), paste0("cov.MAR.b.",60,".",paste0(names.columns)),
                                  paste0("cov.MAR.b.",65,".",paste0(names.columns)), paste0("cov.MAR.b.",70,".",paste0(names.columns)),
                                  paste0("cov.MAR.b.",75,".",paste0(names.columns)), paste0("cov.MAR.b.",80,".",paste0(names.columns)),
                                  paste0("cov.MAR.b.",85,".",paste0(names.columns)), paste0("cov.MAR.b.",90,".",paste0(names.columns)),
                                  paste0("cov.MAR.b.",95,".",paste0(names.columns)))
        
        # MAR
        # (c) 0.5
        
        
        res.clust.MAR.c.cov.35 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 35,
                                                              seq.gender.ratio = 0.5,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        res.clust.MAR.c.cov.40 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 40,
                                                              seq.gender.ratio = 0.5,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        res.clust.MAR.c.cov.45 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 45,
                                                              seq.gender.ratio = 0.5,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        res.clust.MAR.c.cov.50 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 50,
                                                              seq.gender.ratio = 0.5,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        res.clust.MAR.c.cov.55 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 55,
                                                              seq.gender.ratio = 0.5,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        
        res.clust.MAR.c.cov.60 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 60,
                                                              seq.gender.ratio = 0.5,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        
        res.clust.MAR.c.cov.65 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 65,
                                                              seq.gender.ratio = 0.5,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        
        
        
        res.clust.MAR.c.cov.70 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 70,
                                                              seq.gender.ratio = 0.5,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        
        
        res.clust.MAR.c.cov.75 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 75,
                                                              seq.gender.ratio = 0.5,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        
        
        
        res.clust.MAR.c.cov.80 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 80,
                                                              seq.gender.ratio = 0.5,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        res.clust.MAR.c.cov.85 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 85,
                                                              seq.gender.ratio = 0.5,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        
        
        
        res.clust.MAR.c.cov.90 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 90,
                                                              seq.gender.ratio = 0.5,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        res.clust.MAR.c.cov.95 <- tryCatch(age.mixing.MAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(35, 40),
                                                              seq.cov = 95,
                                                              seq.gender.ratio = 0.5,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                           error=function(e) return(rep(NA, 111)))
        
        
        
        results.mar.c <- as.numeric(c(res.clust.MAR.c.cov.35, res.clust.MAR.c.cov.40, 
                                      res.clust.MAR.c.cov.45, res.clust.MAR.c.cov.50, 
                                      res.clust.MAR.c.cov.55, res.clust.MAR.c.cov.60, 
                                      res.clust.MAR.c.cov.65, res.clust.MAR.c.cov.70, 
                                      res.clust.MAR.c.cov.75, res.clust.MAR.c.cov.80, 
                                      res.clust.MAR.c.cov.85,res.clust.MAR.c.cov.90, 
                                      res.clust.MAR.c.cov.95))
        
        
        names(results.mar.c) <- c(paste0("cov.MAR.c.",35,".",paste0(names.columns)), paste0("cov.MAR.c.",40,".",paste0(names.columns)),
                                  paste0("cov.MAR.c.",45,".",paste0(names.columns)), paste0("cov.MAR.c.",50,".",paste0(names.columns)),
                                  paste0("cov.MAR.c.",55,".",paste0(names.columns)), paste0("cov.MAR.c.",60,".",paste0(names.columns)),
                                  paste0("cov.MAR.c.",65,".",paste0(names.columns)), paste0("cov.MAR.c.",70,".",paste0(names.columns)),
                                  paste0("cov.MAR.c.",75,".",paste0(names.columns)), paste0("cov.MAR.c.",80,".",paste0(names.columns)),
                                  paste0("cov.MAR.c.",85,".",paste0(names.columns)), paste0("cov.MAR.c.",90,".",paste0(names.columns)),
                                  paste0("cov.MAR.c.",95,".",paste0(names.columns)))
        
        
        
        results.summary.epidemic.rels.df <- summary.epidemic.rels.df
        
        
        names(results.summary.epidemic.rels.df) <- paste0(names(summary.epidemic.rels.df))
        
        
        # cov.100.stats <- cov.100.stats
        # 
        # names(cov.100.stats) <- names(cov.100.stats)
        
        
        
        
        results.outputvector <- c(results.summary.epidemic.rels.df, 
                                  results.mcar, results.mar.a, results.mar.b, results.mar.c)
        
        
        
      }else{
        
        results.outputvector <- rep(NA, 5908) # 25 + 111 + 111*13*4
        
      }
      
    }
    
  }
  
  unlink(paste0(sub.dir.rename), recursive = TRUE)
  
  return(results.outputvector)
  
  
}


