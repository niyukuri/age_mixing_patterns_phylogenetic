# 
# # Statistical test: compairing missingness scenarios: MCAR & MAR

# Go to Line 10575 and load the comparison table

set.seed(777)

suppressMessages(library(dplyr))


### Load data 


# dr1 <- read.csv("/home/david/Dropbox/RedigerMaThese/Manuscript_1_age_mixing/results_april_2020/results.mcarmar.large.AD_280_seed_1.csv")
# dr2 <- read.csv("/home/david/Dropbox/RedigerMaThese/Manuscript_1_age_mixing/results_april_2020/results.mcarmar.large.AD_280_seed_285.csv")
# dr3 <- read.csv("/home/david/Dropbox/RedigerMaThese/Manuscript_1_age_mixing/results_april_2020/results.mcarmar.large.AD_280_seed_570.csv")
# dr4 <- read.csv("/home/david/Dropbox/RedigerMaThese/Manuscript_1_age_mixing/results_april_2020/results.mcarmar.large.AD_280_seed_855.csv")
# dr5 <- read.csv("/home/david/Dropbox/RedigerMaThese/Manuscript_1_age_mixing/results_april_2020/results.mcarmar.large.AD_280_seed_1135.csv")
# dr6 <- read.csv("/home/david/Dropbox/RedigerMaThese/Manuscript_1_age_mixing/results_april_2020/results.mcarmar.large.AD_280_seed_1420.csv")
# dr7 <- read.csv("/home/david/Dropbox/RedigerMaThese/Manuscript_1_age_mixing/results_april_2020/results.mcarmar.large.AD_280_seed_1705.csv")
# dr8 <- read.csv("/home/david/Dropbox/RedigerMaThese/Manuscript_1_age_mixing/results_april_2020/results.mcarmar.large.AD_280_seed_1990.csv")
# dr9 <- read.csv("/home/david/Dropbox/RedigerMaThese/Manuscript_1_age_mixing/results_april_2020/results.mcarmar.large.AD_280_seed_2275.csv")
# dr10 <- read.csv("/home/david/Dropbox/RedigerMaThese/Manuscript_1_age_mixing/results_april_2020/results.mcarmar.large.AD_280_seed_2560.csv")


# 
# dr <- rbind(dr1, dr2, dr3, dr4, dr5, dr6, dr7, dr8, dr9, dr10)

# write.csv(dr, file = "/home/david/age_mixing_patterns_phylogenetic/sim_outputs/april_2020_results.mcarmar.large.AD.csv")
# 
# dr <- read.csv("/home/david/age_mixing_patterns_phylogenetic/sim_outputs/F_results.mcarmar.large.AD.csv")


dr <- read.csv("/home/david/age_mixing_patterns_phylogenetic/sim_outputs/april_2020_results.mcarmar.large.AD.csv")



# Wilcox function for one by one comaprison of vectors of median values of metric A and B
# E.g.: proportions of women 15 - 24 phylogenetically linked to men between 40 - 49 in all sequence coverage 35 - 95%


wilcox.test.A.B <- function(mcar=mcar, mar=mar) {
  
  d <- data.frame(mcar,mar)
  r <- na.omit(d)
  di <- r$mcar - r$mar
  xx <- r

    mcar_mar <- wilcox.test(xx$mcar, xx$mar, paired = TRUE)
  
  V <- mcar_mar[1]$statistic[[1]]
  p_value <- mcar_mar[[3]]
  
  outputvector <-  mcar_mar[[3]] # c(V, p_value)
  
  return(outputvector)
  
}



# I. Proportions of  pairings ---------

# MCAR

d.MCAR <- dr %>%
  select(contains("MCAR."))

# 35

d.MCAR.cov.35 <- d.MCAR %>%
  select(contains("cov.MCAR.35.")) 
d.MACR.cov.35.cl.prop.men <-  d.MCAR.cov.35 %>%
  select(contains(".cl.prop.men"))
d.MACR.cov.35.cl.prop.women <-  d.MCAR.cov.35 %>%
  select(contains(".cl.prop.women"))

vector.MCAR.cov.35.cl.prop.men15.25.F.15.25 <- d.MACR.cov.35.cl.prop.men[,1]
vector.MCAR.cov.35.cl.prop.women15.25.M.15.25 <- d.MACR.cov.35.cl.prop.women[,1]

vector.MCAR.cov.35.cl.prop.men25.40.F.15.25 <- d.MACR.cov.35.cl.prop.men[,2]
vector.MCAR.cov.35.cl.prop.women15.25.M.25.40 <- d.MACR.cov.35.cl.prop.women[,4]

vector.MCAR.cov.35.cl.prop.men25.40.F.25.40 <- d.MACR.cov.35.cl.prop.men[,5]
vector.MCAR.cov.35.cl.prop.women25.40.M.25.40 <- d.MACR.cov.35.cl.prop.women[,5]

vector.MCAR.cov.35.cl.prop.men40.50.F.15.25 <- d.MACR.cov.35.cl.prop.men[,3]
vector.MCAR.cov.35.cl.prop.women15.25.M.40.50 <- d.MACR.cov.35.cl.prop.women[,7]

vector.MCAR.cov.35.cl.prop.men40.50.F.25.40 <- d.MACR.cov.35.cl.prop.men[,6]
vector.MCAR.cov.35.cl.prop.women25.40.M.40.50 <- d.MACR.cov.35.cl.prop.women[,8]




# 40

d.MCAR.cov.40 <- d.MCAR %>%
  select(contains("cov.MCAR.40.")) 
d.MACR.cov.40.cl.prop.men <-  d.MCAR.cov.40 %>%
  select(contains(".cl.prop.men"))
d.MACR.cov.40.cl.prop.women <-  d.MCAR.cov.40 %>%
  select(contains(".cl.prop.women"))


vector.MCAR.cov.40.cl.prop.men15.25.F.15.25 <- d.MACR.cov.40.cl.prop.men[,1]
vector.MCAR.cov.40.cl.prop.women15.25.M.15.25 <- d.MACR.cov.40.cl.prop.women[,1]

vector.MCAR.cov.40.cl.prop.men25.40.F.15.25 <- d.MACR.cov.40.cl.prop.men[,2]
vector.MCAR.cov.40.cl.prop.women15.25.M.25.40 <- d.MACR.cov.40.cl.prop.women[,4]

vector.MCAR.cov.40.cl.prop.men25.40.F.25.40 <- d.MACR.cov.40.cl.prop.men[,5]
vector.MCAR.cov.40.cl.prop.women25.40.M.25.40 <- d.MACR.cov.40.cl.prop.women[,5]


vector.MCAR.cov.40.cl.prop.men40.50.F.15.25 <- d.MACR.cov.40.cl.prop.men[,3]
vector.MCAR.cov.40.cl.prop.women15.25.M.40.50 <- d.MACR.cov.40.cl.prop.women[,7]

vector.MCAR.cov.40.cl.prop.men40.50.F.25.40 <- d.MACR.cov.40.cl.prop.men[,6]
vector.MCAR.cov.40.cl.prop.women25.40.M.40.50 <- d.MACR.cov.40.cl.prop.women[,8]


# 45

d.MCAR.cov.45 <- d.MCAR %>%
  select(contains("cov.MCAR.45.")) 
d.MACR.cov.45.cl.prop.men <-  d.MCAR.cov.45 %>%
  select(contains(".cl.prop.men"))
d.MACR.cov.45.cl.prop.women <-  d.MCAR.cov.45 %>%
  select(contains(".cl.prop.women"))


vector.MCAR.cov.45.cl.prop.men15.25.F.15.25 <- d.MACR.cov.45.cl.prop.men[,1]
vector.MCAR.cov.45.cl.prop.women15.25.M.15.25 <- d.MACR.cov.45.cl.prop.women[,1]

vector.MCAR.cov.45.cl.prop.men25.40.F.15.25 <- d.MACR.cov.45.cl.prop.men[,2]
vector.MCAR.cov.45.cl.prop.women15.25.M.25.40 <- d.MACR.cov.45.cl.prop.women[,4]

vector.MCAR.cov.45.cl.prop.men25.40.F.25.40 <- d.MACR.cov.45.cl.prop.men[,5]
vector.MCAR.cov.45.cl.prop.women25.40.M.25.40 <- d.MACR.cov.45.cl.prop.women[,5]

vector.MCAR.cov.45.cl.prop.men40.50.F.15.25 <- d.MACR.cov.45.cl.prop.men[,3]
vector.MCAR.cov.45.cl.prop.women15.25.M.40.50 <- d.MACR.cov.45.cl.prop.women[,7]

vector.MCAR.cov.45.cl.prop.men40.50.F.25.40 <- d.MACR.cov.45.cl.prop.men[,6]
vector.MCAR.cov.45.cl.prop.women25.40.M.40.50 <- d.MACR.cov.45.cl.prop.women[,8]


# 50

d.MCAR.cov.50 <- d.MCAR %>%
  select(contains("cov.MCAR.50.")) 
d.MACR.cov.50.cl.prop.men <-  d.MCAR.cov.50 %>%
  select(contains(".cl.prop.men"))
d.MACR.cov.50.cl.prop.women <-  d.MCAR.cov.50 %>%
  select(contains(".cl.prop.women"))


vector.MCAR.cov.50.cl.prop.men15.25.F.15.25 <- d.MACR.cov.50.cl.prop.men[,1]
vector.MCAR.cov.50.cl.prop.women15.25.M.15.25 <- d.MACR.cov.50.cl.prop.women[,1]

vector.MCAR.cov.50.cl.prop.men25.40.F.15.25 <- d.MACR.cov.50.cl.prop.men[,2]
vector.MCAR.cov.50.cl.prop.women15.25.M.25.40 <- d.MACR.cov.50.cl.prop.women[,4]

vector.MCAR.cov.50.cl.prop.men25.40.F.25.40 <- d.MACR.cov.50.cl.prop.men[,5]
vector.MCAR.cov.50.cl.prop.women25.40.M.25.40 <- d.MACR.cov.50.cl.prop.women[,5]

vector.MCAR.cov.50.cl.prop.men40.50.F.15.25 <- d.MACR.cov.50.cl.prop.men[,3]
vector.MCAR.cov.50.cl.prop.women15.25.M.40.50 <- d.MACR.cov.50.cl.prop.women[,7]

vector.MCAR.cov.50.cl.prop.men40.50.F.25.40 <- d.MACR.cov.50.cl.prop.men[,6]
vector.MCAR.cov.50.cl.prop.women25.40.M.40.50 <- d.MACR.cov.50.cl.prop.women[,8]



# 55

d.MCAR.cov.55 <- d.MCAR %>%
  select(contains("cov.MCAR.55.")) 
d.MACR.cov.55.cl.prop.men <-  d.MCAR.cov.55 %>%
  select(contains(".cl.prop.men"))
d.MACR.cov.55.cl.prop.women <-  d.MCAR.cov.55 %>%
  select(contains(".cl.prop.women"))


vector.MCAR.cov.55.cl.prop.men15.25.F.15.25 <- d.MACR.cov.55.cl.prop.men[,1]
vector.MCAR.cov.55.cl.prop.women15.25.M.15.25 <- d.MACR.cov.55.cl.prop.women[,1]

vector.MCAR.cov.55.cl.prop.men25.40.F.15.25 <- d.MACR.cov.55.cl.prop.men[,2]
vector.MCAR.cov.55.cl.prop.women15.25.M.25.40 <- d.MACR.cov.55.cl.prop.women[,4]

vector.MCAR.cov.55.cl.prop.men25.40.F.25.40 <- d.MACR.cov.55.cl.prop.men[,5]
vector.MCAR.cov.55.cl.prop.women25.40.M.25.40 <- d.MACR.cov.55.cl.prop.women[,5]

vector.MCAR.cov.55.cl.prop.men40.50.F.15.25 <- d.MACR.cov.55.cl.prop.men[,3]
vector.MCAR.cov.55.cl.prop.women15.25.M.40.50 <- d.MACR.cov.55.cl.prop.women[,7]

vector.MCAR.cov.55.cl.prop.men40.50.F.25.40 <- d.MACR.cov.55.cl.prop.men[,6]
vector.MCAR.cov.55.cl.prop.women25.40.M.40.50 <- d.MACR.cov.55.cl.prop.women[,8]


# 60

d.MCAR.cov.60 <- d.MCAR %>%
  select(contains("cov.MCAR.60.")) 
d.MACR.cov.60.cl.prop.men <-  d.MCAR.cov.60 %>%
  select(contains(".cl.prop.men"))
d.MACR.cov.60.cl.prop.women <-  d.MCAR.cov.60 %>%
  select(contains(".cl.prop.women"))


vector.MCAR.cov.60.cl.prop.men15.25.F.15.25 <- d.MACR.cov.60.cl.prop.men[,1]
vector.MCAR.cov.60.cl.prop.women15.25.M.15.25 <- d.MACR.cov.60.cl.prop.women[,1]

vector.MCAR.cov.60.cl.prop.men25.40.F.15.25 <- d.MACR.cov.60.cl.prop.men[,2]
vector.MCAR.cov.60.cl.prop.women15.25.M.25.40 <- d.MACR.cov.60.cl.prop.women[,4]

vector.MCAR.cov.60.cl.prop.men25.40.F.25.40 <- d.MACR.cov.60.cl.prop.men[,5]
vector.MCAR.cov.60.cl.prop.women25.40.M.25.40 <- d.MACR.cov.60.cl.prop.women[,5]

vector.MCAR.cov.60.cl.prop.men40.50.F.15.25 <- d.MACR.cov.60.cl.prop.men[,3]
vector.MCAR.cov.60.cl.prop.women15.25.M.40.50 <- d.MACR.cov.60.cl.prop.women[,7]

vector.MCAR.cov.60.cl.prop.men40.50.F.25.40 <- d.MACR.cov.60.cl.prop.men[,6]
vector.MCAR.cov.60.cl.prop.women25.40.M.40.50 <- d.MACR.cov.60.cl.prop.women[,8]



# 65

d.MCAR.cov.65 <- d.MCAR %>%
  select(contains("cov.MCAR.65.")) 
d.MACR.cov.65.cl.prop.men <-  d.MCAR.cov.65 %>%
  select(contains(".cl.prop.men"))
d.MACR.cov.65.cl.prop.women <-  d.MCAR.cov.65 %>%
  select(contains(".cl.prop.women"))


vector.MCAR.cov.65.cl.prop.men15.25.F.15.25 <- d.MACR.cov.65.cl.prop.men[,1]
vector.MCAR.cov.65.cl.prop.women15.25.M.15.25 <- d.MACR.cov.65.cl.prop.women[,1]

vector.MCAR.cov.65.cl.prop.men25.40.F.15.25 <- d.MACR.cov.65.cl.prop.men[,2]
vector.MCAR.cov.65.cl.prop.women15.25.M.25.40 <- d.MACR.cov.65.cl.prop.women[,4]

vector.MCAR.cov.65.cl.prop.men25.40.F.25.40 <- d.MACR.cov.65.cl.prop.men[,5]
vector.MCAR.cov.65.cl.prop.women25.40.M.25.40 <- d.MACR.cov.65.cl.prop.women[,5]

vector.MCAR.cov.65.cl.prop.men40.50.F.15.25 <- d.MACR.cov.65.cl.prop.men[,3]
vector.MCAR.cov.65.cl.prop.women15.25.M.40.50 <- d.MACR.cov.65.cl.prop.women[,7]

vector.MCAR.cov.65.cl.prop.men40.50.F.25.40 <- d.MACR.cov.65.cl.prop.men[,6]
vector.MCAR.cov.65.cl.prop.women25.40.M.40.50 <- d.MACR.cov.65.cl.prop.women[,8]



# 70

d.MCAR.cov.70 <- d.MCAR %>%
  select(contains("cov.MCAR.70.")) 
d.MACR.cov.70.cl.prop.men <-  d.MCAR.cov.70 %>%
  select(contains(".cl.prop.men"))
d.MACR.cov.70.cl.prop.women <-  d.MCAR.cov.70 %>%
  select(contains(".cl.prop.women"))


vector.MCAR.cov.70.cl.prop.men15.25.F.15.25 <- d.MACR.cov.70.cl.prop.men[,1]
vector.MCAR.cov.70.cl.prop.women15.25.M.15.25 <- d.MACR.cov.70.cl.prop.women[,1]

vector.MCAR.cov.70.cl.prop.men25.40.F.15.25 <- d.MACR.cov.70.cl.prop.men[,2]
vector.MCAR.cov.70.cl.prop.women15.25.M.25.40 <- d.MACR.cov.70.cl.prop.women[,4]

vector.MCAR.cov.70.cl.prop.men25.40.F.25.40 <- d.MACR.cov.70.cl.prop.men[,5]
vector.MCAR.cov.70.cl.prop.women25.40.M.25.40 <- d.MACR.cov.70.cl.prop.women[,5]

vector.MCAR.cov.70.cl.prop.men40.50.F.15.25 <- d.MACR.cov.70.cl.prop.men[,3]
vector.MCAR.cov.70.cl.prop.women15.25.M.40.50 <- d.MACR.cov.70.cl.prop.women[,7]

vector.MCAR.cov.70.cl.prop.men40.50.F.25.40 <- d.MACR.cov.70.cl.prop.men[,6]
vector.MCAR.cov.70.cl.prop.women25.40.M.40.50 <- d.MACR.cov.70.cl.prop.women[,8]



# 75

d.MCAR.cov.75 <- d.MCAR %>%
  select(contains("cov.MCAR.75.")) 
d.MACR.cov.75.cl.prop.men <-  d.MCAR.cov.75 %>%
  select(contains(".cl.prop.men"))
d.MACR.cov.75.cl.prop.women <-  d.MCAR.cov.75 %>%
  select(contains(".cl.prop.women"))


vector.MCAR.cov.75.cl.prop.men15.25.F.15.25 <- d.MACR.cov.75.cl.prop.men[,1]
vector.MCAR.cov.75.cl.prop.women15.25.M.15.25 <- d.MACR.cov.75.cl.prop.women[,1]

vector.MCAR.cov.75.cl.prop.men25.40.F.15.25 <- d.MACR.cov.75.cl.prop.men[,2]
vector.MCAR.cov.75.cl.prop.women15.25.M.25.40 <- d.MACR.cov.75.cl.prop.women[,4]

vector.MCAR.cov.75.cl.prop.men25.40.F.25.40 <- d.MACR.cov.75.cl.prop.men[,5]
vector.MCAR.cov.75.cl.prop.women25.40.M.25.40 <- d.MACR.cov.75.cl.prop.women[,5]

vector.MCAR.cov.75.cl.prop.men40.50.F.15.25 <- d.MACR.cov.75.cl.prop.men[,3]
vector.MCAR.cov.75.cl.prop.women15.25.M.40.50 <- d.MACR.cov.75.cl.prop.women[,7]

vector.MCAR.cov.75.cl.prop.men40.50.F.25.40 <- d.MACR.cov.75.cl.prop.men[,6]
vector.MCAR.cov.75.cl.prop.women25.40.M.40.50 <- d.MACR.cov.75.cl.prop.women[,8]



# 80

d.MCAR.cov.80 <- d.MCAR %>%
  select(contains("cov.MCAR.80.")) 
d.MACR.cov.80.cl.prop.men <-  d.MCAR.cov.80 %>%
  select(contains(".cl.prop.men"))
d.MACR.cov.80.cl.prop.women <-  d.MCAR.cov.80 %>%
  select(contains(".cl.prop.women"))


vector.MCAR.cov.80.cl.prop.men15.25.F.15.25 <- d.MACR.cov.80.cl.prop.men[,1]
vector.MCAR.cov.80.cl.prop.women15.25.M.15.25 <- d.MACR.cov.80.cl.prop.women[,1]

vector.MCAR.cov.80.cl.prop.men25.40.F.15.25 <- d.MACR.cov.80.cl.prop.men[,2]
vector.MCAR.cov.80.cl.prop.women15.25.M.25.40 <- d.MACR.cov.80.cl.prop.women[,4]

vector.MCAR.cov.80.cl.prop.men25.40.F.25.40 <- d.MACR.cov.80.cl.prop.men[,5]
vector.MCAR.cov.80.cl.prop.women25.40.M.25.40 <- d.MACR.cov.80.cl.prop.women[,5]

vector.MCAR.cov.80.cl.prop.men40.50.F.15.25 <- d.MACR.cov.80.cl.prop.men[,3]
vector.MCAR.cov.80.cl.prop.women15.25.M.40.50 <- d.MACR.cov.80.cl.prop.women[,7]

vector.MCAR.cov.80.cl.prop.men40.50.F.25.40 <- d.MACR.cov.80.cl.prop.men[,6]
vector.MCAR.cov.80.cl.prop.women25.40.M.40.50 <- d.MACR.cov.80.cl.prop.women[,8]



# 85

d.MCAR.cov.85 <- d.MCAR %>%
  select(contains("cov.MCAR.85.")) 
d.MACR.cov.85.cl.prop.men <-  d.MCAR.cov.85 %>%
  select(contains(".cl.prop.men"))
d.MACR.cov.85.cl.prop.women <-  d.MCAR.cov.85 %>%
  select(contains(".cl.prop.women"))


vector.MCAR.cov.85.cl.prop.men15.25.F.15.25 <- d.MACR.cov.85.cl.prop.men[,1]
vector.MCAR.cov.85.cl.prop.women15.25.M.15.25 <- d.MACR.cov.85.cl.prop.women[,1]

vector.MCAR.cov.85.cl.prop.men25.40.F.15.25 <- d.MACR.cov.85.cl.prop.men[,2]
vector.MCAR.cov.85.cl.prop.women15.25.M.25.40 <- d.MACR.cov.85.cl.prop.women[,4]

vector.MCAR.cov.85.cl.prop.men25.40.F.25.40 <- d.MACR.cov.85.cl.prop.men[,5]
vector.MCAR.cov.85.cl.prop.women25.40.M.25.40 <- d.MACR.cov.85.cl.prop.women[,5]

vector.MCAR.cov.85.cl.prop.men40.50.F.15.25 <- d.MACR.cov.85.cl.prop.men[,3]
vector.MCAR.cov.85.cl.prop.women15.25.M.40.50 <- d.MACR.cov.85.cl.prop.women[,7]

vector.MCAR.cov.85.cl.prop.men40.50.F.25.40 <- d.MACR.cov.85.cl.prop.men[,6]
vector.MCAR.cov.85.cl.prop.women25.40.M.40.50 <- d.MACR.cov.85.cl.prop.women[,8]


# 90

d.MCAR.cov.90 <- d.MCAR %>%
  select(contains("cov.MCAR.90.")) 
d.MACR.cov.90.cl.prop.men <-  d.MCAR.cov.90 %>%
  select(contains(".cl.prop.men"))
d.MACR.cov.90.cl.prop.women <-  d.MCAR.cov.90 %>%
  select(contains(".cl.prop.women"))


vector.MCAR.cov.90.cl.prop.men15.25.F.15.25 <- d.MACR.cov.90.cl.prop.men[,1]
vector.MCAR.cov.90.cl.prop.women15.25.M.15.25 <- d.MACR.cov.90.cl.prop.women[,1]

vector.MCAR.cov.90.cl.prop.men25.40.F.15.25 <- d.MACR.cov.90.cl.prop.men[,2]
vector.MCAR.cov.90.cl.prop.women15.25.M.25.40 <- d.MACR.cov.90.cl.prop.women[,4]

vector.MCAR.cov.90.cl.prop.men25.40.F.25.40 <- d.MACR.cov.90.cl.prop.men[,5]
vector.MCAR.cov.90.cl.prop.women25.40.M.25.40 <- d.MACR.cov.90.cl.prop.women[,5]

vector.MCAR.cov.90.cl.prop.men40.50.F.15.25 <- d.MACR.cov.90.cl.prop.men[,3]
vector.MCAR.cov.90.cl.prop.women15.25.M.40.50 <- d.MACR.cov.90.cl.prop.women[,7]

vector.MCAR.cov.90.cl.prop.men40.50.F.25.40 <- d.MACR.cov.90.cl.prop.men[,6]
vector.MCAR.cov.90.cl.prop.women25.40.M.40.50 <- d.MACR.cov.90.cl.prop.women[,8]


# 95

d.MCAR.cov.95 <- d.MCAR %>%
  select(contains("cov.MCAR.95.")) 
d.MACR.cov.95.cl.prop.men <-  d.MCAR.cov.95 %>%
  select(contains(".cl.prop.men"))
d.MACR.cov.95.cl.prop.women <-  d.MCAR.cov.95 %>%
  select(contains(".cl.prop.women"))


vector.MCAR.cov.95.cl.prop.men15.25.F.15.25 <- d.MACR.cov.95.cl.prop.men[,1]
vector.MCAR.cov.95.cl.prop.women15.25.M.15.25 <- d.MACR.cov.95.cl.prop.women[,1]

vector.MCAR.cov.95.cl.prop.men25.40.F.15.25 <- d.MACR.cov.95.cl.prop.men[,2]
vector.MCAR.cov.95.cl.prop.women15.25.M.25.40 <- d.MACR.cov.95.cl.prop.women[,4]

vector.MCAR.cov.95.cl.prop.men25.40.F.25.40 <- d.MACR.cov.95.cl.prop.men[,5]
vector.MCAR.cov.95.cl.prop.women25.40.M.25.40 <- d.MACR.cov.95.cl.prop.women[,5]

vector.MCAR.cov.95.cl.prop.men40.50.F.15.25 <- d.MACR.cov.95.cl.prop.men[,3]
vector.MCAR.cov.95.cl.prop.women15.25.M.40.50 <- d.MACR.cov.95.cl.prop.women[,7]

vector.MCAR.cov.95.cl.prop.men40.50.F.25.40 <- d.MACR.cov.95.cl.prop.men[,6]
vector.MCAR.cov.95.cl.prop.women25.40.M.40.50 <- d.MACR.cov.95.cl.prop.women[,8]



# MAR_

d.MAR <- dr %>%
  select(contains("MAR."))



# MAR - a

# 35

d.MAR.cov.35 <- d.MAR %>%
  select(contains("cov.MAR.a.35.")) 
d.MAR.cov.35.cl.prop.men <-  d.MAR.cov.35 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.35.cl.prop.women <-  d.MAR.cov.35 %>%
  select(contains(".cl.prop.women"))

vector.MAR.a.cov.35.cl.prop.men15.25.F.15.25 <- d.MAR.cov.35.cl.prop.men[,1]
vector.MAR.a.cov.35.cl.prop.women15.25.M.15.25 <- d.MAR.cov.35.cl.prop.women[,1]

vector.MAR.a.cov.35.cl.prop.men25.40.F.15.25 <- d.MAR.cov.35.cl.prop.men[,2]
vector.MAR.a.cov.35.cl.prop.women15.25.M.25.40 <- d.MAR.cov.35.cl.prop.women[,4]

vector.MAR.a.cov.35.cl.prop.men25.40.F.25.40 <- d.MAR.cov.35.cl.prop.men[,5]
vector.MAR.a.cov.35.cl.prop.women25.40.M.25.40 <- d.MAR.cov.35.cl.prop.women[,5]

vector.MAR.a.cov.35.cl.prop.men40.50.F.15.25 <- d.MAR.cov.35.cl.prop.men[,3]
vector.MAR.a.cov.35.cl.prop.women15.25.M.40.50 <- d.MAR.cov.35.cl.prop.women[,7]

vector.MAR.a.cov.35.cl.prop.men40.50.F.25.40 <- d.MAR.cov.35.cl.prop.men[,6]
vector.MAR.a.cov.35.cl.prop.women25.40.M.40.50 <- d.MAR.cov.35.cl.prop.women[,8]




# 40

d.MAR.cov.40 <- d.MAR %>%
  select(contains("cov.MAR.a.40.")) 
d.MAR.cov.40.cl.prop.men <-  d.MAR.cov.40 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.40.cl.prop.women <-  d.MAR.cov.40 %>%
  select(contains(".cl.prop.women"))

vector.MAR.a.cov.40.cl.prop.men15.25.F.15.25 <- d.MAR.cov.40.cl.prop.men[,1]
vector.MAR.a.cov.40.cl.prop.women15.25.M.15.25 <- d.MAR.cov.40.cl.prop.women[,1]

vector.MAR.a.cov.40.cl.prop.men25.40.F.15.25 <- d.MAR.cov.40.cl.prop.men[,2]
vector.MAR.a.cov.40.cl.prop.women15.25.M.25.40 <- d.MAR.cov.40.cl.prop.women[,4]

vector.MAR.a.cov.40.cl.prop.men25.40.F.25.40 <- d.MAR.cov.40.cl.prop.men[,5]
vector.MAR.a.cov.40.cl.prop.women25.40.M.25.40 <- d.MAR.cov.40.cl.prop.women[,5]

vector.MAR.a.cov.40.cl.prop.men40.50.F.15.25 <- d.MAR.cov.40.cl.prop.men[,3]
vector.MAR.a.cov.40.cl.prop.women15.25.M.40.50 <- d.MAR.cov.40.cl.prop.women[,7]

vector.MAR.a.cov.40.cl.prop.men40.50.F.25.40 <- d.MAR.cov.40.cl.prop.men[,6]
vector.MAR.a.cov.40.cl.prop.women25.40.M.40.50 <- d.MAR.cov.40.cl.prop.women[,8]


# 45

d.MAR.cov.45 <- d.MAR %>%
  select(contains("cov.MAR.a.45.")) 
d.MAR.cov.45.cl.prop.men <-  d.MAR.cov.45 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.45.cl.prop.women <-  d.MAR.cov.45 %>%
  select(contains(".cl.prop.women"))

vector.MAR.a.cov.45.cl.prop.men15.25.F.15.25 <- d.MAR.cov.45.cl.prop.men[,1]
vector.MAR.a.cov.45.cl.prop.women15.25.M.15.25 <- d.MAR.cov.45.cl.prop.women[,1]

vector.MAR.a.cov.45.cl.prop.men25.40.F.15.25 <- d.MAR.cov.45.cl.prop.men[,2]
vector.MAR.a.cov.45.cl.prop.women15.25.M.25.40 <- d.MAR.cov.45.cl.prop.women[,4]

vector.MAR.a.cov.45.cl.prop.men25.40.F.25.40 <- d.MAR.cov.45.cl.prop.men[,5]
vector.MAR.a.cov.45.cl.prop.women25.40.M.25.40 <- d.MAR.cov.45.cl.prop.women[,5]

vector.MAR.a.cov.45.cl.prop.men40.50.F.15.25 <- d.MAR.cov.45.cl.prop.men[,3]
vector.MAR.a.cov.45.cl.prop.women15.25.M.40.50 <- d.MAR.cov.45.cl.prop.women[,7]

vector.MAR.a.cov.45.cl.prop.men40.50.F.25.40 <- d.MAR.cov.45.cl.prop.men[,6]
vector.MAR.a.cov.45.cl.prop.women25.40.M.40.50 <- d.MAR.cov.45.cl.prop.women[,8]


# 50

d.MAR.cov.50 <- d.MAR %>%
  select(contains("cov.MAR.a.50.")) 
d.MAR.cov.50.cl.prop.men <-  d.MAR.cov.50 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.50.cl.prop.women <-  d.MAR.cov.50 %>%
  select(contains(".cl.prop.women"))


vector.MAR.a.cov.50.cl.prop.men15.25.F.15.25 <- d.MAR.cov.50.cl.prop.men[,1]
vector.MAR.a.cov.50.cl.prop.women15.25.M.15.25 <- d.MAR.cov.50.cl.prop.women[,1]

vector.MAR.a.cov.50.cl.prop.men25.40.F.15.25 <- d.MAR.cov.50.cl.prop.men[,2]
vector.MAR.a.cov.50.cl.prop.women15.25.M.25.40 <- d.MAR.cov.50.cl.prop.women[,4]

vector.MAR.a.cov.50.cl.prop.men25.40.F.25.40 <- d.MAR.cov.50.cl.prop.men[,5]
vector.MAR.a.cov.50.cl.prop.women25.40.M.25.40 <- d.MAR.cov.50.cl.prop.women[,5]

vector.MAR.a.cov.50.cl.prop.men40.50.F.15.25 <- d.MAR.cov.50.cl.prop.men[,3]
vector.MAR.a.cov.50.cl.prop.women15.25.M.40.50 <- d.MAR.cov.50.cl.prop.women[,7]

vector.MAR.a.cov.50.cl.prop.men40.50.F.25.40 <- d.MAR.cov.50.cl.prop.men[,6]
vector.MAR.a.cov.50.cl.prop.women25.40.M.40.50 <- d.MAR.cov.50.cl.prop.women[,8]


# 55

d.MAR.cov.55 <- d.MAR %>%
  select(contains("cov.MAR.a.55.")) 
d.MAR.cov.55.cl.prop.men <-  d.MAR.cov.55 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.55.cl.prop.women <-  d.MAR.cov.55 %>%
  select(contains(".cl.prop.women"))

vector.MAR.a.cov.55.cl.prop.men15.25.F.15.25 <- d.MAR.cov.55.cl.prop.men[,1]
vector.MAR.a.cov.55.cl.prop.women15.25.M.15.25 <- d.MAR.cov.55.cl.prop.women[,1]

vector.MAR.a.cov.55.cl.prop.men25.40.F.15.25 <- d.MAR.cov.55.cl.prop.men[,2]
vector.MAR.a.cov.55.cl.prop.women15.25.M.25.40 <- d.MAR.cov.55.cl.prop.women[,4]

vector.MAR.a.cov.55.cl.prop.men25.40.F.25.40 <- d.MAR.cov.55.cl.prop.men[,5]
vector.MAR.a.cov.55.cl.prop.women25.40.M.25.40 <- d.MAR.cov.55.cl.prop.women[,5]

vector.MAR.a.cov.55.cl.prop.men40.50.F.15.25 <- d.MAR.cov.55.cl.prop.men[,3]
vector.MAR.a.cov.55.cl.prop.women15.25.M.40.50 <- d.MAR.cov.55.cl.prop.women[,7]

vector.MAR.a.cov.55.cl.prop.men40.50.F.25.40 <- d.MAR.cov.55.cl.prop.men[,6]
vector.MAR.a.cov.55.cl.prop.women25.40.M.40.50 <- d.MAR.cov.55.cl.prop.women[,8]


# 60

d.MAR.cov.60 <- d.MAR %>%
  select(contains("cov.MAR.a.60.")) 
d.MAR.cov.60.cl.prop.men <-  d.MAR.cov.60 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.60.cl.prop.women <-  d.MAR.cov.60 %>%
  select(contains(".cl.prop.women"))

vector.MAR.a.cov.60.cl.prop.men15.25.F.15.25 <- d.MAR.cov.60.cl.prop.men[,1]
vector.MAR.a.cov.60.cl.prop.women15.25.M.15.25 <- d.MAR.cov.60.cl.prop.women[,1]

vector.MAR.a.cov.60.cl.prop.men25.40.F.15.25 <- d.MAR.cov.60.cl.prop.men[,2]
vector.MAR.a.cov.60.cl.prop.women15.25.M.25.40 <- d.MAR.cov.60.cl.prop.women[,4]

vector.MAR.a.cov.60.cl.prop.men25.40.F.25.40 <- d.MAR.cov.60.cl.prop.men[,5]
vector.MAR.a.cov.60.cl.prop.women25.40.M.25.40 <- d.MAR.cov.60.cl.prop.women[,5]

vector.MAR.a.cov.60.cl.prop.men40.50.F.15.25 <- d.MAR.cov.60.cl.prop.men[,3]
vector.MAR.a.cov.60.cl.prop.women15.25.M.40.50 <- d.MAR.cov.60.cl.prop.women[,7]

vector.MAR.a.cov.60.cl.prop.men40.50.F.25.40 <- d.MAR.cov.60.cl.prop.men[,6]
vector.MAR.a.cov.60.cl.prop.women25.40.M.40.50 <- d.MAR.cov.60.cl.prop.women[,8]


# 65

d.MAR.cov.65 <- d.MAR %>%
  select(contains("cov.MAR.a.65.")) 
d.MAR.cov.65.cl.prop.men <-  d.MAR.cov.65 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.65.cl.prop.women <-  d.MAR.cov.65 %>%
  select(contains(".cl.prop.women"))

vector.MAR.a.cov.65.cl.prop.men15.25.F.15.25 <- d.MAR.cov.65.cl.prop.men[,1]
vector.MAR.a.cov.65.cl.prop.women15.25.M.15.25 <- d.MAR.cov.65.cl.prop.women[,1]

vector.MAR.a.cov.65.cl.prop.men25.40.F.15.25 <- d.MAR.cov.65.cl.prop.men[,2]
vector.MAR.a.cov.65.cl.prop.women15.25.M.25.40 <- d.MAR.cov.65.cl.prop.women[,4]

vector.MAR.a.cov.65.cl.prop.men25.40.F.25.40 <- d.MAR.cov.65.cl.prop.men[,5]
vector.MAR.a.cov.65.cl.prop.women25.40.M.25.40 <- d.MAR.cov.65.cl.prop.women[,5]

vector.MAR.a.cov.65.cl.prop.men40.50.F.15.25 <- d.MAR.cov.65.cl.prop.men[,3]
vector.MAR.a.cov.65.cl.prop.women15.25.M.40.50 <- d.MAR.cov.65.cl.prop.women[,7]

vector.MAR.a.cov.65.cl.prop.men40.50.F.25.40 <- d.MAR.cov.65.cl.prop.men[,6]
vector.MAR.a.cov.65.cl.prop.women25.40.M.40.50 <- d.MAR.cov.65.cl.prop.women[,8]


# 70

d.MAR.cov.70 <- d.MAR %>%
  select(contains("cov.MAR.a.70.")) 
d.MAR.cov.70.cl.prop.men <-  d.MAR.cov.70 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.70.cl.prop.women <-  d.MAR.cov.70 %>%
  select(contains(".cl.prop.women"))


vector.MAR.a.cov.70.cl.prop.men15.25.F.15.25 <- d.MAR.cov.70.cl.prop.men[,1]
vector.MAR.a.cov.70.cl.prop.women15.25.M.15.25 <- d.MAR.cov.70.cl.prop.women[,1]

vector.MAR.a.cov.70.cl.prop.men25.40.F.15.25 <- d.MAR.cov.70.cl.prop.men[,2]
vector.MAR.a.cov.70.cl.prop.women15.25.M.25.40 <- d.MAR.cov.70.cl.prop.women[,4]

vector.MAR.a.cov.70.cl.prop.men25.40.F.25.40 <- d.MAR.cov.70.cl.prop.men[,5]
vector.MAR.a.cov.70.cl.prop.women25.40.M.25.40 <- d.MAR.cov.70.cl.prop.women[,5]

vector.MAR.a.cov.70.cl.prop.men40.50.F.15.25 <- d.MAR.cov.70.cl.prop.men[,3]
vector.MAR.a.cov.70.cl.prop.women15.25.M.40.50 <- d.MAR.cov.70.cl.prop.women[,7]

vector.MAR.a.cov.70.cl.prop.men40.50.F.25.40 <- d.MAR.cov.70.cl.prop.men[,6]
vector.MAR.a.cov.70.cl.prop.women25.40.M.40.50 <- d.MAR.cov.70.cl.prop.women[,8]


# 75

d.MAR.cov.75 <- d.MAR %>%
  select(contains("cov.MAR.a.75.")) 
d.MAR.cov.75.cl.prop.men <-  d.MAR.cov.75 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.75.cl.prop.women <-  d.MAR.cov.75 %>%
  select(contains(".cl.prop.women"))



vector.MAR.a.cov.75.cl.prop.men15.25.F.15.25 <- d.MAR.cov.75.cl.prop.men[,1]
vector.MAR.a.cov.75.cl.prop.women15.25.M.15.25 <- d.MAR.cov.75.cl.prop.women[,1]

vector.MAR.a.cov.75.cl.prop.men25.40.F.15.25 <- d.MAR.cov.75.cl.prop.men[,2]
vector.MAR.a.cov.75.cl.prop.women15.25.M.25.40 <- d.MAR.cov.75.cl.prop.women[,4]

vector.MAR.a.cov.75.cl.prop.men25.40.F.25.40 <- d.MAR.cov.75.cl.prop.men[,5]
vector.MAR.a.cov.75.cl.prop.women25.40.M.25.40 <- d.MAR.cov.75.cl.prop.women[,5]

vector.MAR.a.cov.75.cl.prop.men40.50.F.15.25 <- d.MAR.cov.75.cl.prop.men[,3]
vector.MAR.a.cov.75.cl.prop.women15.25.M.40.50 <- d.MAR.cov.75.cl.prop.women[,7]

vector.MAR.a.cov.75.cl.prop.men40.50.F.25.40 <- d.MAR.cov.75.cl.prop.men[,6]
vector.MAR.a.cov.75.cl.prop.women25.40.M.40.50 <- d.MAR.cov.75.cl.prop.women[,8]


# 80

d.MAR.cov.80 <- d.MAR %>%
  select(contains("cov.MAR.a.80.")) 
d.MAR.cov.80.cl.prop.men <-  d.MAR.cov.80 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.80.cl.prop.women <-  d.MAR.cov.80 %>%
  select(contains(".cl.prop.women"))


vector.MAR.a.cov.80.cl.prop.men15.25.F.15.25 <- d.MAR.cov.80.cl.prop.men[,1]
vector.MAR.a.cov.80.cl.prop.women15.25.M.15.25 <- d.MAR.cov.80.cl.prop.women[,1]

vector.MAR.a.cov.80.cl.prop.men25.40.F.15.25 <- d.MAR.cov.80.cl.prop.men[,2]
vector.MAR.a.cov.80.cl.prop.women15.25.M.25.40 <- d.MAR.cov.80.cl.prop.women[,4]

vector.MAR.a.cov.80.cl.prop.men25.40.F.25.40 <- d.MAR.cov.80.cl.prop.men[,5]
vector.MAR.a.cov.80.cl.prop.women25.40.M.25.40 <- d.MAR.cov.80.cl.prop.women[,5]

vector.MAR.a.cov.80.cl.prop.men40.50.F.15.25 <- d.MAR.cov.80.cl.prop.men[,3]
vector.MAR.a.cov.80.cl.prop.women15.25.M.40.50 <- d.MAR.cov.80.cl.prop.women[,7]

vector.MAR.a.cov.80.cl.prop.men40.50.F.25.40 <- d.MAR.cov.80.cl.prop.men[,6]
vector.MAR.a.cov.80.cl.prop.women25.40.M.40.50 <- d.MAR.cov.80.cl.prop.women[,8]


# 85

d.MAR.cov.85 <- d.MAR %>%
  select(contains("cov.MAR.a.85.")) 
d.MAR.cov.85.cl.prop.men <-  d.MAR.cov.85 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.85.cl.prop.women <-  d.MAR.cov.85 %>%
  select(contains(".cl.prop.women"))


vector.MAR.a.cov.85.cl.prop.men15.25.F.15.25 <- d.MAR.cov.85.cl.prop.men[,1]
vector.MAR.a.cov.85.cl.prop.women15.25.M.15.25 <- d.MAR.cov.85.cl.prop.women[,1]

vector.MAR.a.cov.85.cl.prop.men25.40.F.15.25 <- d.MAR.cov.85.cl.prop.men[,2]
vector.MAR.a.cov.85.cl.prop.women15.25.M.25.40 <- d.MAR.cov.85.cl.prop.women[,4]

vector.MAR.a.cov.85.cl.prop.men25.40.F.25.40 <- d.MAR.cov.85.cl.prop.men[,5]
vector.MAR.a.cov.85.cl.prop.women25.40.M.25.40 <- d.MAR.cov.85.cl.prop.women[,5]

vector.MAR.a.cov.85.cl.prop.men40.50.F.15.25 <- d.MAR.cov.85.cl.prop.men[,3]
vector.MAR.a.cov.85.cl.prop.women15.25.M.40.50 <- d.MAR.cov.85.cl.prop.women[,7]

vector.MAR.a.cov.85.cl.prop.men40.50.F.25.40 <- d.MAR.cov.85.cl.prop.men[,6]
vector.MAR.a.cov.85.cl.prop.women25.40.M.40.50 <- d.MAR.cov.85.cl.prop.women[,8]


# 90

d.MAR.cov.90 <- d.MAR %>%
  select(contains("cov.MAR.a.90.")) 
d.MAR.cov.90.cl.prop.men <-  d.MAR.cov.90 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.90.cl.prop.women <-  d.MAR.cov.90 %>%
  select(contains(".cl.prop.women"))


vector.MAR.a.cov.90.cl.prop.men15.25.F.15.25 <- d.MAR.cov.90.cl.prop.men[,1]
vector.MAR.a.cov.90.cl.prop.women15.25.M.15.25 <- d.MAR.cov.90.cl.prop.women[,1]

vector.MAR.a.cov.90.cl.prop.men25.40.F.15.25 <- d.MAR.cov.90.cl.prop.men[,2]
vector.MAR.a.cov.90.cl.prop.women15.25.M.25.40 <- d.MAR.cov.90.cl.prop.women[,4]

vector.MAR.a.cov.90.cl.prop.men25.40.F.25.40 <- d.MAR.cov.90.cl.prop.men[,5]
vector.MAR.a.cov.90.cl.prop.women25.40.M.25.40 <- d.MAR.cov.90.cl.prop.women[,5]

vector.MAR.a.cov.90.cl.prop.men40.50.F.15.25 <- d.MAR.cov.90.cl.prop.men[,3]
vector.MAR.a.cov.90.cl.prop.women15.25.M.40.50 <- d.MAR.cov.90.cl.prop.women[,7]

vector.MAR.a.cov.90.cl.prop.men40.50.F.25.40 <- d.MAR.cov.90.cl.prop.men[,6]
vector.MAR.a.cov.90.cl.prop.women25.40.M.40.50 <- d.MAR.cov.90.cl.prop.women[,8]


# 95

d.MAR.cov.95 <- d.MAR %>%
  select(contains("cov.MAR.a.95.")) 
d.MAR.cov.95.cl.prop.men <-  d.MAR.cov.95 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.95.cl.prop.women <-  d.MAR.cov.95 %>%
  select(contains(".cl.prop.women"))


vector.MAR.a.cov.95.cl.prop.men15.25.F.15.25 <- d.MAR.cov.95.cl.prop.men[,1]
vector.MAR.a.cov.95.cl.prop.women15.25.M.15.25 <- d.MAR.cov.95.cl.prop.women[,1]

vector.MAR.a.cov.95.cl.prop.men25.40.F.15.25 <- d.MAR.cov.95.cl.prop.men[,2]
vector.MAR.a.cov.95.cl.prop.women15.25.M.25.40 <- d.MAR.cov.95.cl.prop.women[,4]

vector.MAR.a.cov.95.cl.prop.men25.40.F.25.40 <- d.MAR.cov.95.cl.prop.men[,5]
vector.MAR.a.cov.95.cl.prop.women25.40.M.25.40 <- d.MAR.cov.95.cl.prop.women[,5]

vector.MAR.a.cov.95.cl.prop.men40.50.F.15.25 <- d.MAR.cov.95.cl.prop.men[,3]
vector.MAR.a.cov.95.cl.prop.women15.25.M.40.50 <- d.MAR.cov.95.cl.prop.women[,7]

vector.MAR.a.cov.95.cl.prop.men40.50.F.25.40 <- d.MAR.cov.95.cl.prop.men[,6]
vector.MAR.a.cov.95.cl.prop.women25.40.M.40.50 <- d.MAR.cov.95.cl.prop.women[,8]



# MAR - b


# 35

d.MAR.cov.35 <- d.MAR %>%
  select(contains("cov.MAR.b.35.")) 
d.MAR.cov.35.cl.prop.men <-  d.MAR.cov.35 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.35.cl.prop.women <-  d.MAR.cov.35 %>%
  select(contains(".cl.prop.women"))

vector.MAR.b.cov.35.cl.prop.men15.25.F.15.25 <- d.MAR.cov.35.cl.prop.men[,1]
vector.MAR.b.cov.35.cl.prop.women15.25.M.15.25 <- d.MAR.cov.35.cl.prop.women[,1]

vector.MAR.b.cov.35.cl.prop.men25.40.F.15.25 <- d.MAR.cov.35.cl.prop.men[,2]
vector.MAR.b.cov.35.cl.prop.women15.25.M.25.40 <- d.MAR.cov.35.cl.prop.women[,4]

vector.MAR.b.cov.35.cl.prop.men25.40.F.25.40 <- d.MAR.cov.35.cl.prop.men[,5]
vector.MAR.b.cov.35.cl.prop.women25.40.M.25.40 <- d.MAR.cov.35.cl.prop.women[,5]

vector.MAR.b.cov.35.cl.prop.men40.50.F.15.25 <- d.MAR.cov.35.cl.prop.men[,3]
vector.MAR.b.cov.35.cl.prop.women15.25.M.40.50 <- d.MAR.cov.35.cl.prop.women[,7]

vector.MAR.b.cov.35.cl.prop.men40.50.F.25.40 <- d.MAR.cov.35.cl.prop.men[,6]
vector.MAR.b.cov.35.cl.prop.women25.40.M.40.50 <- d.MAR.cov.35.cl.prop.women[,8]




# 40

d.MAR.cov.40 <- d.MAR %>%
  select(contains("cov.MAR.b.40.")) 
d.MAR.cov.40.cl.prop.men <-  d.MAR.cov.40 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.40.cl.prop.women <-  d.MAR.cov.40 %>%
  select(contains(".cl.prop.women"))

vector.MAR.b.cov.40.cl.prop.men15.25.F.15.25 <- d.MAR.cov.40.cl.prop.men[,1]
vector.MAR.b.cov.40.cl.prop.women15.25.M.15.25 <- d.MAR.cov.40.cl.prop.women[,1]

vector.MAR.b.cov.40.cl.prop.men25.40.F.15.25 <- d.MAR.cov.40.cl.prop.men[,2]
vector.MAR.b.cov.40.cl.prop.women15.25.M.25.40 <- d.MAR.cov.40.cl.prop.women[,4]

vector.MAR.b.cov.40.cl.prop.men25.40.F.25.40 <- d.MAR.cov.40.cl.prop.men[,5]
vector.MAR.b.cov.40.cl.prop.women25.40.M.25.40 <- d.MAR.cov.40.cl.prop.women[,5]

vector.MAR.b.cov.40.cl.prop.men40.50.F.15.25 <- d.MAR.cov.40.cl.prop.men[,3]
vector.MAR.b.cov.40.cl.prop.women15.25.M.40.50 <- d.MAR.cov.40.cl.prop.women[,7]

vector.MAR.b.cov.40.cl.prop.men40.50.F.25.40 <- d.MAR.cov.40.cl.prop.men[,6]
vector.MAR.b.cov.40.cl.prop.women25.40.M.40.50 <- d.MAR.cov.40.cl.prop.women[,8]


# 45

d.MAR.cov.45 <- d.MAR %>%
  select(contains("cov.MAR.b.45.")) 
d.MAR.cov.45.cl.prop.men <-  d.MAR.cov.45 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.45.cl.prop.women <-  d.MAR.cov.45 %>%
  select(contains(".cl.prop.women"))

vector.MAR.b.cov.45.cl.prop.men15.25.F.15.25 <- d.MAR.cov.45.cl.prop.men[,1]
vector.MAR.b.cov.45.cl.prop.women15.25.M.15.25 <- d.MAR.cov.45.cl.prop.women[,1]

vector.MAR.b.cov.45.cl.prop.men25.40.F.15.25 <- d.MAR.cov.45.cl.prop.men[,2]
vector.MAR.b.cov.45.cl.prop.women15.25.M.25.40 <- d.MAR.cov.45.cl.prop.women[,4]

vector.MAR.b.cov.45.cl.prop.men25.40.F.25.40 <- d.MAR.cov.45.cl.prop.men[,5]
vector.MAR.b.cov.45.cl.prop.women25.40.M.25.40 <- d.MAR.cov.45.cl.prop.women[,5]

vector.MAR.b.cov.45.cl.prop.men40.50.F.15.25 <- d.MAR.cov.45.cl.prop.men[,3]
vector.MAR.b.cov.45.cl.prop.women15.25.M.40.50 <- d.MAR.cov.45.cl.prop.women[,7]

vector.MAR.b.cov.45.cl.prop.men40.50.F.25.40 <- d.MAR.cov.45.cl.prop.men[,6]
vector.MAR.b.cov.45.cl.prop.women25.40.M.40.50 <- d.MAR.cov.45.cl.prop.women[,8]


# 50

d.MAR.cov.50 <- d.MAR %>%
  select(contains("cov.MAR.b.50.")) 
d.MAR.cov.50.cl.prop.men <-  d.MAR.cov.50 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.50.cl.prop.women <-  d.MAR.cov.50 %>%
  select(contains(".cl.prop.women"))


vector.MAR.b.cov.50.cl.prop.men15.25.F.15.25 <- d.MAR.cov.50.cl.prop.men[,1]
vector.MAR.b.cov.50.cl.prop.women15.25.M.15.25 <- d.MAR.cov.50.cl.prop.women[,1]

vector.MAR.b.cov.50.cl.prop.men25.40.F.15.25 <- d.MAR.cov.50.cl.prop.men[,2]
vector.MAR.b.cov.50.cl.prop.women15.25.M.25.40 <- d.MAR.cov.50.cl.prop.women[,4]

vector.MAR.b.cov.50.cl.prop.men25.40.F.25.40 <- d.MAR.cov.50.cl.prop.men[,5]
vector.MAR.b.cov.50.cl.prop.women25.40.M.25.40 <- d.MAR.cov.50.cl.prop.women[,5]

vector.MAR.b.cov.50.cl.prop.men40.50.F.15.25 <- d.MAR.cov.50.cl.prop.men[,3]
vector.MAR.b.cov.50.cl.prop.women15.25.M.40.50 <- d.MAR.cov.50.cl.prop.women[,7]

vector.MAR.b.cov.50.cl.prop.men40.50.F.25.40 <- d.MAR.cov.50.cl.prop.men[,6]
vector.MAR.b.cov.50.cl.prop.women25.40.M.40.50 <- d.MAR.cov.50.cl.prop.women[,8]


# 55

d.MAR.cov.55 <- d.MAR %>%
  select(contains("cov.MAR.b.55.")) 
d.MAR.cov.55.cl.prop.men <-  d.MAR.cov.55 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.55.cl.prop.women <-  d.MAR.cov.55 %>%
  select(contains(".cl.prop.women"))

vector.MAR.b.cov.55.cl.prop.men15.25.F.15.25 <- d.MAR.cov.55.cl.prop.men[,1]
vector.MAR.b.cov.55.cl.prop.women15.25.M.15.25 <- d.MAR.cov.55.cl.prop.women[,1]

vector.MAR.b.cov.55.cl.prop.men25.40.F.15.25 <- d.MAR.cov.55.cl.prop.men[,2]
vector.MAR.b.cov.55.cl.prop.women15.25.M.25.40 <- d.MAR.cov.55.cl.prop.women[,4]

vector.MAR.b.cov.55.cl.prop.men25.40.F.25.40 <- d.MAR.cov.55.cl.prop.men[,5]
vector.MAR.b.cov.55.cl.prop.women25.40.M.25.40 <- d.MAR.cov.55.cl.prop.women[,5]

vector.MAR.b.cov.55.cl.prop.men40.50.F.15.25 <- d.MAR.cov.55.cl.prop.men[,3]
vector.MAR.b.cov.55.cl.prop.women15.25.M.40.50 <- d.MAR.cov.55.cl.prop.women[,7]

vector.MAR.b.cov.55.cl.prop.men40.50.F.25.40 <- d.MAR.cov.55.cl.prop.men[,6]
vector.MAR.b.cov.55.cl.prop.women25.40.M.40.50 <- d.MAR.cov.55.cl.prop.women[,8]


# 60

d.MAR.cov.60 <- d.MAR %>%
  select(contains("cov.MAR.b.60.")) 
d.MAR.cov.60.cl.prop.men <-  d.MAR.cov.60 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.60.cl.prop.women <-  d.MAR.cov.60 %>%
  select(contains(".cl.prop.women"))

vector.MAR.b.cov.60.cl.prop.men15.25.F.15.25 <- d.MAR.cov.60.cl.prop.men[,1]
vector.MAR.b.cov.60.cl.prop.women15.25.M.15.25 <- d.MAR.cov.60.cl.prop.women[,1]

vector.MAR.b.cov.60.cl.prop.men25.40.F.15.25 <- d.MAR.cov.60.cl.prop.men[,2]
vector.MAR.b.cov.60.cl.prop.women15.25.M.25.40 <- d.MAR.cov.60.cl.prop.women[,4]

vector.MAR.b.cov.60.cl.prop.men25.40.F.25.40 <- d.MAR.cov.60.cl.prop.men[,5]
vector.MAR.b.cov.60.cl.prop.women25.40.M.25.40 <- d.MAR.cov.60.cl.prop.women[,5]

vector.MAR.b.cov.60.cl.prop.men40.50.F.15.25 <- d.MAR.cov.60.cl.prop.men[,3]
vector.MAR.b.cov.60.cl.prop.women15.25.M.40.50 <- d.MAR.cov.60.cl.prop.women[,7]

vector.MAR.b.cov.60.cl.prop.men40.50.F.25.40 <- d.MAR.cov.60.cl.prop.men[,6]
vector.MAR.b.cov.60.cl.prop.women25.40.M.40.50 <- d.MAR.cov.60.cl.prop.women[,8]


# 65

d.MAR.cov.65 <- d.MAR %>%
  select(contains("cov.MAR.b.65.")) 
d.MAR.cov.65.cl.prop.men <-  d.MAR.cov.65 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.65.cl.prop.women <-  d.MAR.cov.65 %>%
  select(contains(".cl.prop.women"))

vector.MAR.b.cov.65.cl.prop.men15.25.F.15.25 <- d.MAR.cov.65.cl.prop.men[,1]
vector.MAR.b.cov.65.cl.prop.women15.25.M.15.25 <- d.MAR.cov.65.cl.prop.women[,1]

vector.MAR.b.cov.65.cl.prop.men25.40.F.15.25 <- d.MAR.cov.65.cl.prop.men[,2]
vector.MAR.b.cov.65.cl.prop.women15.25.M.25.40 <- d.MAR.cov.65.cl.prop.women[,4]

vector.MAR.b.cov.65.cl.prop.men25.40.F.25.40 <- d.MAR.cov.65.cl.prop.men[,5]
vector.MAR.b.cov.65.cl.prop.women25.40.M.25.40 <- d.MAR.cov.65.cl.prop.women[,5]

vector.MAR.b.cov.65.cl.prop.men40.50.F.15.25 <- d.MAR.cov.65.cl.prop.men[,3]
vector.MAR.b.cov.65.cl.prop.women15.25.M.40.50 <- d.MAR.cov.65.cl.prop.women[,7]

vector.MAR.b.cov.65.cl.prop.men40.50.F.25.40 <- d.MAR.cov.65.cl.prop.men[,6]
vector.MAR.b.cov.65.cl.prop.women25.40.M.40.50 <- d.MAR.cov.65.cl.prop.women[,8]


# 70

d.MAR.cov.70 <- d.MAR %>%
  select(contains("cov.MAR.b.70.")) 
d.MAR.cov.70.cl.prop.men <-  d.MAR.cov.70 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.70.cl.prop.women <-  d.MAR.cov.70 %>%
  select(contains(".cl.prop.women"))


vector.MAR.b.cov.70.cl.prop.men15.25.F.15.25 <- d.MAR.cov.70.cl.prop.men[,1]
vector.MAR.b.cov.70.cl.prop.women15.25.M.15.25 <- d.MAR.cov.70.cl.prop.women[,1]

vector.MAR.b.cov.70.cl.prop.men25.40.F.15.25 <- d.MAR.cov.70.cl.prop.men[,2]
vector.MAR.b.cov.70.cl.prop.women15.25.M.25.40 <- d.MAR.cov.70.cl.prop.women[,4]

vector.MAR.b.cov.70.cl.prop.men25.40.F.25.40 <- d.MAR.cov.70.cl.prop.men[,5]
vector.MAR.b.cov.70.cl.prop.women25.40.M.25.40 <- d.MAR.cov.70.cl.prop.women[,5]

vector.MAR.b.cov.70.cl.prop.men40.50.F.15.25 <- d.MAR.cov.70.cl.prop.men[,3]
vector.MAR.b.cov.70.cl.prop.women15.25.M.40.50 <- d.MAR.cov.70.cl.prop.women[,7]

vector.MAR.b.cov.70.cl.prop.men40.50.F.25.40 <- d.MAR.cov.70.cl.prop.men[,6]
vector.MAR.b.cov.70.cl.prop.women25.40.M.40.50 <- d.MAR.cov.70.cl.prop.women[,8]


# 75

d.MAR.cov.75 <- d.MAR %>%
  select(contains("cov.MAR.b.75.")) 
d.MAR.cov.75.cl.prop.men <-  d.MAR.cov.75 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.75.cl.prop.women <-  d.MAR.cov.75 %>%
  select(contains(".cl.prop.women"))



vector.MAR.b.cov.75.cl.prop.men15.25.F.15.25 <- d.MAR.cov.75.cl.prop.men[,1]
vector.MAR.b.cov.75.cl.prop.women15.25.M.15.25 <- d.MAR.cov.75.cl.prop.women[,1]

vector.MAR.b.cov.75.cl.prop.men25.40.F.15.25 <- d.MAR.cov.75.cl.prop.men[,2]
vector.MAR.b.cov.75.cl.prop.women15.25.M.25.40 <- d.MAR.cov.75.cl.prop.women[,4]

vector.MAR.b.cov.75.cl.prop.men25.40.F.25.40 <- d.MAR.cov.75.cl.prop.men[,5]
vector.MAR.b.cov.75.cl.prop.women25.40.M.25.40 <- d.MAR.cov.75.cl.prop.women[,5]

vector.MAR.b.cov.75.cl.prop.men40.50.F.15.25 <- d.MAR.cov.75.cl.prop.men[,3]
vector.MAR.b.cov.75.cl.prop.women15.25.M.40.50 <- d.MAR.cov.75.cl.prop.women[,7]

vector.MAR.b.cov.75.cl.prop.men40.50.F.25.40 <- d.MAR.cov.75.cl.prop.men[,6]
vector.MAR.b.cov.75.cl.prop.women25.40.M.40.50 <- d.MAR.cov.75.cl.prop.women[,8]


# 80

d.MAR.cov.80 <- d.MAR %>%
  select(contains("cov.MAR.b.80.")) 
d.MAR.cov.80.cl.prop.men <-  d.MAR.cov.80 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.80.cl.prop.women <-  d.MAR.cov.80 %>%
  select(contains(".cl.prop.women"))


vector.MAR.b.cov.80.cl.prop.men15.25.F.15.25 <- d.MAR.cov.80.cl.prop.men[,1]
vector.MAR.b.cov.80.cl.prop.women15.25.M.15.25 <- d.MAR.cov.80.cl.prop.women[,1]

vector.MAR.b.cov.80.cl.prop.men25.40.F.15.25 <- d.MAR.cov.80.cl.prop.men[,2]
vector.MAR.b.cov.80.cl.prop.women15.25.M.25.40 <- d.MAR.cov.80.cl.prop.women[,4]

vector.MAR.b.cov.80.cl.prop.men25.40.F.25.40 <- d.MAR.cov.80.cl.prop.men[,5]
vector.MAR.b.cov.80.cl.prop.women25.40.M.25.40 <- d.MAR.cov.80.cl.prop.women[,5]

vector.MAR.b.cov.80.cl.prop.men40.50.F.15.25 <- d.MAR.cov.80.cl.prop.men[,3]
vector.MAR.b.cov.80.cl.prop.women15.25.M.40.50 <- d.MAR.cov.80.cl.prop.women[,7]

vector.MAR.b.cov.80.cl.prop.men40.50.F.25.40 <- d.MAR.cov.80.cl.prop.men[,6]
vector.MAR.b.cov.80.cl.prop.women25.40.M.40.50 <- d.MAR.cov.80.cl.prop.women[,8]


# 85

d.MAR.cov.85 <- d.MAR %>%
  select(contains("cov.MAR.b.85.")) 
d.MAR.cov.85.cl.prop.men <-  d.MAR.cov.85 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.85.cl.prop.women <-  d.MAR.cov.85 %>%
  select(contains(".cl.prop.women"))


vector.MAR.b.cov.85.cl.prop.men15.25.F.15.25 <- d.MAR.cov.85.cl.prop.men[,1]
vector.MAR.b.cov.85.cl.prop.women15.25.M.15.25 <- d.MAR.cov.85.cl.prop.women[,1]

vector.MAR.b.cov.85.cl.prop.men25.40.F.15.25 <- d.MAR.cov.85.cl.prop.men[,2]
vector.MAR.b.cov.85.cl.prop.women15.25.M.25.40 <- d.MAR.cov.85.cl.prop.women[,4]

vector.MAR.b.cov.85.cl.prop.men25.40.F.25.40 <- d.MAR.cov.85.cl.prop.men[,5]
vector.MAR.b.cov.85.cl.prop.women25.40.M.25.40 <- d.MAR.cov.85.cl.prop.women[,5]

vector.MAR.b.cov.85.cl.prop.men40.50.F.15.25 <- d.MAR.cov.85.cl.prop.men[,3]
vector.MAR.b.cov.85.cl.prop.women15.25.M.40.50 <- d.MAR.cov.85.cl.prop.women[,7]

vector.MAR.b.cov.85.cl.prop.men40.50.F.25.40 <- d.MAR.cov.85.cl.prop.men[,6]
vector.MAR.b.cov.85.cl.prop.women25.40.M.40.50 <- d.MAR.cov.85.cl.prop.women[,8]


# 90

d.MAR.cov.90 <- d.MAR %>%
  select(contains("cov.MAR.b.90.")) 
d.MAR.cov.90.cl.prop.men <-  d.MAR.cov.90 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.90.cl.prop.women <-  d.MAR.cov.90 %>%
  select(contains(".cl.prop.women"))


vector.MAR.b.cov.90.cl.prop.men15.25.F.15.25 <- d.MAR.cov.90.cl.prop.men[,1]
vector.MAR.b.cov.90.cl.prop.women15.25.M.15.25 <- d.MAR.cov.90.cl.prop.women[,1]

vector.MAR.b.cov.90.cl.prop.men25.40.F.15.25 <- d.MAR.cov.90.cl.prop.men[,2]
vector.MAR.b.cov.90.cl.prop.women15.25.M.25.40 <- d.MAR.cov.90.cl.prop.women[,4]

vector.MAR.b.cov.90.cl.prop.men25.40.F.25.40 <- d.MAR.cov.90.cl.prop.men[,5]
vector.MAR.b.cov.90.cl.prop.women25.40.M.25.40 <- d.MAR.cov.90.cl.prop.women[,5]

vector.MAR.b.cov.90.cl.prop.men40.50.F.15.25 <- d.MAR.cov.90.cl.prop.men[,3]
vector.MAR.b.cov.90.cl.prop.women15.25.M.40.50 <- d.MAR.cov.90.cl.prop.women[,7]

vector.MAR.b.cov.90.cl.prop.men40.50.F.25.40 <- d.MAR.cov.90.cl.prop.men[,6]
vector.MAR.b.cov.90.cl.prop.women25.40.M.40.50 <- d.MAR.cov.90.cl.prop.women[,8]


# 95

d.MAR.cov.95 <- d.MAR %>%
  select(contains("cov.MAR.b.95.")) 
d.MAR.cov.95.cl.prop.men <-  d.MAR.cov.95 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.95.cl.prop.women <-  d.MAR.cov.95 %>%
  select(contains(".cl.prop.women"))


vector.MAR.b.cov.95.cl.prop.men15.25.F.15.25 <- d.MAR.cov.95.cl.prop.men[,1]
vector.MAR.b.cov.95.cl.prop.women15.25.M.15.25 <- d.MAR.cov.95.cl.prop.women[,1]

vector.MAR.b.cov.95.cl.prop.men25.40.F.15.25 <- d.MAR.cov.95.cl.prop.men[,2]
vector.MAR.b.cov.95.cl.prop.women15.25.M.25.40 <- d.MAR.cov.95.cl.prop.women[,4]

vector.MAR.b.cov.95.cl.prop.men25.40.F.25.40 <- d.MAR.cov.95.cl.prop.men[,5]
vector.MAR.b.cov.95.cl.prop.women25.40.M.25.40 <- d.MAR.cov.95.cl.prop.women[,5]

vector.MAR.b.cov.95.cl.prop.men40.50.F.15.25 <- d.MAR.cov.95.cl.prop.men[,3]
vector.MAR.b.cov.95.cl.prop.women15.25.M.40.50 <- d.MAR.cov.95.cl.prop.women[,7]

vector.MAR.b.cov.95.cl.prop.men40.50.F.25.40 <- d.MAR.cov.95.cl.prop.men[,6]
vector.MAR.b.cov.95.cl.prop.women25.40.M.40.50 <- d.MAR.cov.95.cl.prop.women[,8]



# MAR - c


# 35

d.MAR.cov.35 <- d.MAR %>%
  select(contains("cov.MAR.c.35.")) 
d.MAR.cov.35.cl.prop.men <-  d.MAR.cov.35 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.35.cl.prop.women <-  d.MAR.cov.35 %>%
  select(contains(".cl.prop.women"))

vector.MAR.c.cov.35.cl.prop.men15.25.F.15.25 <- d.MAR.cov.35.cl.prop.men[,1]
vector.MAR.c.cov.35.cl.prop.women15.25.M.15.25 <- d.MAR.cov.35.cl.prop.women[,1]

vector.MAR.c.cov.35.cl.prop.men25.40.F.15.25 <- d.MAR.cov.35.cl.prop.men[,2]
vector.MAR.c.cov.35.cl.prop.women15.25.M.25.40 <- d.MAR.cov.35.cl.prop.women[,4]

vector.MAR.c.cov.35.cl.prop.men25.40.F.25.40 <- d.MAR.cov.35.cl.prop.men[,5]
vector.MAR.c.cov.35.cl.prop.women25.40.M.25.40 <- d.MAR.cov.35.cl.prop.women[,5]

vector.MAR.c.cov.35.cl.prop.men40.50.F.15.25 <- d.MAR.cov.35.cl.prop.men[,3]
vector.MAR.c.cov.35.cl.prop.women15.25.M.40.50 <- d.MAR.cov.35.cl.prop.women[,7]

vector.MAR.c.cov.35.cl.prop.men40.50.F.25.40 <- d.MAR.cov.35.cl.prop.men[,6]
vector.MAR.c.cov.35.cl.prop.women25.40.M.40.50 <- d.MAR.cov.35.cl.prop.women[,8]




# 40

d.MAR.cov.40 <- d.MAR %>%
  select(contains("cov.MAR.c.40.")) 
d.MAR.cov.40.cl.prop.men <-  d.MAR.cov.40 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.40.cl.prop.women <-  d.MAR.cov.40 %>%
  select(contains(".cl.prop.women"))

vector.MAR.c.cov.40.cl.prop.men15.25.F.15.25 <- d.MAR.cov.40.cl.prop.men[,1]
vector.MAR.c.cov.40.cl.prop.women15.25.M.15.25 <- d.MAR.cov.40.cl.prop.women[,1]

vector.MAR.c.cov.40.cl.prop.men25.40.F.15.25 <- d.MAR.cov.40.cl.prop.men[,2]
vector.MAR.c.cov.40.cl.prop.women15.25.M.25.40 <- d.MAR.cov.40.cl.prop.women[,4]

vector.MAR.c.cov.40.cl.prop.men25.40.F.25.40 <- d.MAR.cov.40.cl.prop.men[,5]
vector.MAR.c.cov.40.cl.prop.women25.40.M.25.40 <- d.MAR.cov.40.cl.prop.women[,5]

vector.MAR.c.cov.40.cl.prop.men40.50.F.15.25 <- d.MAR.cov.40.cl.prop.men[,3]
vector.MAR.c.cov.40.cl.prop.women15.25.M.40.50 <- d.MAR.cov.40.cl.prop.women[,7]

vector.MAR.c.cov.40.cl.prop.men40.50.F.25.40 <- d.MAR.cov.40.cl.prop.men[,6]
vector.MAR.c.cov.40.cl.prop.women25.40.M.40.50 <- d.MAR.cov.40.cl.prop.women[,8]


# 45

d.MAR.cov.45 <- d.MAR %>%
  select(contains("cov.MAR.c.45.")) 
d.MAR.cov.45.cl.prop.men <-  d.MAR.cov.45 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.45.cl.prop.women <-  d.MAR.cov.45 %>%
  select(contains(".cl.prop.women"))

vector.MAR.c.cov.45.cl.prop.men15.25.F.15.25 <- d.MAR.cov.45.cl.prop.men[,1]
vector.MAR.c.cov.45.cl.prop.women15.25.M.15.25 <- d.MAR.cov.45.cl.prop.women[,1]

vector.MAR.c.cov.45.cl.prop.men25.40.F.15.25 <- d.MAR.cov.45.cl.prop.men[,2]
vector.MAR.c.cov.45.cl.prop.women15.25.M.25.40 <- d.MAR.cov.45.cl.prop.women[,4]

vector.MAR.c.cov.45.cl.prop.men25.40.F.25.40 <- d.MAR.cov.45.cl.prop.men[,5]
vector.MAR.c.cov.45.cl.prop.women25.40.M.25.40 <- d.MAR.cov.45.cl.prop.women[,5]

vector.MAR.c.cov.45.cl.prop.men40.50.F.15.25 <- d.MAR.cov.45.cl.prop.men[,3]
vector.MAR.c.cov.45.cl.prop.women15.25.M.40.50 <- d.MAR.cov.45.cl.prop.women[,7]

vector.MAR.c.cov.45.cl.prop.men40.50.F.25.40 <- d.MAR.cov.45.cl.prop.men[,6]
vector.MAR.c.cov.45.cl.prop.women25.40.M.40.50 <- d.MAR.cov.45.cl.prop.women[,8]


# 50

d.MAR.cov.50 <- d.MAR %>%
  select(contains("cov.MAR.c.50.")) 
d.MAR.cov.50.cl.prop.men <-  d.MAR.cov.50 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.50.cl.prop.women <-  d.MAR.cov.50 %>%
  select(contains(".cl.prop.women"))


vector.MAR.c.cov.50.cl.prop.men15.25.F.15.25 <- d.MAR.cov.50.cl.prop.men[,1]
vector.MAR.c.cov.50.cl.prop.women15.25.M.15.25 <- d.MAR.cov.50.cl.prop.women[,1]

vector.MAR.c.cov.50.cl.prop.men25.40.F.15.25 <- d.MAR.cov.50.cl.prop.men[,2]
vector.MAR.c.cov.50.cl.prop.women15.25.M.25.40 <- d.MAR.cov.50.cl.prop.women[,4]

vector.MAR.c.cov.50.cl.prop.men25.40.F.25.40 <- d.MAR.cov.50.cl.prop.men[,5]
vector.MAR.c.cov.50.cl.prop.women25.40.M.25.40 <- d.MAR.cov.50.cl.prop.women[,5]

vector.MAR.c.cov.50.cl.prop.men40.50.F.15.25 <- d.MAR.cov.50.cl.prop.men[,3]
vector.MAR.c.cov.50.cl.prop.women15.25.M.40.50 <- d.MAR.cov.50.cl.prop.women[,7]

vector.MAR.c.cov.50.cl.prop.men40.50.F.25.40 <- d.MAR.cov.50.cl.prop.men[,6]
vector.MAR.c.cov.50.cl.prop.women25.40.M.40.50 <- d.MAR.cov.50.cl.prop.women[,8]


# 55

d.MAR.cov.55 <- d.MAR %>%
  select(contains("cov.MAR.c.55.")) 
d.MAR.cov.55.cl.prop.men <-  d.MAR.cov.55 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.55.cl.prop.women <-  d.MAR.cov.55 %>%
  select(contains(".cl.prop.women"))

vector.MAR.c.cov.55.cl.prop.men15.25.F.15.25 <- d.MAR.cov.55.cl.prop.men[,1]
vector.MAR.c.cov.55.cl.prop.women15.25.M.15.25 <- d.MAR.cov.55.cl.prop.women[,1]

vector.MAR.c.cov.55.cl.prop.men25.40.F.15.25 <- d.MAR.cov.55.cl.prop.men[,2]
vector.MAR.c.cov.55.cl.prop.women15.25.M.25.40 <- d.MAR.cov.55.cl.prop.women[,4]

vector.MAR.c.cov.55.cl.prop.men25.40.F.25.40 <- d.MAR.cov.55.cl.prop.men[,5]
vector.MAR.c.cov.55.cl.prop.women25.40.M.25.40 <- d.MAR.cov.55.cl.prop.women[,5]

vector.MAR.c.cov.55.cl.prop.men40.50.F.15.25 <- d.MAR.cov.55.cl.prop.men[,3]
vector.MAR.c.cov.55.cl.prop.women15.25.M.40.50 <- d.MAR.cov.55.cl.prop.women[,7]

vector.MAR.c.cov.55.cl.prop.men40.50.F.25.40 <- d.MAR.cov.55.cl.prop.men[,6]
vector.MAR.c.cov.55.cl.prop.women25.40.M.40.50 <- d.MAR.cov.55.cl.prop.women[,8]


# 60

d.MAR.cov.60 <- d.MAR %>%
  select(contains("cov.MAR.c.60.")) 
d.MAR.cov.60.cl.prop.men <-  d.MAR.cov.60 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.60.cl.prop.women <-  d.MAR.cov.60 %>%
  select(contains(".cl.prop.women"))

vector.MAR.c.cov.60.cl.prop.men15.25.F.15.25 <- d.MAR.cov.60.cl.prop.men[,1]
vector.MAR.c.cov.60.cl.prop.women15.25.M.15.25 <- d.MAR.cov.60.cl.prop.women[,1]

vector.MAR.c.cov.60.cl.prop.men25.40.F.15.25 <- d.MAR.cov.60.cl.prop.men[,2]
vector.MAR.c.cov.60.cl.prop.women15.25.M.25.40 <- d.MAR.cov.60.cl.prop.women[,4]

vector.MAR.c.cov.60.cl.prop.men25.40.F.25.40 <- d.MAR.cov.60.cl.prop.men[,5]
vector.MAR.c.cov.60.cl.prop.women25.40.M.25.40 <- d.MAR.cov.60.cl.prop.women[,5]

vector.MAR.c.cov.60.cl.prop.men40.50.F.15.25 <- d.MAR.cov.60.cl.prop.men[,3]
vector.MAR.c.cov.60.cl.prop.women15.25.M.40.50 <- d.MAR.cov.60.cl.prop.women[,7]

vector.MAR.c.cov.60.cl.prop.men40.50.F.25.40 <- d.MAR.cov.60.cl.prop.men[,6]
vector.MAR.c.cov.60.cl.prop.women25.40.M.40.50 <- d.MAR.cov.60.cl.prop.women[,8]


# 65

d.MAR.cov.65 <- d.MAR %>%
  select(contains("cov.MAR.c.65.")) 
d.MAR.cov.65.cl.prop.men <-  d.MAR.cov.65 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.65.cl.prop.women <-  d.MAR.cov.65 %>%
  select(contains(".cl.prop.women"))

vector.MAR.c.cov.65.cl.prop.men15.25.F.15.25 <- d.MAR.cov.65.cl.prop.men[,1]
vector.MAR.c.cov.65.cl.prop.women15.25.M.15.25 <- d.MAR.cov.65.cl.prop.women[,1]

vector.MAR.c.cov.65.cl.prop.men25.40.F.15.25 <- d.MAR.cov.65.cl.prop.men[,2]
vector.MAR.c.cov.65.cl.prop.women15.25.M.25.40 <- d.MAR.cov.65.cl.prop.women[,4]

vector.MAR.c.cov.65.cl.prop.men25.40.F.25.40 <- d.MAR.cov.65.cl.prop.men[,5]
vector.MAR.c.cov.65.cl.prop.women25.40.M.25.40 <- d.MAR.cov.65.cl.prop.women[,5]

vector.MAR.c.cov.65.cl.prop.men40.50.F.15.25 <- d.MAR.cov.65.cl.prop.men[,3]
vector.MAR.c.cov.65.cl.prop.women15.25.M.40.50 <- d.MAR.cov.65.cl.prop.women[,7]

vector.MAR.c.cov.65.cl.prop.men40.50.F.25.40 <- d.MAR.cov.65.cl.prop.men[,6]
vector.MAR.c.cov.65.cl.prop.women25.40.M.40.50 <- d.MAR.cov.65.cl.prop.women[,8]


# 70

d.MAR.cov.70 <- d.MAR %>%
  select(contains("cov.MAR.c.70.")) 
d.MAR.cov.70.cl.prop.men <-  d.MAR.cov.70 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.70.cl.prop.women <-  d.MAR.cov.70 %>%
  select(contains(".cl.prop.women"))


vector.MAR.c.cov.70.cl.prop.men15.25.F.15.25 <- d.MAR.cov.70.cl.prop.men[,1]
vector.MAR.c.cov.70.cl.prop.women15.25.M.15.25 <- d.MAR.cov.70.cl.prop.women[,1]

vector.MAR.c.cov.70.cl.prop.men25.40.F.15.25 <- d.MAR.cov.70.cl.prop.men[,2]
vector.MAR.c.cov.70.cl.prop.women15.25.M.25.40 <- d.MAR.cov.70.cl.prop.women[,4]

vector.MAR.c.cov.70.cl.prop.men25.40.F.25.40 <- d.MAR.cov.70.cl.prop.men[,5]
vector.MAR.c.cov.70.cl.prop.women25.40.M.25.40 <- d.MAR.cov.70.cl.prop.women[,5]

vector.MAR.c.cov.70.cl.prop.men40.50.F.15.25 <- d.MAR.cov.70.cl.prop.men[,3]
vector.MAR.c.cov.70.cl.prop.women15.25.M.40.50 <- d.MAR.cov.70.cl.prop.women[,7]

vector.MAR.c.cov.70.cl.prop.men40.50.F.25.40 <- d.MAR.cov.70.cl.prop.men[,6]
vector.MAR.c.cov.70.cl.prop.women25.40.M.40.50 <- d.MAR.cov.70.cl.prop.women[,8]


# 75

d.MAR.cov.75 <- d.MAR %>%
  select(contains("cov.MAR.c.75.")) 
d.MAR.cov.75.cl.prop.men <-  d.MAR.cov.75 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.75.cl.prop.women <-  d.MAR.cov.75 %>%
  select(contains(".cl.prop.women"))



vector.MAR.c.cov.75.cl.prop.men15.25.F.15.25 <- d.MAR.cov.75.cl.prop.men[,1]
vector.MAR.c.cov.75.cl.prop.women15.25.M.15.25 <- d.MAR.cov.75.cl.prop.women[,1]

vector.MAR.c.cov.75.cl.prop.men25.40.F.15.25 <- d.MAR.cov.75.cl.prop.men[,2]
vector.MAR.c.cov.75.cl.prop.women15.25.M.25.40 <- d.MAR.cov.75.cl.prop.women[,4]

vector.MAR.c.cov.75.cl.prop.men25.40.F.25.40 <- d.MAR.cov.75.cl.prop.men[,5]
vector.MAR.c.cov.75.cl.prop.women25.40.M.25.40 <- d.MAR.cov.75.cl.prop.women[,5]

vector.MAR.c.cov.75.cl.prop.men40.50.F.15.25 <- d.MAR.cov.75.cl.prop.men[,3]
vector.MAR.c.cov.75.cl.prop.women15.25.M.40.50 <- d.MAR.cov.75.cl.prop.women[,7]

vector.MAR.c.cov.75.cl.prop.men40.50.F.25.40 <- d.MAR.cov.75.cl.prop.men[,6]
vector.MAR.c.cov.75.cl.prop.women25.40.M.40.50 <- d.MAR.cov.75.cl.prop.women[,8]


# 80

d.MAR.cov.80 <- d.MAR %>%
  select(contains("cov.MAR.c.80.")) 
d.MAR.cov.80.cl.prop.men <-  d.MAR.cov.80 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.80.cl.prop.women <-  d.MAR.cov.80 %>%
  select(contains(".cl.prop.women"))


vector.MAR.c.cov.80.cl.prop.men15.25.F.15.25 <- d.MAR.cov.80.cl.prop.men[,1]
vector.MAR.c.cov.80.cl.prop.women15.25.M.15.25 <- d.MAR.cov.80.cl.prop.women[,1]

vector.MAR.c.cov.80.cl.prop.men25.40.F.15.25 <- d.MAR.cov.80.cl.prop.men[,2]
vector.MAR.c.cov.80.cl.prop.women15.25.M.25.40 <- d.MAR.cov.80.cl.prop.women[,4]

vector.MAR.c.cov.80.cl.prop.men25.40.F.25.40 <- d.MAR.cov.80.cl.prop.men[,5]
vector.MAR.c.cov.80.cl.prop.women25.40.M.25.40 <- d.MAR.cov.80.cl.prop.women[,5]

vector.MAR.c.cov.80.cl.prop.men40.50.F.15.25 <- d.MAR.cov.80.cl.prop.men[,3]
vector.MAR.c.cov.80.cl.prop.women15.25.M.40.50 <- d.MAR.cov.80.cl.prop.women[,7]

vector.MAR.c.cov.80.cl.prop.men40.50.F.25.40 <- d.MAR.cov.80.cl.prop.men[,6]
vector.MAR.c.cov.80.cl.prop.women25.40.M.40.50 <- d.MAR.cov.80.cl.prop.women[,8]


# 85

d.MAR.cov.85 <- d.MAR %>%
  select(contains("cov.MAR.c.85.")) 
d.MAR.cov.85.cl.prop.men <-  d.MAR.cov.85 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.85.cl.prop.women <-  d.MAR.cov.85 %>%
  select(contains(".cl.prop.women"))


vector.MAR.c.cov.85.cl.prop.men15.25.F.15.25 <- d.MAR.cov.85.cl.prop.men[,1]
vector.MAR.c.cov.85.cl.prop.women15.25.M.15.25 <- d.MAR.cov.85.cl.prop.women[,1]

vector.MAR.c.cov.85.cl.prop.men25.40.F.15.25 <- d.MAR.cov.85.cl.prop.men[,2]
vector.MAR.c.cov.85.cl.prop.women15.25.M.25.40 <- d.MAR.cov.85.cl.prop.women[,4]

vector.MAR.c.cov.85.cl.prop.men25.40.F.25.40 <- d.MAR.cov.85.cl.prop.men[,5]
vector.MAR.c.cov.85.cl.prop.women25.40.M.25.40 <- d.MAR.cov.85.cl.prop.women[,5]

vector.MAR.c.cov.85.cl.prop.men40.50.F.15.25 <- d.MAR.cov.85.cl.prop.men[,3]
vector.MAR.c.cov.85.cl.prop.women15.25.M.40.50 <- d.MAR.cov.85.cl.prop.women[,7]

vector.MAR.c.cov.85.cl.prop.men40.50.F.25.40 <- d.MAR.cov.85.cl.prop.men[,6]
vector.MAR.c.cov.85.cl.prop.women25.40.M.40.50 <- d.MAR.cov.85.cl.prop.women[,8]


# 90

d.MAR.cov.90 <- d.MAR %>%
  select(contains("cov.MAR.c.90.")) 
d.MAR.cov.90.cl.prop.men <-  d.MAR.cov.90 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.90.cl.prop.women <-  d.MAR.cov.90 %>%
  select(contains(".cl.prop.women"))


vector.MAR.c.cov.90.cl.prop.men15.25.F.15.25 <- d.MAR.cov.90.cl.prop.men[,1]
vector.MAR.c.cov.90.cl.prop.women15.25.M.15.25 <- d.MAR.cov.90.cl.prop.women[,1]

vector.MAR.c.cov.90.cl.prop.men25.40.F.15.25 <- d.MAR.cov.90.cl.prop.men[,2]
vector.MAR.c.cov.90.cl.prop.women15.25.M.25.40 <- d.MAR.cov.90.cl.prop.women[,4]

vector.MAR.c.cov.90.cl.prop.men25.40.F.25.40 <- d.MAR.cov.90.cl.prop.men[,5]
vector.MAR.c.cov.90.cl.prop.women25.40.M.25.40 <- d.MAR.cov.90.cl.prop.women[,5]

vector.MAR.c.cov.90.cl.prop.men40.50.F.15.25 <- d.MAR.cov.90.cl.prop.men[,3]
vector.MAR.c.cov.90.cl.prop.women15.25.M.40.50 <- d.MAR.cov.90.cl.prop.women[,7]

vector.MAR.c.cov.90.cl.prop.men40.50.F.25.40 <- d.MAR.cov.90.cl.prop.men[,6]
vector.MAR.c.cov.90.cl.prop.women25.40.M.40.50 <- d.MAR.cov.90.cl.prop.women[,8]


# 95

d.MAR.cov.95 <- d.MAR %>%
  select(contains("cov.MAR.c.95.")) 
d.MAR.cov.95.cl.prop.men <-  d.MAR.cov.95 %>%
  select(contains(".cl.prop.men")) #
d.MAR.cov.95.cl.prop.women <-  d.MAR.cov.95 %>%
  select(contains(".cl.prop.women"))


vector.MAR.c.cov.95.cl.prop.men15.25.F.15.25 <- d.MAR.cov.95.cl.prop.men[,1]
vector.MAR.c.cov.95.cl.prop.women15.25.M.15.25 <- d.MAR.cov.95.cl.prop.women[,1]

vector.MAR.c.cov.95.cl.prop.men25.40.F.15.25 <- d.MAR.cov.95.cl.prop.men[,2]
vector.MAR.c.cov.95.cl.prop.women15.25.M.25.40 <- d.MAR.cov.95.cl.prop.women[,4]

vector.MAR.c.cov.95.cl.prop.men25.40.F.25.40 <- d.MAR.cov.95.cl.prop.men[,5]
vector.MAR.c.cov.95.cl.prop.women25.40.M.25.40 <- d.MAR.cov.95.cl.prop.women[,5]

vector.MAR.c.cov.95.cl.prop.men40.50.F.15.25 <- d.MAR.cov.95.cl.prop.men[,3]
vector.MAR.c.cov.95.cl.prop.women15.25.M.40.50 <- d.MAR.cov.95.cl.prop.women[,7]

vector.MAR.c.cov.95.cl.prop.men40.50.F.25.40 <- d.MAR.cov.95.cl.prop.men[,6]
vector.MAR.c.cov.95.cl.prop.women25.40.M.40.50 <- d.MAR.cov.95.cl.prop.women[,8]





# II. Age Difference ----------------

# MCAR

d.MCAR <- dr %>%
  select(contains("MCAR."))


# 35


d.MCAR.cov.35 <- d.MCAR %>%
  select(contains("cov.MCAR.35.")) # MAR - a
AD.MCAR.cov.35 <- d.MCAR.cov.35 %>%
  select(contains(".AD.")) 

# Mean

vector.mean.MCAR.cov.35.AD.women.cl.15.25 <- AD.MCAR.cov.35[,1]
vector.mean.MCAR.cov.35.AD.men.cl.15.25 <- AD.MCAR.cov.35[,2]

vector.mean.MCAR.cov.35.AD.women.cl.25.40 <- AD.MCAR.cov.35[,3]
vector.mean.MCAR.cov.35.AD.men.cl.25.40 <- AD.MCAR.cov.35[,4]

vector.mean.MCAR.cov.35.AD.women.cl.40.50 <- AD.MCAR.cov.35[,5]
vector.mean.MCAR.cov.35.AD.men.cl.40.50 <- AD.MCAR.cov.35[,6]

# Median

vector.med.MCAR.cov.35.AD.women.cl.15.25 <- AD.MCAR.cov.35[,7]
vector.med.MCAR.cov.35.AD.men.cl.15.25 <- AD.MCAR.cov.35[,8]

vector.med.MCAR.cov.35.AD.women.cl.25.40 <- AD.MCAR.cov.35[,9]
vector.med.MCAR.cov.35.AD.men.cl.25.40 <- AD.MCAR.cov.35[,10]

vector.med.MCAR.cov.35.AD.women.cl.40.50 <- AD.MCAR.cov.35[,11]
vector.med.MCAR.cov.35.AD.men.cl.40.50 <- AD.MCAR.cov.35[,12]

# Standard deviation

vector.sd.MCAR.cov.35.AD.women.cl.15.25 <- AD.MCAR.cov.35[,13]
vector.sd.MCAR.cov.35.AD.men.cl.15.25 <- AD.MCAR.cov.35[,14]

vector.sd.MCAR.cov.35.AD.women.cl.25.40 <- AD.MCAR.cov.35[,15]
vector.sd.MCAR.cov.35.AD.men.cl.25.40 <- AD.MCAR.cov.35[,16]

vector.sd.MCAR.cov.35.AD.women.cl.40.50 <- AD.MCAR.cov.35[,17]
vector.sd.MCAR.cov.35.AD.men.cl.40.50 <- AD.MCAR.cov.35[,18]


# 40

d.MCAR.cov.40 <- d.MCAR %>%
  select(contains("cov.MCAR.40.")) # MAR - a
AD.MCAR.cov.40 <- d.MCAR.cov.40 %>%
  select(contains(".AD.")) 

# Mean

vector.mean.MCAR.cov.40.AD.women.cl.15.25 <- AD.MCAR.cov.40[,1]
vector.mean.MCAR.cov.40.AD.men.cl.15.25 <- AD.MCAR.cov.40[,2]

vector.mean.MCAR.cov.40.AD.women.cl.25.40 <- AD.MCAR.cov.40[,3]
vector.mean.MCAR.cov.40.AD.men.cl.25.40 <- AD.MCAR.cov.40[,4]

vector.mean.MCAR.cov.40.AD.women.cl.40.50 <- AD.MCAR.cov.40[,5]
vector.mean.MCAR.cov.40.AD.men.cl.40.50 <- AD.MCAR.cov.40[,6]

# Median

vector.med.MCAR.cov.40.AD.women.cl.15.25 <- AD.MCAR.cov.40[,7]
vector.med.MCAR.cov.40.AD.men.cl.15.25 <- AD.MCAR.cov.40[,8]

vector.med.MCAR.cov.40.AD.women.cl.25.40 <- AD.MCAR.cov.40[,9]
vector.med.MCAR.cov.40.AD.men.cl.25.40 <- AD.MCAR.cov.40[,10]

vector.med.MCAR.cov.40.AD.women.cl.40.50 <- AD.MCAR.cov.40[,11]
vector.med.MCAR.cov.40.AD.men.cl.40.50 <- AD.MCAR.cov.40[,12]

# Standard deviation

vector.sd.MCAR.cov.40.AD.women.cl.15.25 <- AD.MCAR.cov.40[,13]
vector.sd.MCAR.cov.40.AD.men.cl.15.25 <- AD.MCAR.cov.40[,14]

vector.sd.MCAR.cov.40.AD.women.cl.25.40 <- AD.MCAR.cov.40[,15]
vector.sd.MCAR.cov.40.AD.men.cl.25.40 <- AD.MCAR.cov.40[,16]

vector.sd.MCAR.cov.40.AD.women.cl.40.50 <- AD.MCAR.cov.40[,17]
vector.sd.MCAR.cov.40.AD.men.cl.40.50 <- AD.MCAR.cov.40[,18]


# 45

d.MCAR.cov.45 <- d.MCAR %>%
  select(contains("cov.MCAR.45.")) # MAR - a
AD.MCAR.cov.45 <- d.MCAR.cov.45 %>%
  select(contains(".AD.")) 

# Mean

vector.mean.MCAR.cov.45.AD.women.cl.15.25 <- AD.MCAR.cov.45[,1]
vector.mean.MCAR.cov.45.AD.men.cl.15.25 <- AD.MCAR.cov.45[,2]

vector.mean.MCAR.cov.45.AD.women.cl.25.40 <- AD.MCAR.cov.45[,3]
vector.mean.MCAR.cov.45.AD.men.cl.25.40 <- AD.MCAR.cov.45[,4]

vector.mean.MCAR.cov.45.AD.women.cl.40.50 <- AD.MCAR.cov.45[,5]
vector.mean.MCAR.cov.45.AD.men.cl.40.50 <- AD.MCAR.cov.45[,6]

# Median

vector.med.MCAR.cov.45.AD.women.cl.15.25 <- AD.MCAR.cov.45[,7]
vector.med.MCAR.cov.45.AD.men.cl.15.25 <- AD.MCAR.cov.45[,8]

vector.med.MCAR.cov.45.AD.women.cl.25.40 <- AD.MCAR.cov.45[,9]
vector.med.MCAR.cov.45.AD.men.cl.25.40 <- AD.MCAR.cov.45[,10]

vector.med.MCAR.cov.45.AD.women.cl.40.50 <- AD.MCAR.cov.45[,11]
vector.med.MCAR.cov.45.AD.men.cl.40.50 <- AD.MCAR.cov.45[,12]

# Standard deviation

vector.sd.MCAR.cov.45.AD.women.cl.15.25 <- AD.MCAR.cov.45[,13]
vector.sd.MCAR.cov.45.AD.men.cl.15.25 <- AD.MCAR.cov.45[,14]

vector.sd.MCAR.cov.45.AD.women.cl.25.40 <- AD.MCAR.cov.45[,15]
vector.sd.MCAR.cov.45.AD.men.cl.25.40 <- AD.MCAR.cov.45[,16]

vector.sd.MCAR.cov.45.AD.women.cl.40.50 <- AD.MCAR.cov.45[,17]
vector.sd.MCAR.cov.45.AD.men.cl.40.50 <- AD.MCAR.cov.45[,18]



# 50

d.MCAR.cov.50 <- d.MCAR %>%
  select(contains("cov.MCAR.50.")) # MAR - a
AD.MCAR.cov.50 <- d.MCAR.cov.50 %>%
  select(contains(".AD.")) 

# Mean

vector.mean.MCAR.cov.50.AD.women.cl.15.25 <- AD.MCAR.cov.50[,1]
vector.mean.MCAR.cov.50.AD.men.cl.15.25 <- AD.MCAR.cov.50[,2]

vector.mean.MCAR.cov.50.AD.women.cl.25.40 <- AD.MCAR.cov.50[,3]
vector.mean.MCAR.cov.50.AD.men.cl.25.40 <- AD.MCAR.cov.50[,4]

vector.mean.MCAR.cov.50.AD.women.cl.40.50 <- AD.MCAR.cov.50[,5]
vector.mean.MCAR.cov.50.AD.men.cl.40.50 <- AD.MCAR.cov.50[,6]

# Median

vector.med.MCAR.cov.50.AD.women.cl.15.25 <- AD.MCAR.cov.50[,7]
vector.med.MCAR.cov.50.AD.men.cl.15.25 <- AD.MCAR.cov.50[,8]

vector.med.MCAR.cov.50.AD.women.cl.25.40 <- AD.MCAR.cov.50[,9]
vector.med.MCAR.cov.50.AD.men.cl.25.40 <- AD.MCAR.cov.50[,10]

vector.med.MCAR.cov.50.AD.women.cl.40.50 <- AD.MCAR.cov.50[,11]
vector.med.MCAR.cov.50.AD.men.cl.40.50 <- AD.MCAR.cov.50[,12]

# Standard deviation

vector.sd.MCAR.cov.50.AD.women.cl.15.25 <- AD.MCAR.cov.50[,13]
vector.sd.MCAR.cov.50.AD.men.cl.15.25 <- AD.MCAR.cov.50[,14]

vector.sd.MCAR.cov.50.AD.women.cl.25.40 <- AD.MCAR.cov.50[,15]
vector.sd.MCAR.cov.50.AD.men.cl.25.40 <- AD.MCAR.cov.50[,16]

vector.sd.MCAR.cov.50.AD.women.cl.40.50 <- AD.MCAR.cov.50[,17]
vector.sd.MCAR.cov.50.AD.men.cl.40.50 <- AD.MCAR.cov.50[,18]



# 55

d.MCAR.cov.55 <- d.MCAR %>%
  select(contains("cov.MCAR.55.")) # MAR - a
AD.MCAR.cov.55 <- d.MCAR.cov.55 %>%
  select(contains(".AD.")) 

# Mean

vector.mean.MCAR.cov.55.AD.women.cl.15.25 <- AD.MCAR.cov.55[,1]
vector.mean.MCAR.cov.55.AD.men.cl.15.25 <- AD.MCAR.cov.55[,2]

vector.mean.MCAR.cov.55.AD.women.cl.25.40 <- AD.MCAR.cov.55[,3]
vector.mean.MCAR.cov.55.AD.men.cl.25.40 <- AD.MCAR.cov.55[,4]

vector.mean.MCAR.cov.55.AD.women.cl.40.50 <- AD.MCAR.cov.55[,5]
vector.mean.MCAR.cov.55.AD.men.cl.40.50 <- AD.MCAR.cov.55[,6]

# Median

vector.med.MCAR.cov.55.AD.women.cl.15.25 <- AD.MCAR.cov.55[,7]
vector.med.MCAR.cov.55.AD.men.cl.15.25 <- AD.MCAR.cov.55[,8]

vector.med.MCAR.cov.55.AD.women.cl.25.40 <- AD.MCAR.cov.55[,9]
vector.med.MCAR.cov.55.AD.men.cl.25.40 <- AD.MCAR.cov.55[,10]

vector.med.MCAR.cov.55.AD.women.cl.40.50 <- AD.MCAR.cov.55[,11]
vector.med.MCAR.cov.55.AD.men.cl.40.50 <- AD.MCAR.cov.55[,12]

# Standard deviation

vector.sd.MCAR.cov.55.AD.women.cl.15.25 <- AD.MCAR.cov.55[,13]
vector.sd.MCAR.cov.55.AD.men.cl.15.25 <- AD.MCAR.cov.55[,14]

vector.sd.MCAR.cov.55.AD.women.cl.25.40 <- AD.MCAR.cov.55[,15]
vector.sd.MCAR.cov.55.AD.men.cl.25.40 <- AD.MCAR.cov.55[,16]

vector.sd.MCAR.cov.55.AD.women.cl.40.50 <- AD.MCAR.cov.55[,17]
vector.sd.MCAR.cov.55.AD.men.cl.40.50 <- AD.MCAR.cov.55[,18]


# 60

d.MCAR.cov.60 <- d.MCAR %>%
  select(contains("cov.MCAR.60.")) # MAR - a
AD.MCAR.cov.60 <- d.MCAR.cov.60 %>%
  select(contains(".AD.")) 

# Mean

vector.mean.MCAR.cov.60.AD.women.cl.15.25 <- AD.MCAR.cov.60[,1]
vector.mean.MCAR.cov.60.AD.men.cl.15.25 <- AD.MCAR.cov.60[,2]

vector.mean.MCAR.cov.60.AD.women.cl.25.40 <- AD.MCAR.cov.60[,3]
vector.mean.MCAR.cov.60.AD.men.cl.25.40 <- AD.MCAR.cov.60[,4]

vector.mean.MCAR.cov.60.AD.women.cl.40.50 <- AD.MCAR.cov.60[,5]
vector.mean.MCAR.cov.60.AD.men.cl.40.50 <- AD.MCAR.cov.60[,6]

# Median

vector.med.MCAR.cov.60.AD.women.cl.15.25 <- AD.MCAR.cov.60[,7]
vector.med.MCAR.cov.60.AD.men.cl.15.25 <- AD.MCAR.cov.60[,8]

vector.med.MCAR.cov.60.AD.women.cl.25.40 <- AD.MCAR.cov.60[,9]
vector.med.MCAR.cov.60.AD.men.cl.25.40 <- AD.MCAR.cov.60[,10]

vector.med.MCAR.cov.60.AD.women.cl.40.50 <- AD.MCAR.cov.60[,11]
vector.med.MCAR.cov.60.AD.men.cl.40.50 <- AD.MCAR.cov.60[,12]

# Standard deviation

vector.sd.MCAR.cov.60.AD.women.cl.15.25 <- AD.MCAR.cov.60[,13]
vector.sd.MCAR.cov.60.AD.men.cl.15.25 <- AD.MCAR.cov.60[,14]

vector.sd.MCAR.cov.60.AD.women.cl.25.40 <- AD.MCAR.cov.60[,15]
vector.sd.MCAR.cov.60.AD.men.cl.25.40 <- AD.MCAR.cov.60[,16]

vector.sd.MCAR.cov.60.AD.women.cl.40.50 <- AD.MCAR.cov.60[,17]
vector.sd.MCAR.cov.60.AD.men.cl.40.50 <- AD.MCAR.cov.60[,18]


# 65

d.MCAR.cov.65 <- d.MCAR %>%
  select(contains("cov.MCAR.65.")) # MAR - a
AD.MCAR.cov.65 <- d.MCAR.cov.65 %>%
  select(contains(".AD.")) 

# Mean

vector.mean.MCAR.cov.65.AD.women.cl.15.25 <- AD.MCAR.cov.65[,1]
vector.mean.MCAR.cov.65.AD.men.cl.15.25 <- AD.MCAR.cov.65[,2]

vector.mean.MCAR.cov.65.AD.women.cl.25.40 <- AD.MCAR.cov.65[,3]
vector.mean.MCAR.cov.65.AD.men.cl.25.40 <- AD.MCAR.cov.65[,4]

vector.mean.MCAR.cov.65.AD.women.cl.40.50 <- AD.MCAR.cov.65[,5]
vector.mean.MCAR.cov.65.AD.men.cl.40.50 <- AD.MCAR.cov.65[,6]

# Median

vector.med.MCAR.cov.65.AD.women.cl.15.25 <- AD.MCAR.cov.65[,7]
vector.med.MCAR.cov.65.AD.men.cl.15.25 <- AD.MCAR.cov.65[,8]

vector.med.MCAR.cov.65.AD.women.cl.25.40 <- AD.MCAR.cov.65[,9]
vector.med.MCAR.cov.65.AD.men.cl.25.40 <- AD.MCAR.cov.65[,10]

vector.med.MCAR.cov.65.AD.women.cl.40.50 <- AD.MCAR.cov.65[,11]
vector.med.MCAR.cov.65.AD.men.cl.40.50 <- AD.MCAR.cov.65[,12]

# Standard deviation

vector.sd.MCAR.cov.65.AD.women.cl.15.25 <- AD.MCAR.cov.65[,13]
vector.sd.MCAR.cov.65.AD.men.cl.15.25 <- AD.MCAR.cov.65[,14]

vector.sd.MCAR.cov.65.AD.women.cl.25.40 <- AD.MCAR.cov.65[,15]
vector.sd.MCAR.cov.65.AD.men.cl.25.40 <- AD.MCAR.cov.65[,16]

vector.sd.MCAR.cov.65.AD.women.cl.40.50 <- AD.MCAR.cov.65[,17]
vector.sd.MCAR.cov.65.AD.men.cl.40.50 <- AD.MCAR.cov.65[,18]


# 70

d.MCAR.cov.70 <- d.MCAR %>%
  select(contains("cov.MCAR.70.")) # MAR - a
AD.MCAR.cov.70 <- d.MCAR.cov.70 %>%
  select(contains(".AD.")) 

# Mean

vector.mean.MCAR.cov.70.AD.women.cl.15.25 <- AD.MCAR.cov.70[,1]
vector.mean.MCAR.cov.70.AD.men.cl.15.25 <- AD.MCAR.cov.70[,2]

vector.mean.MCAR.cov.70.AD.women.cl.25.40 <- AD.MCAR.cov.70[,3]
vector.mean.MCAR.cov.70.AD.men.cl.25.40 <- AD.MCAR.cov.70[,4]

vector.mean.MCAR.cov.70.AD.women.cl.40.50 <- AD.MCAR.cov.70[,5]
vector.mean.MCAR.cov.70.AD.men.cl.40.50 <- AD.MCAR.cov.70[,6]

# Median

vector.med.MCAR.cov.70.AD.women.cl.15.25 <- AD.MCAR.cov.70[,7]
vector.med.MCAR.cov.70.AD.men.cl.15.25 <- AD.MCAR.cov.70[,8]

vector.med.MCAR.cov.70.AD.women.cl.25.40 <- AD.MCAR.cov.70[,9]
vector.med.MCAR.cov.70.AD.men.cl.25.40 <- AD.MCAR.cov.70[,10]

vector.med.MCAR.cov.70.AD.women.cl.40.50 <- AD.MCAR.cov.70[,11]
vector.med.MCAR.cov.70.AD.men.cl.40.50 <- AD.MCAR.cov.70[,12]

# Standard deviation

vector.sd.MCAR.cov.70.AD.women.cl.15.25 <- AD.MCAR.cov.70[,13]
vector.sd.MCAR.cov.70.AD.men.cl.15.25 <- AD.MCAR.cov.70[,14]

vector.sd.MCAR.cov.70.AD.women.cl.25.40 <- AD.MCAR.cov.70[,15]
vector.sd.MCAR.cov.70.AD.men.cl.25.40 <- AD.MCAR.cov.70[,16]

vector.sd.MCAR.cov.70.AD.women.cl.40.50 <- AD.MCAR.cov.70[,17]
vector.sd.MCAR.cov.70.AD.men.cl.40.50 <- AD.MCAR.cov.70[,18]


# 75

d.MCAR.cov.75 <- d.MCAR %>%
  select(contains("cov.MCAR.75.")) # MAR - a
AD.MCAR.cov.75 <- d.MCAR.cov.75 %>%
  select(contains(".AD.")) 

# Mean

vector.mean.MCAR.cov.75.AD.women.cl.15.25 <- AD.MCAR.cov.75[,1]
vector.mean.MCAR.cov.75.AD.men.cl.15.25 <- AD.MCAR.cov.75[,2]

vector.mean.MCAR.cov.75.AD.women.cl.25.40 <- AD.MCAR.cov.75[,3]
vector.mean.MCAR.cov.75.AD.men.cl.25.40 <- AD.MCAR.cov.75[,4]

vector.mean.MCAR.cov.75.AD.women.cl.40.50 <- AD.MCAR.cov.75[,5]
vector.mean.MCAR.cov.75.AD.men.cl.40.50 <- AD.MCAR.cov.75[,6]

# Median

vector.med.MCAR.cov.75.AD.women.cl.15.25 <- AD.MCAR.cov.75[,7]
vector.med.MCAR.cov.75.AD.men.cl.15.25 <- AD.MCAR.cov.75[,8]

vector.med.MCAR.cov.75.AD.women.cl.25.40 <- AD.MCAR.cov.75[,9]
vector.med.MCAR.cov.75.AD.men.cl.25.40 <- AD.MCAR.cov.75[,10]

vector.med.MCAR.cov.75.AD.women.cl.40.50 <- AD.MCAR.cov.75[,11]
vector.med.MCAR.cov.75.AD.men.cl.40.50 <- AD.MCAR.cov.75[,12]

# Standard deviation

vector.sd.MCAR.cov.75.AD.women.cl.15.25 <- AD.MCAR.cov.75[,13]
vector.sd.MCAR.cov.75.AD.men.cl.15.25 <- AD.MCAR.cov.75[,14]

vector.sd.MCAR.cov.75.AD.women.cl.25.40 <- AD.MCAR.cov.75[,15]
vector.sd.MCAR.cov.75.AD.men.cl.25.40 <- AD.MCAR.cov.75[,16]

vector.sd.MCAR.cov.75.AD.women.cl.40.50 <- AD.MCAR.cov.75[,17]
vector.sd.MCAR.cov.75.AD.men.cl.40.50 <- AD.MCAR.cov.75[,18]


# 80

d.MCAR.cov.80 <- d.MCAR %>%
  select(contains("cov.MCAR.80.")) # MAR - a
AD.MCAR.cov.80 <- d.MCAR.cov.80 %>%
  select(contains(".AD.")) 

# Mean

vector.mean.MCAR.cov.80.AD.women.cl.15.25 <- AD.MCAR.cov.80[,1]
vector.mean.MCAR.cov.80.AD.men.cl.15.25 <- AD.MCAR.cov.80[,2]

vector.mean.MCAR.cov.80.AD.women.cl.25.40 <- AD.MCAR.cov.80[,3]
vector.mean.MCAR.cov.80.AD.men.cl.25.40 <- AD.MCAR.cov.80[,4]

vector.mean.MCAR.cov.80.AD.women.cl.40.50 <- AD.MCAR.cov.80[,5]
vector.mean.MCAR.cov.80.AD.men.cl.40.50 <- AD.MCAR.cov.80[,6]

# Median

vector.med.MCAR.cov.80.AD.women.cl.15.25 <- AD.MCAR.cov.80[,7]
vector.med.MCAR.cov.80.AD.men.cl.15.25 <- AD.MCAR.cov.80[,8]

vector.med.MCAR.cov.80.AD.women.cl.25.40 <- AD.MCAR.cov.80[,9]
vector.med.MCAR.cov.80.AD.men.cl.25.40 <- AD.MCAR.cov.80[,10]

vector.med.MCAR.cov.80.AD.women.cl.40.50 <- AD.MCAR.cov.80[,11]
vector.med.MCAR.cov.80.AD.men.cl.40.50 <- AD.MCAR.cov.80[,12]

# Standard deviation

vector.sd.MCAR.cov.80.AD.women.cl.15.25 <- AD.MCAR.cov.80[,13]
vector.sd.MCAR.cov.80.AD.men.cl.15.25 <- AD.MCAR.cov.80[,14]

vector.sd.MCAR.cov.80.AD.women.cl.25.40 <- AD.MCAR.cov.80[,15]
vector.sd.MCAR.cov.80.AD.men.cl.25.40 <- AD.MCAR.cov.80[,16]

vector.sd.MCAR.cov.80.AD.women.cl.40.50 <- AD.MCAR.cov.80[,17]
vector.sd.MCAR.cov.80.AD.men.cl.40.50 <- AD.MCAR.cov.80[,18]


# 85

d.MCAR.cov.85 <- d.MCAR %>%
  select(contains("cov.MCAR.85.")) # MAR - a
AD.MCAR.cov.85 <- d.MCAR.cov.85 %>%
  select(contains(".AD.")) 

# Mean

vector.mean.MCAR.cov.85.AD.women.cl.15.25 <- AD.MCAR.cov.85[,1]
vector.mean.MCAR.cov.85.AD.men.cl.15.25 <- AD.MCAR.cov.85[,2]

vector.mean.MCAR.cov.85.AD.women.cl.25.40 <- AD.MCAR.cov.85[,3]
vector.mean.MCAR.cov.85.AD.men.cl.25.40 <- AD.MCAR.cov.85[,4]

vector.mean.MCAR.cov.85.AD.women.cl.40.50 <- AD.MCAR.cov.85[,5]
vector.mean.MCAR.cov.85.AD.men.cl.40.50 <- AD.MCAR.cov.85[,6]

# Median

vector.med.MCAR.cov.85.AD.women.cl.15.25 <- AD.MCAR.cov.85[,7]
vector.med.MCAR.cov.85.AD.men.cl.15.25 <- AD.MCAR.cov.85[,8]

vector.med.MCAR.cov.85.AD.women.cl.25.40 <- AD.MCAR.cov.85[,9]
vector.med.MCAR.cov.85.AD.men.cl.25.40 <- AD.MCAR.cov.85[,10]

vector.med.MCAR.cov.85.AD.women.cl.40.50 <- AD.MCAR.cov.85[,11]
vector.med.MCAR.cov.85.AD.men.cl.40.50 <- AD.MCAR.cov.85[,12]

# Standard deviation

vector.sd.MCAR.cov.85.AD.women.cl.15.25 <- AD.MCAR.cov.85[,13]
vector.sd.MCAR.cov.85.AD.men.cl.15.25 <- AD.MCAR.cov.85[,14]

vector.sd.MCAR.cov.85.AD.women.cl.25.40 <- AD.MCAR.cov.85[,15]
vector.sd.MCAR.cov.85.AD.men.cl.25.40 <- AD.MCAR.cov.85[,16]

vector.sd.MCAR.cov.85.AD.women.cl.40.50 <- AD.MCAR.cov.85[,17]
vector.sd.MCAR.cov.85.AD.men.cl.40.50 <- AD.MCAR.cov.85[,18]


# 90

d.MCAR.cov.90 <- d.MCAR %>%
  select(contains("cov.MCAR.90.")) # MAR - a
AD.MCAR.cov.90 <- d.MCAR.cov.90 %>%
  select(contains(".AD.")) 

# Mean

vector.mean.MCAR.cov.90.AD.women.cl.15.25 <- AD.MCAR.cov.90[,1]
vector.mean.MCAR.cov.90.AD.men.cl.15.25 <- AD.MCAR.cov.90[,2]

vector.mean.MCAR.cov.90.AD.women.cl.25.40 <- AD.MCAR.cov.90[,3]
vector.mean.MCAR.cov.90.AD.men.cl.25.40 <- AD.MCAR.cov.90[,4]

vector.mean.MCAR.cov.90.AD.women.cl.40.50 <- AD.MCAR.cov.90[,5]
vector.mean.MCAR.cov.90.AD.men.cl.40.50 <- AD.MCAR.cov.90[,6]

# Median

vector.med.MCAR.cov.90.AD.women.cl.15.25 <- AD.MCAR.cov.90[,7]
vector.med.MCAR.cov.90.AD.men.cl.15.25 <- AD.MCAR.cov.90[,8]

vector.med.MCAR.cov.90.AD.women.cl.25.40 <- AD.MCAR.cov.90[,9]
vector.med.MCAR.cov.90.AD.men.cl.25.40 <- AD.MCAR.cov.90[,10]

vector.med.MCAR.cov.90.AD.women.cl.40.50 <- AD.MCAR.cov.90[,11]
vector.med.MCAR.cov.90.AD.men.cl.40.50 <- AD.MCAR.cov.90[,12]

# Standard deviation

vector.sd.MCAR.cov.90.AD.women.cl.15.25 <- AD.MCAR.cov.90[,13]
vector.sd.MCAR.cov.90.AD.men.cl.15.25 <- AD.MCAR.cov.90[,14]

vector.sd.MCAR.cov.90.AD.women.cl.25.40 <- AD.MCAR.cov.90[,15]
vector.sd.MCAR.cov.90.AD.men.cl.25.40 <- AD.MCAR.cov.90[,16]

vector.sd.MCAR.cov.90.AD.women.cl.40.50 <- AD.MCAR.cov.90[,17]
vector.sd.MCAR.cov.90.AD.men.cl.40.50 <- AD.MCAR.cov.90[,18]


# 95

d.MCAR.cov.95 <- d.MCAR %>%
  select(contains("cov.MCAR.95.")) # MAR - a
AD.MCAR.cov.95 <- d.MCAR.cov.95 %>%
  select(contains(".AD.")) 

# Mean

vector.mean.MCAR.cov.95.AD.women.cl.15.25 <- AD.MCAR.cov.95[,1]
vector.mean.MCAR.cov.95.AD.men.cl.15.25 <- AD.MCAR.cov.95[,2]

vector.mean.MCAR.cov.95.AD.women.cl.25.40 <- AD.MCAR.cov.95[,3]
vector.mean.MCAR.cov.95.AD.men.cl.25.40 <- AD.MCAR.cov.95[,4]

vector.mean.MCAR.cov.95.AD.women.cl.40.50 <- AD.MCAR.cov.95[,5]
vector.mean.MCAR.cov.95.AD.men.cl.40.50 <- AD.MCAR.cov.95[,6]

# Median

vector.med.MCAR.cov.95.AD.women.cl.15.25 <- AD.MCAR.cov.95[,7]
vector.med.MCAR.cov.95.AD.men.cl.15.25 <- AD.MCAR.cov.95[,8]

vector.med.MCAR.cov.95.AD.women.cl.25.40 <- AD.MCAR.cov.95[,9]
vector.med.MCAR.cov.95.AD.men.cl.25.40 <- AD.MCAR.cov.95[,10]

vector.med.MCAR.cov.95.AD.women.cl.40.50 <- AD.MCAR.cov.95[,11]
vector.med.MCAR.cov.95.AD.men.cl.40.50 <- AD.MCAR.cov.95[,12]

# Standard deviation

vector.sd.MCAR.cov.95.AD.women.cl.15.25 <- AD.MCAR.cov.95[,13]
vector.sd.MCAR.cov.95.AD.men.cl.15.25 <- AD.MCAR.cov.95[,14]

vector.sd.MCAR.cov.95.AD.women.cl.25.40 <- AD.MCAR.cov.95[,15]
vector.sd.MCAR.cov.95.AD.men.cl.25.40 <- AD.MCAR.cov.95[,16]

vector.sd.MCAR.cov.95.AD.women.cl.40.50 <- AD.MCAR.cov.95[,17]
vector.sd.MCAR.cov.95.AD.men.cl.40.50 <- AD.MCAR.cov.95[,18]






# MAR - a

d.MAR <- dr %>%
  select(contains("MAR."))

# 35

d.MAR.a.cov.35 <- d.MAR %>%
  select(contains("cov.MAR.a.35.")) # MAR - a
AD.MAR.a.cov.35 <- d.MAR.a.cov.35 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.a.cov.35.AD.women.cl.15.25 <- AD.MAR.a.cov.35[,1]
vector.mean.MAR.a.cov.35.AD.men.cl.15.25 <- AD.MAR.a.cov.35[,2]

vector.mean.MAR.a.cov.35.AD.women.cl.25.40 <- AD.MAR.a.cov.35[,3]
vector.mean.MAR.a.cov.35.AD.men.cl.25.40 <- AD.MAR.a.cov.35[,4]

vector.mean.MAR.a.cov.35.AD.women.cl.40.50 <- AD.MAR.a.cov.35[,5]
vector.mean.MAR.a.cov.35.AD.men.cl.40.50 <- AD.MAR.a.cov.35[,6]

# Median

vector.med.MAR.a.cov.35.AD.women.cl.15.25 <- AD.MAR.a.cov.35[,7]
vector.med.MAR.a.cov.35.AD.men.cl.15.25 <- AD.MAR.a.cov.35[,8]

vector.med.MAR.a.cov.35.AD.women.cl.25.40 <- AD.MAR.a.cov.35[,9]
vector.med.MAR.a.cov.35.AD.men.cl.25.40 <- AD.MAR.a.cov.35[,10]

vector.med.MAR.a.cov.35.AD.women.cl.40.50 <- AD.MAR.a.cov.35[,11]
vector.med.MAR.a.cov.35.AD.men.cl.40.50 <- AD.MAR.a.cov.35[,12]

# Standard deviation

vector.sd.MAR.a.cov.35.AD.women.cl.15.25 <- AD.MAR.a.cov.35[,13]
vector.sd.MAR.a.cov.35.AD.men.cl.15.25 <- AD.MAR.a.cov.35[,14]

vector.sd.MAR.a.cov.35.AD.women.cl.25.40 <- AD.MAR.a.cov.35[,15]
vector.sd.MAR.a.cov.35.AD.men.cl.25.40 <- AD.MAR.a.cov.35[,16]

vector.sd.MAR.a.cov.35.AD.women.cl.40.50 <- AD.MAR.a.cov.35[,17]
vector.sd.MAR.a.cov.35.AD.men.cl.40.50 <- AD.MAR.a.cov.35[,18]


# 40



d.MAR.a.cov.40 <- d.MAR %>%
  select(contains("cov.MAR.a.40.")) # MAR - a
AD.MAR.a.cov.40 <- d.MAR.a.cov.40 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.a.cov.40.AD.women.cl.15.25 <- AD.MAR.a.cov.40[,1]
vector.mean.MAR.a.cov.40.AD.men.cl.15.25 <- AD.MAR.a.cov.40[,2]

vector.mean.MAR.a.cov.40.AD.women.cl.25.40 <- AD.MAR.a.cov.40[,3]
vector.mean.MAR.a.cov.40.AD.men.cl.25.40 <- AD.MAR.a.cov.40[,4]

vector.mean.MAR.a.cov.40.AD.women.cl.40.50 <- AD.MAR.a.cov.40[,5]
vector.mean.MAR.a.cov.40.AD.men.cl.40.50 <- AD.MAR.a.cov.40[,6]

# Median

vector.med.MAR.a.cov.40.AD.women.cl.15.25 <- AD.MAR.a.cov.40[,7]
vector.med.MAR.a.cov.40.AD.men.cl.15.25 <- AD.MAR.a.cov.40[,8]

vector.med.MAR.a.cov.40.AD.women.cl.25.40 <- AD.MAR.a.cov.40[,9]
vector.med.MAR.a.cov.40.AD.men.cl.25.40 <- AD.MAR.a.cov.40[,10]

vector.med.MAR.a.cov.40.AD.women.cl.40.50 <- AD.MAR.a.cov.40[,11]
vector.med.MAR.a.cov.40.AD.men.cl.40.50 <- AD.MAR.a.cov.40[,12]

# Standard deviation

vector.sd.MAR.a.cov.40.AD.women.cl.15.25 <- AD.MAR.a.cov.40[,13]
vector.sd.MAR.a.cov.40.AD.men.cl.15.25 <- AD.MAR.a.cov.40[,14]

vector.sd.MAR.a.cov.40.AD.women.cl.25.40 <- AD.MAR.a.cov.40[,15]
vector.sd.MAR.a.cov.40.AD.men.cl.25.40 <- AD.MAR.a.cov.40[,16]

vector.sd.MAR.a.cov.40.AD.women.cl.40.50 <- AD.MAR.a.cov.40[,17]
vector.sd.MAR.a.cov.40.AD.men.cl.40.50 <- AD.MAR.a.cov.40[,18]



# 45



d.MAR.a.cov.45 <- d.MAR %>%
  select(contains("cov.MAR.a.45."))
AD.MAR.a.cov.45 <- d.MAR.a.cov.45 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.a.cov.45.AD.women.cl.15.25 <- AD.MAR.a.cov.45[,1]
vector.mean.MAR.a.cov.45.AD.men.cl.15.25 <- AD.MAR.a.cov.45[,2]

vector.mean.MAR.a.cov.45.AD.women.cl.25.40 <- AD.MAR.a.cov.45[,3]
vector.mean.MAR.a.cov.45.AD.men.cl.25.40 <- AD.MAR.a.cov.45[,4]

vector.mean.MAR.a.cov.45.AD.women.cl.40.50 <- AD.MAR.a.cov.45[,5]
vector.mean.MAR.a.cov.45.AD.men.cl.40.50 <- AD.MAR.a.cov.45[,6]

# Median

vector.med.MAR.a.cov.45.AD.women.cl.15.25 <- AD.MAR.a.cov.45[,7]
vector.med.MAR.a.cov.45.AD.men.cl.15.25 <- AD.MAR.a.cov.45[,8]

vector.med.MAR.a.cov.45.AD.women.cl.25.40 <- AD.MAR.a.cov.45[,9]
vector.med.MAR.a.cov.45.AD.men.cl.25.40 <- AD.MAR.a.cov.45[,10]

vector.med.MAR.a.cov.45.AD.women.cl.40.50 <- AD.MAR.a.cov.45[,11]
vector.med.MAR.a.cov.45.AD.men.cl.40.50 <- AD.MAR.a.cov.45[,12]

# Standard deviation

vector.sd.MAR.a.cov.45.AD.women.cl.15.25 <- AD.MAR.a.cov.45[,13]
vector.sd.MAR.a.cov.45.AD.men.cl.15.25 <- AD.MAR.a.cov.45[,14]

vector.sd.MAR.a.cov.45.AD.women.cl.25.40 <- AD.MAR.a.cov.45[,15]
vector.sd.MAR.a.cov.45.AD.men.cl.25.40 <- AD.MAR.a.cov.45[,16]

vector.sd.MAR.a.cov.45.AD.women.cl.40.50 <- AD.MAR.a.cov.45[,17]
vector.sd.MAR.a.cov.45.AD.men.cl.40.50 <- AD.MAR.a.cov.45[,18]



# 50



d.MAR.a.cov.50 <- d.MAR %>%
  select(contains("cov.MAR.a.50."))
AD.MAR.a.cov.50 <- d.MAR.a.cov.50 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.a.cov.50.AD.women.cl.15.25 <- AD.MAR.a.cov.50[,1]
vector.mean.MAR.a.cov.50.AD.men.cl.15.25 <- AD.MAR.a.cov.50[,2]

vector.mean.MAR.a.cov.50.AD.women.cl.25.40 <- AD.MAR.a.cov.50[,3]
vector.mean.MAR.a.cov.50.AD.men.cl.25.40 <- AD.MAR.a.cov.50[,4]

vector.mean.MAR.a.cov.50.AD.women.cl.40.50 <- AD.MAR.a.cov.50[,5]
vector.mean.MAR.a.cov.50.AD.men.cl.40.50 <- AD.MAR.a.cov.50[,6]

# Median

vector.med.MAR.a.cov.50.AD.women.cl.15.25 <- AD.MAR.a.cov.50[,7]
vector.med.MAR.a.cov.50.AD.men.cl.15.25 <- AD.MAR.a.cov.50[,8]

vector.med.MAR.a.cov.50.AD.women.cl.25.40 <- AD.MAR.a.cov.50[,9]
vector.med.MAR.a.cov.50.AD.men.cl.25.40 <- AD.MAR.a.cov.50[,10]

vector.med.MAR.a.cov.50.AD.women.cl.40.50 <- AD.MAR.a.cov.50[,11]
vector.med.MAR.a.cov.50.AD.men.cl.40.50 <- AD.MAR.a.cov.50[,12]

# Standard deviation

vector.sd.MAR.a.cov.50.AD.women.cl.15.25 <- AD.MAR.a.cov.50[,13]
vector.sd.MAR.a.cov.50.AD.men.cl.15.25 <- AD.MAR.a.cov.50[,14]

vector.sd.MAR.a.cov.50.AD.women.cl.25.40 <- AD.MAR.a.cov.50[,15]
vector.sd.MAR.a.cov.50.AD.men.cl.25.40 <- AD.MAR.a.cov.50[,16]

vector.sd.MAR.a.cov.50.AD.women.cl.40.50 <- AD.MAR.a.cov.50[,17]
vector.sd.MAR.a.cov.50.AD.men.cl.40.50 <- AD.MAR.a.cov.50[,18]



# 55



d.MAR.a.cov.55 <- d.MAR %>%
  select(contains("cov.MAR.a.55."))
AD.MAR.a.cov.55 <- d.MAR.a.cov.55 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.a.cov.55.AD.women.cl.15.25 <- AD.MAR.a.cov.55[,1]
vector.mean.MAR.a.cov.55.AD.men.cl.15.25 <- AD.MAR.a.cov.55[,2]

vector.mean.MAR.a.cov.55.AD.women.cl.25.40 <- AD.MAR.a.cov.55[,3]
vector.mean.MAR.a.cov.55.AD.men.cl.25.40 <- AD.MAR.a.cov.55[,4]

vector.mean.MAR.a.cov.55.AD.women.cl.40.50 <- AD.MAR.a.cov.55[,5]
vector.mean.MAR.a.cov.55.AD.men.cl.40.50 <- AD.MAR.a.cov.55[,6]

# Median

vector.med.MAR.a.cov.55.AD.women.cl.15.25 <- AD.MAR.a.cov.55[,7]
vector.med.MAR.a.cov.55.AD.men.cl.15.25 <- AD.MAR.a.cov.55[,8]

vector.med.MAR.a.cov.55.AD.women.cl.25.40 <- AD.MAR.a.cov.55[,9]
vector.med.MAR.a.cov.55.AD.men.cl.25.40 <- AD.MAR.a.cov.55[,10]

vector.med.MAR.a.cov.55.AD.women.cl.40.50 <- AD.MAR.a.cov.55[,11]
vector.med.MAR.a.cov.55.AD.men.cl.40.50 <- AD.MAR.a.cov.55[,12]

# Standard deviation

vector.sd.MAR.a.cov.55.AD.women.cl.15.25 <- AD.MAR.a.cov.55[,13]
vector.sd.MAR.a.cov.55.AD.men.cl.15.25 <- AD.MAR.a.cov.55[,14]

vector.sd.MAR.a.cov.55.AD.women.cl.25.40 <- AD.MAR.a.cov.55[,15]
vector.sd.MAR.a.cov.55.AD.men.cl.25.40 <- AD.MAR.a.cov.55[,16]

vector.sd.MAR.a.cov.55.AD.women.cl.40.50 <- AD.MAR.a.cov.55[,17]
vector.sd.MAR.a.cov.55.AD.men.cl.40.50 <- AD.MAR.a.cov.55[,18]



# 60



d.MAR.a.cov.60 <- d.MAR %>%
  select(contains("cov.MAR.a.60."))
AD.MAR.a.cov.60 <- d.MAR.a.cov.60 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.a.cov.60.AD.women.cl.15.25 <- AD.MAR.a.cov.60[,1]
vector.mean.MAR.a.cov.60.AD.men.cl.15.25 <- AD.MAR.a.cov.60[,2]

vector.mean.MAR.a.cov.60.AD.women.cl.25.40 <- AD.MAR.a.cov.60[,3]
vector.mean.MAR.a.cov.60.AD.men.cl.25.40 <- AD.MAR.a.cov.60[,4]

vector.mean.MAR.a.cov.60.AD.women.cl.40.50 <- AD.MAR.a.cov.60[,5]
vector.mean.MAR.a.cov.60.AD.men.cl.40.50 <- AD.MAR.a.cov.60[,6]

# Median

vector.med.MAR.a.cov.60.AD.women.cl.15.25 <- AD.MAR.a.cov.60[,7]
vector.med.MAR.a.cov.60.AD.men.cl.15.25 <- AD.MAR.a.cov.60[,8]

vector.med.MAR.a.cov.60.AD.women.cl.25.40 <- AD.MAR.a.cov.60[,9]
vector.med.MAR.a.cov.60.AD.men.cl.25.40 <- AD.MAR.a.cov.60[,10]

vector.med.MAR.a.cov.60.AD.women.cl.40.50 <- AD.MAR.a.cov.60[,11]
vector.med.MAR.a.cov.60.AD.men.cl.40.50 <- AD.MAR.a.cov.60[,12]

# Standard deviation

vector.sd.MAR.a.cov.60.AD.women.cl.15.25 <- AD.MAR.a.cov.60[,13]
vector.sd.MAR.a.cov.60.AD.men.cl.15.25 <- AD.MAR.a.cov.60[,14]

vector.sd.MAR.a.cov.60.AD.women.cl.25.40 <- AD.MAR.a.cov.60[,15]
vector.sd.MAR.a.cov.60.AD.men.cl.25.40 <- AD.MAR.a.cov.60[,16]

vector.sd.MAR.a.cov.60.AD.women.cl.40.50 <- AD.MAR.a.cov.60[,17]
vector.sd.MAR.a.cov.60.AD.men.cl.40.50 <- AD.MAR.a.cov.60[,18]



# 65



d.MAR.a.cov.65 <- d.MAR %>%
  select(contains("cov.MAR.a.65."))
AD.MAR.a.cov.65 <- d.MAR.a.cov.65 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.a.cov.65.AD.women.cl.15.25 <- AD.MAR.a.cov.65[,1]
vector.mean.MAR.a.cov.65.AD.men.cl.15.25 <- AD.MAR.a.cov.65[,2]

vector.mean.MAR.a.cov.65.AD.women.cl.25.40 <- AD.MAR.a.cov.65[,3]
vector.mean.MAR.a.cov.65.AD.men.cl.25.40 <- AD.MAR.a.cov.65[,4]

vector.mean.MAR.a.cov.65.AD.women.cl.40.50 <- AD.MAR.a.cov.65[,5]
vector.mean.MAR.a.cov.65.AD.men.cl.40.50 <- AD.MAR.a.cov.65[,6]

# Median

vector.med.MAR.a.cov.65.AD.women.cl.15.25 <- AD.MAR.a.cov.65[,7]
vector.med.MAR.a.cov.65.AD.men.cl.15.25 <- AD.MAR.a.cov.65[,8]

vector.med.MAR.a.cov.65.AD.women.cl.25.40 <- AD.MAR.a.cov.65[,9]
vector.med.MAR.a.cov.65.AD.men.cl.25.40 <- AD.MAR.a.cov.65[,10]

vector.med.MAR.a.cov.65.AD.women.cl.40.50 <- AD.MAR.a.cov.65[,11]
vector.med.MAR.a.cov.65.AD.men.cl.40.50 <- AD.MAR.a.cov.65[,12]

# Standard deviation

vector.sd.MAR.a.cov.65.AD.women.cl.15.25 <- AD.MAR.a.cov.65[,13]
vector.sd.MAR.a.cov.65.AD.men.cl.15.25 <- AD.MAR.a.cov.65[,14]

vector.sd.MAR.a.cov.65.AD.women.cl.25.40 <- AD.MAR.a.cov.65[,15]
vector.sd.MAR.a.cov.65.AD.men.cl.25.40 <- AD.MAR.a.cov.65[,16]

vector.sd.MAR.a.cov.65.AD.women.cl.40.50 <- AD.MAR.a.cov.65[,17]
vector.sd.MAR.a.cov.65.AD.men.cl.40.50 <- AD.MAR.a.cov.65[,18]



# 70



d.MAR.a.cov.70 <- d.MAR %>%
  select(contains("cov.MAR.a.70."))
AD.MAR.a.cov.70 <- d.MAR.a.cov.70 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.a.cov.70.AD.women.cl.15.25 <- AD.MAR.a.cov.70[,1]
vector.mean.MAR.a.cov.70.AD.men.cl.15.25 <- AD.MAR.a.cov.70[,2]

vector.mean.MAR.a.cov.70.AD.women.cl.25.40 <- AD.MAR.a.cov.70[,3]
vector.mean.MAR.a.cov.70.AD.men.cl.25.40 <- AD.MAR.a.cov.70[,4]

vector.mean.MAR.a.cov.70.AD.women.cl.40.50 <- AD.MAR.a.cov.70[,5]
vector.mean.MAR.a.cov.70.AD.men.cl.40.50 <- AD.MAR.a.cov.70[,6]

# Median

vector.med.MAR.a.cov.70.AD.women.cl.15.25 <- AD.MAR.a.cov.70[,7]
vector.med.MAR.a.cov.70.AD.men.cl.15.25 <- AD.MAR.a.cov.70[,8]

vector.med.MAR.a.cov.70.AD.women.cl.25.40 <- AD.MAR.a.cov.70[,9]
vector.med.MAR.a.cov.70.AD.men.cl.25.40 <- AD.MAR.a.cov.70[,10]

vector.med.MAR.a.cov.70.AD.women.cl.40.50 <- AD.MAR.a.cov.70[,11]
vector.med.MAR.a.cov.70.AD.men.cl.40.50 <- AD.MAR.a.cov.70[,12]

# Standard deviation

vector.sd.MAR.a.cov.70.AD.women.cl.15.25 <- AD.MAR.a.cov.70[,13]
vector.sd.MAR.a.cov.70.AD.men.cl.15.25 <- AD.MAR.a.cov.70[,14]

vector.sd.MAR.a.cov.70.AD.women.cl.25.40 <- AD.MAR.a.cov.70[,15]
vector.sd.MAR.a.cov.70.AD.men.cl.25.40 <- AD.MAR.a.cov.70[,16]

vector.sd.MAR.a.cov.70.AD.women.cl.40.50 <- AD.MAR.a.cov.70[,17]
vector.sd.MAR.a.cov.70.AD.men.cl.40.50 <- AD.MAR.a.cov.70[,18]



# 75




d.MAR.a.cov.75 <- d.MAR %>%
  select(contains("cov.MAR.a.75."))
AD.MAR.a.cov.75 <- d.MAR.a.cov.75 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.a.cov.75.AD.women.cl.15.25 <- AD.MAR.a.cov.75[,1]
vector.mean.MAR.a.cov.75.AD.men.cl.15.25 <- AD.MAR.a.cov.75[,2]

vector.mean.MAR.a.cov.75.AD.women.cl.25.40 <- AD.MAR.a.cov.75[,3]
vector.mean.MAR.a.cov.75.AD.men.cl.25.40 <- AD.MAR.a.cov.75[,4]

vector.mean.MAR.a.cov.75.AD.women.cl.40.50 <- AD.MAR.a.cov.75[,5]
vector.mean.MAR.a.cov.75.AD.men.cl.40.50 <- AD.MAR.a.cov.75[,6]

# Median

vector.med.MAR.a.cov.75.AD.women.cl.15.25 <- AD.MAR.a.cov.75[,7]
vector.med.MAR.a.cov.75.AD.men.cl.15.25 <- AD.MAR.a.cov.75[,8]

vector.med.MAR.a.cov.75.AD.women.cl.25.40 <- AD.MAR.a.cov.75[,9]
vector.med.MAR.a.cov.75.AD.men.cl.25.40 <- AD.MAR.a.cov.75[,10]

vector.med.MAR.a.cov.75.AD.women.cl.40.50 <- AD.MAR.a.cov.75[,11]
vector.med.MAR.a.cov.75.AD.men.cl.40.50 <- AD.MAR.a.cov.75[,12]

# Standard deviation

vector.sd.MAR.a.cov.75.AD.women.cl.15.25 <- AD.MAR.a.cov.75[,13]
vector.sd.MAR.a.cov.75.AD.men.cl.15.25 <- AD.MAR.a.cov.75[,14]

vector.sd.MAR.a.cov.75.AD.women.cl.25.40 <- AD.MAR.a.cov.75[,15]
vector.sd.MAR.a.cov.75.AD.men.cl.25.40 <- AD.MAR.a.cov.75[,16]

vector.sd.MAR.a.cov.75.AD.women.cl.40.50 <- AD.MAR.a.cov.75[,17]
vector.sd.MAR.a.cov.75.AD.men.cl.40.50 <- AD.MAR.a.cov.75[,18]



# 80



d.MAR.a.cov.80 <- d.MAR %>%
  select(contains("cov.MAR.a.80."))
AD.MAR.a.cov.80 <- d.MAR.a.cov.80 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.a.cov.80.AD.women.cl.15.25 <- AD.MAR.a.cov.80[,1]
vector.mean.MAR.a.cov.80.AD.men.cl.15.25 <- AD.MAR.a.cov.80[,2]

vector.mean.MAR.a.cov.80.AD.women.cl.25.40 <- AD.MAR.a.cov.80[,3]
vector.mean.MAR.a.cov.80.AD.men.cl.25.40 <- AD.MAR.a.cov.80[,4]

vector.mean.MAR.a.cov.80.AD.women.cl.40.50 <- AD.MAR.a.cov.80[,5]
vector.mean.MAR.a.cov.80.AD.men.cl.40.50 <- AD.MAR.a.cov.80[,6]

# Median

vector.med.MAR.a.cov.80.AD.women.cl.15.25 <- AD.MAR.a.cov.80[,7]
vector.med.MAR.a.cov.80.AD.men.cl.15.25 <- AD.MAR.a.cov.80[,8]

vector.med.MAR.a.cov.80.AD.women.cl.25.40 <- AD.MAR.a.cov.80[,9]
vector.med.MAR.a.cov.80.AD.men.cl.25.40 <- AD.MAR.a.cov.80[,10]

vector.med.MAR.a.cov.80.AD.women.cl.40.50 <- AD.MAR.a.cov.80[,11]
vector.med.MAR.a.cov.80.AD.men.cl.40.50 <- AD.MAR.a.cov.80[,12]

# Standard deviation

vector.sd.MAR.a.cov.80.AD.women.cl.15.25 <- AD.MAR.a.cov.80[,13]
vector.sd.MAR.a.cov.80.AD.men.cl.15.25 <- AD.MAR.a.cov.80[,14]

vector.sd.MAR.a.cov.80.AD.women.cl.25.40 <- AD.MAR.a.cov.80[,15]
vector.sd.MAR.a.cov.80.AD.men.cl.25.40 <- AD.MAR.a.cov.80[,16]

vector.sd.MAR.a.cov.80.AD.women.cl.40.50 <- AD.MAR.a.cov.80[,17]
vector.sd.MAR.a.cov.80.AD.men.cl.40.50 <- AD.MAR.a.cov.80[,18]



# 85



d.MAR.a.cov.85 <- d.MAR %>%
  select(contains("cov.MAR.a.85."))
AD.MAR.a.cov.85 <- d.MAR.a.cov.85 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.a.cov.85.AD.women.cl.15.25 <- AD.MAR.a.cov.85[,1]
vector.mean.MAR.a.cov.85.AD.men.cl.15.25 <- AD.MAR.a.cov.85[,2]

vector.mean.MAR.a.cov.85.AD.women.cl.25.40 <- AD.MAR.a.cov.85[,3]
vector.mean.MAR.a.cov.85.AD.men.cl.25.40 <- AD.MAR.a.cov.85[,4]

vector.mean.MAR.a.cov.85.AD.women.cl.40.50 <- AD.MAR.a.cov.85[,5]
vector.mean.MAR.a.cov.85.AD.men.cl.40.50 <- AD.MAR.a.cov.85[,6]

# Median

vector.med.MAR.a.cov.85.AD.women.cl.15.25 <- AD.MAR.a.cov.85[,7]
vector.med.MAR.a.cov.85.AD.men.cl.15.25 <- AD.MAR.a.cov.85[,8]

vector.med.MAR.a.cov.85.AD.women.cl.25.40 <- AD.MAR.a.cov.85[,9]
vector.med.MAR.a.cov.85.AD.men.cl.25.40 <- AD.MAR.a.cov.85[,10]

vector.med.MAR.a.cov.85.AD.women.cl.40.50 <- AD.MAR.a.cov.85[,11]
vector.med.MAR.a.cov.85.AD.men.cl.40.50 <- AD.MAR.a.cov.85[,12]

# Standard deviation

vector.sd.MAR.a.cov.85.AD.women.cl.15.25 <- AD.MAR.a.cov.85[,13]
vector.sd.MAR.a.cov.85.AD.men.cl.15.25 <- AD.MAR.a.cov.85[,14]

vector.sd.MAR.a.cov.85.AD.women.cl.25.40 <- AD.MAR.a.cov.85[,15]
vector.sd.MAR.a.cov.85.AD.men.cl.25.40 <- AD.MAR.a.cov.85[,16]

vector.sd.MAR.a.cov.85.AD.women.cl.40.50 <- AD.MAR.a.cov.85[,17]
vector.sd.MAR.a.cov.85.AD.men.cl.40.50 <- AD.MAR.a.cov.85[,18]



# 90



d.MAR.a.cov.90 <- d.MAR %>%
  select(contains("cov.MAR.a.90."))
AD.MAR.a.cov.90 <- d.MAR.a.cov.90 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.a.cov.90.AD.women.cl.15.25 <- AD.MAR.a.cov.90[,1]
vector.mean.MAR.a.cov.90.AD.men.cl.15.25 <- AD.MAR.a.cov.90[,2]

vector.mean.MAR.a.cov.90.AD.women.cl.25.40 <- AD.MAR.a.cov.90[,3]
vector.mean.MAR.a.cov.90.AD.men.cl.25.40 <- AD.MAR.a.cov.90[,4]

vector.mean.MAR.a.cov.90.AD.women.cl.40.50 <- AD.MAR.a.cov.90[,5]
vector.mean.MAR.a.cov.90.AD.men.cl.40.50 <- AD.MAR.a.cov.90[,6]

# Median

vector.med.MAR.a.cov.90.AD.women.cl.15.25 <- AD.MAR.a.cov.90[,7]
vector.med.MAR.a.cov.90.AD.men.cl.15.25 <- AD.MAR.a.cov.90[,8]

vector.med.MAR.a.cov.90.AD.women.cl.25.40 <- AD.MAR.a.cov.90[,9]
vector.med.MAR.a.cov.90.AD.men.cl.25.40 <- AD.MAR.a.cov.90[,10]

vector.med.MAR.a.cov.90.AD.women.cl.40.50 <- AD.MAR.a.cov.90[,11]
vector.med.MAR.a.cov.90.AD.men.cl.40.50 <- AD.MAR.a.cov.90[,12]

# Standard deviation

vector.sd.MAR.a.cov.90.AD.women.cl.15.25 <- AD.MAR.a.cov.90[,13]
vector.sd.MAR.a.cov.90.AD.men.cl.15.25 <- AD.MAR.a.cov.90[,14]

vector.sd.MAR.a.cov.90.AD.women.cl.25.40 <- AD.MAR.a.cov.90[,15]
vector.sd.MAR.a.cov.90.AD.men.cl.25.40 <- AD.MAR.a.cov.90[,16]

vector.sd.MAR.a.cov.90.AD.women.cl.40.50 <- AD.MAR.a.cov.90[,17]
vector.sd.MAR.a.cov.90.AD.men.cl.40.50 <- AD.MAR.a.cov.90[,18]



# 95



d.MAR.a.cov.95 <- d.MAR %>%
  select(contains("cov.MAR.a.95."))
AD.MAR.a.cov.95 <- d.MAR.a.cov.95 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.a.cov.95.AD.women.cl.15.25 <- AD.MAR.a.cov.95[,1]
vector.mean.MAR.a.cov.95.AD.men.cl.15.25 <- AD.MAR.a.cov.95[,2]

vector.mean.MAR.a.cov.95.AD.women.cl.25.40 <- AD.MAR.a.cov.95[,3]
vector.mean.MAR.a.cov.95.AD.men.cl.25.40 <- AD.MAR.a.cov.95[,4]

vector.mean.MAR.a.cov.95.AD.women.cl.40.50 <- AD.MAR.a.cov.95[,5]
vector.mean.MAR.a.cov.95.AD.men.cl.40.50 <- AD.MAR.a.cov.95[,6]

# Median

vector.med.MAR.a.cov.95.AD.women.cl.15.25 <- AD.MAR.a.cov.95[,7]
vector.med.MAR.a.cov.95.AD.men.cl.15.25 <- AD.MAR.a.cov.95[,8]

vector.med.MAR.a.cov.95.AD.women.cl.25.40 <- AD.MAR.a.cov.95[,9]
vector.med.MAR.a.cov.95.AD.men.cl.25.40 <- AD.MAR.a.cov.95[,10]

vector.med.MAR.a.cov.95.AD.women.cl.40.50 <- AD.MAR.a.cov.95[,11]
vector.med.MAR.a.cov.95.AD.men.cl.40.50 <- AD.MAR.a.cov.95[,12]

# Standard deviation

vector.sd.MAR.a.cov.95.AD.women.cl.15.25 <- AD.MAR.a.cov.95[,13]
vector.sd.MAR.a.cov.95.AD.men.cl.15.25 <- AD.MAR.a.cov.95[,14]

vector.sd.MAR.a.cov.95.AD.women.cl.25.40 <- AD.MAR.a.cov.95[,15]
vector.sd.MAR.a.cov.95.AD.men.cl.25.40 <- AD.MAR.a.cov.95[,16]

vector.sd.MAR.a.cov.95.AD.women.cl.40.50 <- AD.MAR.a.cov.95[,17]
vector.sd.MAR.a.cov.95.AD.men.cl.40.50 <- AD.MAR.a.cov.95[,18]



# MAR - b


# 35

d.MAR.b.cov.35 <- d.MAR %>%
  select(contains("cov.MAR.b.35.")) # MAR - a
AD.MAR.b.cov.35 <- d.MAR.b.cov.35 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.b.cov.35.AD.women.cl.15.25 <- AD.MAR.b.cov.35[,1]
vector.mean.MAR.b.cov.35.AD.men.cl.15.25 <- AD.MAR.b.cov.35[,2]

vector.mean.MAR.b.cov.35.AD.women.cl.25.40 <- AD.MAR.b.cov.35[,3]
vector.mean.MAR.b.cov.35.AD.men.cl.25.40 <- AD.MAR.b.cov.35[,4]

vector.mean.MAR.b.cov.35.AD.women.cl.40.50 <- AD.MAR.b.cov.35[,5]
vector.mean.MAR.b.cov.35.AD.men.cl.40.50 <- AD.MAR.b.cov.35[,6]

# Median

vector.med.MAR.b.cov.35.AD.women.cl.15.25 <- AD.MAR.b.cov.35[,7]
vector.med.MAR.b.cov.35.AD.men.cl.15.25 <- AD.MAR.b.cov.35[,8]

vector.med.MAR.b.cov.35.AD.women.cl.25.40 <- AD.MAR.b.cov.35[,9]
vector.med.MAR.b.cov.35.AD.men.cl.25.40 <- AD.MAR.b.cov.35[,10]

vector.med.MAR.b.cov.35.AD.women.cl.40.50 <- AD.MAR.b.cov.35[,11]
vector.med.MAR.b.cov.35.AD.men.cl.40.50 <- AD.MAR.b.cov.35[,12]

# Standard deviation

vector.sd.MAR.b.cov.35.AD.women.cl.15.25 <- AD.MAR.b.cov.35[,13]
vector.sd.MAR.b.cov.35.AD.men.cl.15.25 <- AD.MAR.b.cov.35[,14]

vector.sd.MAR.b.cov.35.AD.women.cl.25.40 <- AD.MAR.b.cov.35[,15]
vector.sd.MAR.b.cov.35.AD.men.cl.25.40 <- AD.MAR.b.cov.35[,16]

vector.sd.MAR.b.cov.35.AD.women.cl.40.50 <- AD.MAR.b.cov.35[,17]
vector.sd.MAR.b.cov.35.AD.men.cl.40.50 <- AD.MAR.b.cov.35[,18]


# 40



d.MAR.b.cov.40 <- d.MAR %>%
  select(contains("cov.MAR.b.40.")) # MAR - a
AD.MAR.b.cov.40 <- d.MAR.b.cov.40 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.b.cov.40.AD.women.cl.15.25 <- AD.MAR.b.cov.40[,1]
vector.mean.MAR.b.cov.40.AD.men.cl.15.25 <- AD.MAR.b.cov.40[,2]

vector.mean.MAR.b.cov.40.AD.women.cl.25.40 <- AD.MAR.b.cov.40[,3]
vector.mean.MAR.b.cov.40.AD.men.cl.25.40 <- AD.MAR.b.cov.40[,4]

vector.mean.MAR.b.cov.40.AD.women.cl.40.50 <- AD.MAR.b.cov.40[,5]
vector.mean.MAR.b.cov.40.AD.men.cl.40.50 <- AD.MAR.b.cov.40[,6]

# Median

vector.med.MAR.b.cov.40.AD.women.cl.15.25 <- AD.MAR.b.cov.40[,7]
vector.med.MAR.b.cov.40.AD.men.cl.15.25 <- AD.MAR.b.cov.40[,8]

vector.med.MAR.b.cov.40.AD.women.cl.25.40 <- AD.MAR.b.cov.40[,9]
vector.med.MAR.b.cov.40.AD.men.cl.25.40 <- AD.MAR.b.cov.40[,10]

vector.med.MAR.b.cov.40.AD.women.cl.40.50 <- AD.MAR.b.cov.40[,11]
vector.med.MAR.b.cov.40.AD.men.cl.40.50 <- AD.MAR.b.cov.40[,12]

# Standard deviation

vector.sd.MAR.b.cov.40.AD.women.cl.15.25 <- AD.MAR.b.cov.40[,13]
vector.sd.MAR.b.cov.40.AD.men.cl.15.25 <- AD.MAR.b.cov.40[,14]

vector.sd.MAR.b.cov.40.AD.women.cl.25.40 <- AD.MAR.b.cov.40[,15]
vector.sd.MAR.b.cov.40.AD.men.cl.25.40 <- AD.MAR.b.cov.40[,16]

vector.sd.MAR.b.cov.40.AD.women.cl.40.50 <- AD.MAR.b.cov.40[,17]
vector.sd.MAR.b.cov.40.AD.men.cl.40.50 <- AD.MAR.b.cov.40[,18]



# 45



d.MAR.b.cov.45 <- d.MAR %>%
  select(contains("cov.MAR.b.45."))
AD.MAR.b.cov.45 <- d.MAR.b.cov.45 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.b.cov.45.AD.women.cl.15.25 <- AD.MAR.b.cov.45[,1]
vector.mean.MAR.b.cov.45.AD.men.cl.15.25 <- AD.MAR.b.cov.45[,2]

vector.mean.MAR.b.cov.45.AD.women.cl.25.40 <- AD.MAR.b.cov.45[,3]
vector.mean.MAR.b.cov.45.AD.men.cl.25.40 <- AD.MAR.b.cov.45[,4]

vector.mean.MAR.b.cov.45.AD.women.cl.40.50 <- AD.MAR.b.cov.45[,5]
vector.mean.MAR.b.cov.45.AD.men.cl.40.50 <- AD.MAR.b.cov.45[,6]

# Median

vector.med.MAR.b.cov.45.AD.women.cl.15.25 <- AD.MAR.b.cov.45[,7]
vector.med.MAR.b.cov.45.AD.men.cl.15.25 <- AD.MAR.b.cov.45[,8]

vector.med.MAR.b.cov.45.AD.women.cl.25.40 <- AD.MAR.b.cov.45[,9]
vector.med.MAR.b.cov.45.AD.men.cl.25.40 <- AD.MAR.b.cov.45[,10]

vector.med.MAR.b.cov.45.AD.women.cl.40.50 <- AD.MAR.b.cov.45[,11]
vector.med.MAR.b.cov.45.AD.men.cl.40.50 <- AD.MAR.b.cov.45[,12]

# Standard deviation

vector.sd.MAR.b.cov.45.AD.women.cl.15.25 <- AD.MAR.b.cov.45[,13]
vector.sd.MAR.b.cov.45.AD.men.cl.15.25 <- AD.MAR.b.cov.45[,14]

vector.sd.MAR.b.cov.45.AD.women.cl.25.40 <- AD.MAR.b.cov.45[,15]
vector.sd.MAR.b.cov.45.AD.men.cl.25.40 <- AD.MAR.b.cov.45[,16]

vector.sd.MAR.b.cov.45.AD.women.cl.40.50 <- AD.MAR.b.cov.45[,17]
vector.sd.MAR.b.cov.45.AD.men.cl.40.50 <- AD.MAR.b.cov.45[,18]



# 50



d.MAR.b.cov.50 <- d.MAR %>%
  select(contains("cov.MAR.b.50."))
AD.MAR.b.cov.50 <- d.MAR.b.cov.50 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.b.cov.50.AD.women.cl.15.25 <- AD.MAR.b.cov.50[,1]
vector.mean.MAR.b.cov.50.AD.men.cl.15.25 <- AD.MAR.b.cov.50[,2]

vector.mean.MAR.b.cov.50.AD.women.cl.25.40 <- AD.MAR.b.cov.50[,3]
vector.mean.MAR.b.cov.50.AD.men.cl.25.40 <- AD.MAR.b.cov.50[,4]

vector.mean.MAR.b.cov.50.AD.women.cl.40.50 <- AD.MAR.b.cov.50[,5]
vector.mean.MAR.b.cov.50.AD.men.cl.40.50 <- AD.MAR.b.cov.50[,6]

# Median

vector.med.MAR.b.cov.50.AD.women.cl.15.25 <- AD.MAR.b.cov.50[,7]
vector.med.MAR.b.cov.50.AD.men.cl.15.25 <- AD.MAR.b.cov.50[,8]

vector.med.MAR.b.cov.50.AD.women.cl.25.40 <- AD.MAR.b.cov.50[,9]
vector.med.MAR.b.cov.50.AD.men.cl.25.40 <- AD.MAR.b.cov.50[,10]

vector.med.MAR.b.cov.50.AD.women.cl.40.50 <- AD.MAR.b.cov.50[,11]
vector.med.MAR.b.cov.50.AD.men.cl.40.50 <- AD.MAR.b.cov.50[,12]

# Standard deviation

vector.sd.MAR.b.cov.50.AD.women.cl.15.25 <- AD.MAR.b.cov.50[,13]
vector.sd.MAR.b.cov.50.AD.men.cl.15.25 <- AD.MAR.b.cov.50[,14]

vector.sd.MAR.b.cov.50.AD.women.cl.25.40 <- AD.MAR.b.cov.50[,15]
vector.sd.MAR.b.cov.50.AD.men.cl.25.40 <- AD.MAR.b.cov.50[,16]

vector.sd.MAR.b.cov.50.AD.women.cl.40.50 <- AD.MAR.b.cov.50[,17]
vector.sd.MAR.b.cov.50.AD.men.cl.40.50 <- AD.MAR.b.cov.50[,18]



# 55



d.MAR.b.cov.55 <- d.MAR %>%
  select(contains("cov.MAR.b.55."))
AD.MAR.b.cov.55 <- d.MAR.b.cov.55 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.b.cov.55.AD.women.cl.15.25 <- AD.MAR.b.cov.55[,1]
vector.mean.MAR.b.cov.55.AD.men.cl.15.25 <- AD.MAR.b.cov.55[,2]

vector.mean.MAR.b.cov.55.AD.women.cl.25.40 <- AD.MAR.b.cov.55[,3]
vector.mean.MAR.b.cov.55.AD.men.cl.25.40 <- AD.MAR.b.cov.55[,4]

vector.mean.MAR.b.cov.55.AD.women.cl.40.50 <- AD.MAR.b.cov.55[,5]
vector.mean.MAR.b.cov.55.AD.men.cl.40.50 <- AD.MAR.b.cov.55[,6]

# Median

vector.med.MAR.b.cov.55.AD.women.cl.15.25 <- AD.MAR.b.cov.55[,7]
vector.med.MAR.b.cov.55.AD.men.cl.15.25 <- AD.MAR.b.cov.55[,8]

vector.med.MAR.b.cov.55.AD.women.cl.25.40 <- AD.MAR.b.cov.55[,9]
vector.med.MAR.b.cov.55.AD.men.cl.25.40 <- AD.MAR.b.cov.55[,10]

vector.med.MAR.b.cov.55.AD.women.cl.40.50 <- AD.MAR.b.cov.55[,11]
vector.med.MAR.b.cov.55.AD.men.cl.40.50 <- AD.MAR.b.cov.55[,12]

# Standard deviation

vector.sd.MAR.b.cov.55.AD.women.cl.15.25 <- AD.MAR.b.cov.55[,13]
vector.sd.MAR.b.cov.55.AD.men.cl.15.25 <- AD.MAR.b.cov.55[,14]

vector.sd.MAR.b.cov.55.AD.women.cl.25.40 <- AD.MAR.b.cov.55[,15]
vector.sd.MAR.b.cov.55.AD.men.cl.25.40 <- AD.MAR.b.cov.55[,16]

vector.sd.MAR.b.cov.55.AD.women.cl.40.50 <- AD.MAR.b.cov.55[,17]
vector.sd.MAR.b.cov.55.AD.men.cl.40.50 <- AD.MAR.b.cov.55[,18]



# 60



d.MAR.b.cov.60 <- d.MAR %>%
  select(contains("cov.MAR.b.60."))
AD.MAR.b.cov.60 <- d.MAR.b.cov.60 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.b.cov.60.AD.women.cl.15.25 <- AD.MAR.b.cov.60[,1]
vector.mean.MAR.b.cov.60.AD.men.cl.15.25 <- AD.MAR.b.cov.60[,2]

vector.mean.MAR.b.cov.60.AD.women.cl.25.40 <- AD.MAR.b.cov.60[,3]
vector.mean.MAR.b.cov.60.AD.men.cl.25.40 <- AD.MAR.b.cov.60[,4]

vector.mean.MAR.b.cov.60.AD.women.cl.40.50 <- AD.MAR.b.cov.60[,5]
vector.mean.MAR.b.cov.60.AD.men.cl.40.50 <- AD.MAR.b.cov.60[,6]

# Median

vector.med.MAR.b.cov.60.AD.women.cl.15.25 <- AD.MAR.b.cov.60[,7]
vector.med.MAR.b.cov.60.AD.men.cl.15.25 <- AD.MAR.b.cov.60[,8]

vector.med.MAR.b.cov.60.AD.women.cl.25.40 <- AD.MAR.b.cov.60[,9]
vector.med.MAR.b.cov.60.AD.men.cl.25.40 <- AD.MAR.b.cov.60[,10]

vector.med.MAR.b.cov.60.AD.women.cl.40.50 <- AD.MAR.b.cov.60[,11]
vector.med.MAR.b.cov.60.AD.men.cl.40.50 <- AD.MAR.b.cov.60[,12]

# Standard deviation

vector.sd.MAR.b.cov.60.AD.women.cl.15.25 <- AD.MAR.b.cov.60[,13]
vector.sd.MAR.b.cov.60.AD.men.cl.15.25 <- AD.MAR.b.cov.60[,14]

vector.sd.MAR.b.cov.60.AD.women.cl.25.40 <- AD.MAR.b.cov.60[,15]
vector.sd.MAR.b.cov.60.AD.men.cl.25.40 <- AD.MAR.b.cov.60[,16]

vector.sd.MAR.b.cov.60.AD.women.cl.40.50 <- AD.MAR.b.cov.60[,17]
vector.sd.MAR.b.cov.60.AD.men.cl.40.50 <- AD.MAR.b.cov.60[,18]



# 65



d.MAR.b.cov.65 <- d.MAR %>%
  select(contains("cov.MAR.b.65."))
AD.MAR.b.cov.65 <- d.MAR.b.cov.65 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.b.cov.65.AD.women.cl.15.25 <- AD.MAR.b.cov.65[,1]
vector.mean.MAR.b.cov.65.AD.men.cl.15.25 <- AD.MAR.b.cov.65[,2]

vector.mean.MAR.b.cov.65.AD.women.cl.25.40 <- AD.MAR.b.cov.65[,3]
vector.mean.MAR.b.cov.65.AD.men.cl.25.40 <- AD.MAR.b.cov.65[,4]

vector.mean.MAR.b.cov.65.AD.women.cl.40.50 <- AD.MAR.b.cov.65[,5]
vector.mean.MAR.b.cov.65.AD.men.cl.40.50 <- AD.MAR.b.cov.65[,6]

# Median

vector.med.MAR.b.cov.65.AD.women.cl.15.25 <- AD.MAR.b.cov.65[,7]
vector.med.MAR.b.cov.65.AD.men.cl.15.25 <- AD.MAR.b.cov.65[,8]

vector.med.MAR.b.cov.65.AD.women.cl.25.40 <- AD.MAR.b.cov.65[,9]
vector.med.MAR.b.cov.65.AD.men.cl.25.40 <- AD.MAR.b.cov.65[,10]

vector.med.MAR.b.cov.65.AD.women.cl.40.50 <- AD.MAR.b.cov.65[,11]
vector.med.MAR.b.cov.65.AD.men.cl.40.50 <- AD.MAR.b.cov.65[,12]

# Standard deviation

vector.sd.MAR.b.cov.65.AD.women.cl.15.25 <- AD.MAR.b.cov.65[,13]
vector.sd.MAR.b.cov.65.AD.men.cl.15.25 <- AD.MAR.b.cov.65[,14]

vector.sd.MAR.b.cov.65.AD.women.cl.25.40 <- AD.MAR.b.cov.65[,15]
vector.sd.MAR.b.cov.65.AD.men.cl.25.40 <- AD.MAR.b.cov.65[,16]

vector.sd.MAR.b.cov.65.AD.women.cl.40.50 <- AD.MAR.b.cov.65[,17]
vector.sd.MAR.b.cov.65.AD.men.cl.40.50 <- AD.MAR.b.cov.65[,18]



# 70



d.MAR.b.cov.70 <- d.MAR %>%
  select(contains("cov.MAR.b.70."))
AD.MAR.b.cov.70 <- d.MAR.b.cov.70 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.b.cov.70.AD.women.cl.15.25 <- AD.MAR.b.cov.70[,1]
vector.mean.MAR.b.cov.70.AD.men.cl.15.25 <- AD.MAR.b.cov.70[,2]

vector.mean.MAR.b.cov.70.AD.women.cl.25.40 <- AD.MAR.b.cov.70[,3]
vector.mean.MAR.b.cov.70.AD.men.cl.25.40 <- AD.MAR.b.cov.70[,4]

vector.mean.MAR.b.cov.70.AD.women.cl.40.50 <- AD.MAR.b.cov.70[,5]
vector.mean.MAR.b.cov.70.AD.men.cl.40.50 <- AD.MAR.b.cov.70[,6]

# Median

vector.med.MAR.b.cov.70.AD.women.cl.15.25 <- AD.MAR.b.cov.70[,7]
vector.med.MAR.b.cov.70.AD.men.cl.15.25 <- AD.MAR.b.cov.70[,8]

vector.med.MAR.b.cov.70.AD.women.cl.25.40 <- AD.MAR.b.cov.70[,9]
vector.med.MAR.b.cov.70.AD.men.cl.25.40 <- AD.MAR.b.cov.70[,10]

vector.med.MAR.b.cov.70.AD.women.cl.40.50 <- AD.MAR.b.cov.70[,11]
vector.med.MAR.b.cov.70.AD.men.cl.40.50 <- AD.MAR.b.cov.70[,12]

# Standard deviation

vector.sd.MAR.b.cov.70.AD.women.cl.15.25 <- AD.MAR.b.cov.70[,13]
vector.sd.MAR.b.cov.70.AD.men.cl.15.25 <- AD.MAR.b.cov.70[,14]

vector.sd.MAR.b.cov.70.AD.women.cl.25.40 <- AD.MAR.b.cov.70[,15]
vector.sd.MAR.b.cov.70.AD.men.cl.25.40 <- AD.MAR.b.cov.70[,16]

vector.sd.MAR.b.cov.70.AD.women.cl.40.50 <- AD.MAR.b.cov.70[,17]
vector.sd.MAR.b.cov.70.AD.men.cl.40.50 <- AD.MAR.b.cov.70[,18]



# 75




d.MAR.b.cov.75 <- d.MAR %>%
  select(contains("cov.MAR.b.75."))
AD.MAR.b.cov.75 <- d.MAR.b.cov.75 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.b.cov.75.AD.women.cl.15.25 <- AD.MAR.b.cov.75[,1]
vector.mean.MAR.b.cov.75.AD.men.cl.15.25 <- AD.MAR.b.cov.75[,2]

vector.mean.MAR.b.cov.75.AD.women.cl.25.40 <- AD.MAR.b.cov.75[,3]
vector.mean.MAR.b.cov.75.AD.men.cl.25.40 <- AD.MAR.b.cov.75[,4]

vector.mean.MAR.b.cov.75.AD.women.cl.40.50 <- AD.MAR.b.cov.75[,5]
vector.mean.MAR.b.cov.75.AD.men.cl.40.50 <- AD.MAR.b.cov.75[,6]

# Median

vector.med.MAR.b.cov.75.AD.women.cl.15.25 <- AD.MAR.b.cov.75[,7]
vector.med.MAR.b.cov.75.AD.men.cl.15.25 <- AD.MAR.b.cov.75[,8]

vector.med.MAR.b.cov.75.AD.women.cl.25.40 <- AD.MAR.b.cov.75[,9]
vector.med.MAR.b.cov.75.AD.men.cl.25.40 <- AD.MAR.b.cov.75[,10]

vector.med.MAR.b.cov.75.AD.women.cl.40.50 <- AD.MAR.b.cov.75[,11]
vector.med.MAR.b.cov.75.AD.men.cl.40.50 <- AD.MAR.b.cov.75[,12]

# Standard deviation

vector.sd.MAR.b.cov.75.AD.women.cl.15.25 <- AD.MAR.b.cov.75[,13]
vector.sd.MAR.b.cov.75.AD.men.cl.15.25 <- AD.MAR.b.cov.75[,14]

vector.sd.MAR.b.cov.75.AD.women.cl.25.40 <- AD.MAR.b.cov.75[,15]
vector.sd.MAR.b.cov.75.AD.men.cl.25.40 <- AD.MAR.b.cov.75[,16]

vector.sd.MAR.b.cov.75.AD.women.cl.40.50 <- AD.MAR.b.cov.75[,17]
vector.sd.MAR.b.cov.75.AD.men.cl.40.50 <- AD.MAR.b.cov.75[,18]



# 80



d.MAR.b.cov.80 <- d.MAR %>%
  select(contains("cov.MAR.b.80."))
AD.MAR.b.cov.80 <- d.MAR.b.cov.80 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.b.cov.80.AD.women.cl.15.25 <- AD.MAR.b.cov.80[,1]
vector.mean.MAR.b.cov.80.AD.men.cl.15.25 <- AD.MAR.b.cov.80[,2]

vector.mean.MAR.b.cov.80.AD.women.cl.25.40 <- AD.MAR.b.cov.80[,3]
vector.mean.MAR.b.cov.80.AD.men.cl.25.40 <- AD.MAR.b.cov.80[,4]

vector.mean.MAR.b.cov.80.AD.women.cl.40.50 <- AD.MAR.b.cov.80[,5]
vector.mean.MAR.b.cov.80.AD.men.cl.40.50 <- AD.MAR.b.cov.80[,6]

# Median

vector.med.MAR.b.cov.80.AD.women.cl.15.25 <- AD.MAR.b.cov.80[,7]
vector.med.MAR.b.cov.80.AD.men.cl.15.25 <- AD.MAR.b.cov.80[,8]

vector.med.MAR.b.cov.80.AD.women.cl.25.40 <- AD.MAR.b.cov.80[,9]
vector.med.MAR.b.cov.80.AD.men.cl.25.40 <- AD.MAR.b.cov.80[,10]

vector.med.MAR.b.cov.80.AD.women.cl.40.50 <- AD.MAR.b.cov.80[,11]
vector.med.MAR.b.cov.80.AD.men.cl.40.50 <- AD.MAR.b.cov.80[,12]

# Standard deviation

vector.sd.MAR.b.cov.80.AD.women.cl.15.25 <- AD.MAR.b.cov.80[,13]
vector.sd.MAR.b.cov.80.AD.men.cl.15.25 <- AD.MAR.b.cov.80[,14]

vector.sd.MAR.b.cov.80.AD.women.cl.25.40 <- AD.MAR.b.cov.80[,15]
vector.sd.MAR.b.cov.80.AD.men.cl.25.40 <- AD.MAR.b.cov.80[,16]

vector.sd.MAR.b.cov.80.AD.women.cl.40.50 <- AD.MAR.b.cov.80[,17]
vector.sd.MAR.b.cov.80.AD.men.cl.40.50 <- AD.MAR.b.cov.80[,18]



# 85



d.MAR.b.cov.85 <- d.MAR %>%
  select(contains("cov.MAR.b.85."))
AD.MAR.b.cov.85 <- d.MAR.b.cov.85 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.b.cov.85.AD.women.cl.15.25 <- AD.MAR.b.cov.85[,1]
vector.mean.MAR.b.cov.85.AD.men.cl.15.25 <- AD.MAR.b.cov.85[,2]

vector.mean.MAR.b.cov.85.AD.women.cl.25.40 <- AD.MAR.b.cov.85[,3]
vector.mean.MAR.b.cov.85.AD.men.cl.25.40 <- AD.MAR.b.cov.85[,4]

vector.mean.MAR.b.cov.85.AD.women.cl.40.50 <- AD.MAR.b.cov.85[,5]
vector.mean.MAR.b.cov.85.AD.men.cl.40.50 <- AD.MAR.b.cov.85[,6]

# Median

vector.med.MAR.b.cov.85.AD.women.cl.15.25 <- AD.MAR.b.cov.85[,7]
vector.med.MAR.b.cov.85.AD.men.cl.15.25 <- AD.MAR.b.cov.85[,8]

vector.med.MAR.b.cov.85.AD.women.cl.25.40 <- AD.MAR.b.cov.85[,9]
vector.med.MAR.b.cov.85.AD.men.cl.25.40 <- AD.MAR.b.cov.85[,10]

vector.med.MAR.b.cov.85.AD.women.cl.40.50 <- AD.MAR.b.cov.85[,11]
vector.med.MAR.b.cov.85.AD.men.cl.40.50 <- AD.MAR.b.cov.85[,12]

# Standard deviation

vector.sd.MAR.b.cov.85.AD.women.cl.15.25 <- AD.MAR.b.cov.85[,13]
vector.sd.MAR.b.cov.85.AD.men.cl.15.25 <- AD.MAR.b.cov.85[,14]

vector.sd.MAR.b.cov.85.AD.women.cl.25.40 <- AD.MAR.b.cov.85[,15]
vector.sd.MAR.b.cov.85.AD.men.cl.25.40 <- AD.MAR.b.cov.85[,16]

vector.sd.MAR.b.cov.85.AD.women.cl.40.50 <- AD.MAR.b.cov.85[,17]
vector.sd.MAR.b.cov.85.AD.men.cl.40.50 <- AD.MAR.b.cov.85[,18]



# 90



d.MAR.b.cov.90 <- d.MAR %>%
  select(contains("cov.MAR.b.90."))
AD.MAR.b.cov.90 <- d.MAR.b.cov.90 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.b.cov.90.AD.women.cl.15.25 <- AD.MAR.b.cov.90[,1]
vector.mean.MAR.b.cov.90.AD.men.cl.15.25 <- AD.MAR.b.cov.90[,2]

vector.mean.MAR.b.cov.90.AD.women.cl.25.40 <- AD.MAR.b.cov.90[,3]
vector.mean.MAR.b.cov.90.AD.men.cl.25.40 <- AD.MAR.b.cov.90[,4]

vector.mean.MAR.b.cov.90.AD.women.cl.40.50 <- AD.MAR.b.cov.90[,5]
vector.mean.MAR.b.cov.90.AD.men.cl.40.50 <- AD.MAR.b.cov.90[,6]

# Median

vector.med.MAR.b.cov.90.AD.women.cl.15.25 <- AD.MAR.b.cov.90[,7]
vector.med.MAR.b.cov.90.AD.men.cl.15.25 <- AD.MAR.b.cov.90[,8]

vector.med.MAR.b.cov.90.AD.women.cl.25.40 <- AD.MAR.b.cov.90[,9]
vector.med.MAR.b.cov.90.AD.men.cl.25.40 <- AD.MAR.b.cov.90[,10]

vector.med.MAR.b.cov.90.AD.women.cl.40.50 <- AD.MAR.b.cov.90[,11]
vector.med.MAR.b.cov.90.AD.men.cl.40.50 <- AD.MAR.b.cov.90[,12]

# Standard deviation

vector.sd.MAR.b.cov.90.AD.women.cl.15.25 <- AD.MAR.b.cov.90[,13]
vector.sd.MAR.b.cov.90.AD.men.cl.15.25 <- AD.MAR.b.cov.90[,14]

vector.sd.MAR.b.cov.90.AD.women.cl.25.40 <- AD.MAR.b.cov.90[,15]
vector.sd.MAR.b.cov.90.AD.men.cl.25.40 <- AD.MAR.b.cov.90[,16]

vector.sd.MAR.b.cov.90.AD.women.cl.40.50 <- AD.MAR.b.cov.90[,17]
vector.sd.MAR.b.cov.90.AD.men.cl.40.50 <- AD.MAR.b.cov.90[,18]



# 95



d.MAR.b.cov.95 <- d.MAR %>%
  select(contains("cov.MAR.b.95."))
AD.MAR.b.cov.95 <- d.MAR.b.cov.95 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.b.cov.95.AD.women.cl.15.25 <- AD.MAR.b.cov.95[,1]
vector.mean.MAR.b.cov.95.AD.men.cl.15.25 <- AD.MAR.b.cov.95[,2]

vector.mean.MAR.b.cov.95.AD.women.cl.25.40 <- AD.MAR.b.cov.95[,3]
vector.mean.MAR.b.cov.95.AD.men.cl.25.40 <- AD.MAR.b.cov.95[,4]

vector.mean.MAR.b.cov.95.AD.women.cl.40.50 <- AD.MAR.b.cov.95[,5]
vector.mean.MAR.b.cov.95.AD.men.cl.40.50 <- AD.MAR.b.cov.95[,6]

# Median

vector.med.MAR.b.cov.95.AD.women.cl.15.25 <- AD.MAR.b.cov.95[,7]
vector.med.MAR.b.cov.95.AD.men.cl.15.25 <- AD.MAR.b.cov.95[,8]

vector.med.MAR.b.cov.95.AD.women.cl.25.40 <- AD.MAR.b.cov.95[,9]
vector.med.MAR.b.cov.95.AD.men.cl.25.40 <- AD.MAR.b.cov.95[,10]

vector.med.MAR.b.cov.95.AD.women.cl.40.50 <- AD.MAR.b.cov.95[,11]
vector.med.MAR.b.cov.95.AD.men.cl.40.50 <- AD.MAR.b.cov.95[,12]

# Standard deviation

vector.sd.MAR.b.cov.95.AD.women.cl.15.25 <- AD.MAR.b.cov.95[,13]
vector.sd.MAR.b.cov.95.AD.men.cl.15.25 <- AD.MAR.b.cov.95[,14]

vector.sd.MAR.b.cov.95.AD.women.cl.25.40 <- AD.MAR.b.cov.95[,15]
vector.sd.MAR.b.cov.95.AD.men.cl.25.40 <- AD.MAR.b.cov.95[,16]

vector.sd.MAR.b.cov.95.AD.women.cl.40.50 <- AD.MAR.b.cov.95[,17]
vector.sd.MAR.b.cov.95.AD.men.cl.40.50 <- AD.MAR.b.cov.95[,18]




# MAR - c



# 35

d.MAR.c.cov.35 <- d.MAR %>%
  select(contains("cov.MAR.c.35.")) # MAR - a
AD.MAR.c.cov.35 <- d.MAR.c.cov.35 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.c.cov.35.AD.women.cl.15.25 <- AD.MAR.c.cov.35[,1]
vector.mean.MAR.c.cov.35.AD.men.cl.15.25 <- AD.MAR.c.cov.35[,2]

vector.mean.MAR.c.cov.35.AD.women.cl.25.40 <- AD.MAR.c.cov.35[,3]
vector.mean.MAR.c.cov.35.AD.men.cl.25.40 <- AD.MAR.c.cov.35[,4]

vector.mean.MAR.c.cov.35.AD.women.cl.40.50 <- AD.MAR.c.cov.35[,5]
vector.mean.MAR.c.cov.35.AD.men.cl.40.50 <- AD.MAR.c.cov.35[,6]

# Median

vector.med.MAR.c.cov.35.AD.women.cl.15.25 <- AD.MAR.c.cov.35[,7]
vector.med.MAR.c.cov.35.AD.men.cl.15.25 <- AD.MAR.c.cov.35[,8]

vector.med.MAR.c.cov.35.AD.women.cl.25.40 <- AD.MAR.c.cov.35[,9]
vector.med.MAR.c.cov.35.AD.men.cl.25.40 <- AD.MAR.c.cov.35[,10]

vector.med.MAR.c.cov.35.AD.women.cl.40.50 <- AD.MAR.c.cov.35[,11]
vector.med.MAR.c.cov.35.AD.men.cl.40.50 <- AD.MAR.c.cov.35[,12]

# Standard deviation

vector.sd.MAR.c.cov.35.AD.women.cl.15.25 <- AD.MAR.c.cov.35[,13]
vector.sd.MAR.c.cov.35.AD.men.cl.15.25 <- AD.MAR.c.cov.35[,14]

vector.sd.MAR.c.cov.35.AD.women.cl.25.40 <- AD.MAR.c.cov.35[,15]
vector.sd.MAR.c.cov.35.AD.men.cl.25.40 <- AD.MAR.c.cov.35[,16]

vector.sd.MAR.c.cov.35.AD.women.cl.40.50 <- AD.MAR.c.cov.35[,17]
vector.sd.MAR.c.cov.35.AD.men.cl.40.50 <- AD.MAR.c.cov.35[,18]


# 40



d.MAR.c.cov.40 <- d.MAR %>%
  select(contains("cov.MAR.c.40.")) # MAR - a
AD.MAR.c.cov.40 <- d.MAR.c.cov.40 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.c.cov.40.AD.women.cl.15.25 <- AD.MAR.c.cov.40[,1]
vector.mean.MAR.c.cov.40.AD.men.cl.15.25 <- AD.MAR.c.cov.40[,2]

vector.mean.MAR.c.cov.40.AD.women.cl.25.40 <- AD.MAR.c.cov.40[,3]
vector.mean.MAR.c.cov.40.AD.men.cl.25.40 <- AD.MAR.c.cov.40[,4]

vector.mean.MAR.c.cov.40.AD.women.cl.40.50 <- AD.MAR.c.cov.40[,5]
vector.mean.MAR.c.cov.40.AD.men.cl.40.50 <- AD.MAR.c.cov.40[,6]

# Median

vector.med.MAR.c.cov.40.AD.women.cl.15.25 <- AD.MAR.c.cov.40[,7]
vector.med.MAR.c.cov.40.AD.men.cl.15.25 <- AD.MAR.c.cov.40[,8]

vector.med.MAR.c.cov.40.AD.women.cl.25.40 <- AD.MAR.c.cov.40[,9]
vector.med.MAR.c.cov.40.AD.men.cl.25.40 <- AD.MAR.c.cov.40[,10]

vector.med.MAR.c.cov.40.AD.women.cl.40.50 <- AD.MAR.c.cov.40[,11]
vector.med.MAR.c.cov.40.AD.men.cl.40.50 <- AD.MAR.c.cov.40[,12]

# Standard deviation

vector.sd.MAR.c.cov.40.AD.women.cl.15.25 <- AD.MAR.c.cov.40[,13]
vector.sd.MAR.c.cov.40.AD.men.cl.15.25 <- AD.MAR.c.cov.40[,14]

vector.sd.MAR.c.cov.40.AD.women.cl.25.40 <- AD.MAR.c.cov.40[,15]
vector.sd.MAR.c.cov.40.AD.men.cl.25.40 <- AD.MAR.c.cov.40[,16]

vector.sd.MAR.c.cov.40.AD.women.cl.40.50 <- AD.MAR.c.cov.40[,17]
vector.sd.MAR.c.cov.40.AD.men.cl.40.50 <- AD.MAR.c.cov.40[,18]



# 45



d.MAR.c.cov.45 <- d.MAR %>%
  select(contains("cov.MAR.c.45."))
AD.MAR.c.cov.45 <- d.MAR.c.cov.45 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.c.cov.45.AD.women.cl.15.25 <- AD.MAR.c.cov.45[,1]
vector.mean.MAR.c.cov.45.AD.men.cl.15.25 <- AD.MAR.c.cov.45[,2]

vector.mean.MAR.c.cov.45.AD.women.cl.25.40 <- AD.MAR.c.cov.45[,3]
vector.mean.MAR.c.cov.45.AD.men.cl.25.40 <- AD.MAR.c.cov.45[,4]

vector.mean.MAR.c.cov.45.AD.women.cl.40.50 <- AD.MAR.c.cov.45[,5]
vector.mean.MAR.c.cov.45.AD.men.cl.40.50 <- AD.MAR.c.cov.45[,6]

# Median

vector.med.MAR.c.cov.45.AD.women.cl.15.25 <- AD.MAR.c.cov.45[,7]
vector.med.MAR.c.cov.45.AD.men.cl.15.25 <- AD.MAR.c.cov.45[,8]

vector.med.MAR.c.cov.45.AD.women.cl.25.40 <- AD.MAR.c.cov.45[,9]
vector.med.MAR.c.cov.45.AD.men.cl.25.40 <- AD.MAR.c.cov.45[,10]

vector.med.MAR.c.cov.45.AD.women.cl.40.50 <- AD.MAR.c.cov.45[,11]
vector.med.MAR.c.cov.45.AD.men.cl.40.50 <- AD.MAR.c.cov.45[,12]

# Standard deviation

vector.sd.MAR.c.cov.45.AD.women.cl.15.25 <- AD.MAR.c.cov.45[,13]
vector.sd.MAR.c.cov.45.AD.men.cl.15.25 <- AD.MAR.c.cov.45[,14]

vector.sd.MAR.c.cov.45.AD.women.cl.25.40 <- AD.MAR.c.cov.45[,15]
vector.sd.MAR.c.cov.45.AD.men.cl.25.40 <- AD.MAR.c.cov.45[,16]

vector.sd.MAR.c.cov.45.AD.women.cl.40.50 <- AD.MAR.c.cov.45[,17]
vector.sd.MAR.c.cov.45.AD.men.cl.40.50 <- AD.MAR.c.cov.45[,18]



# 50



d.MAR.c.cov.50 <- d.MAR %>%
  select(contains("cov.MAR.c.50."))
AD.MAR.c.cov.50 <- d.MAR.c.cov.50 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.c.cov.50.AD.women.cl.15.25 <- AD.MAR.c.cov.50[,1]
vector.mean.MAR.c.cov.50.AD.men.cl.15.25 <- AD.MAR.c.cov.50[,2]

vector.mean.MAR.c.cov.50.AD.women.cl.25.40 <- AD.MAR.c.cov.50[,3]
vector.mean.MAR.c.cov.50.AD.men.cl.25.40 <- AD.MAR.c.cov.50[,4]

vector.mean.MAR.c.cov.50.AD.women.cl.40.50 <- AD.MAR.c.cov.50[,5]
vector.mean.MAR.c.cov.50.AD.men.cl.40.50 <- AD.MAR.c.cov.50[,6]

# Median

vector.med.MAR.c.cov.50.AD.women.cl.15.25 <- AD.MAR.c.cov.50[,7]
vector.med.MAR.c.cov.50.AD.men.cl.15.25 <- AD.MAR.c.cov.50[,8]

vector.med.MAR.c.cov.50.AD.women.cl.25.40 <- AD.MAR.c.cov.50[,9]
vector.med.MAR.c.cov.50.AD.men.cl.25.40 <- AD.MAR.c.cov.50[,10]

vector.med.MAR.c.cov.50.AD.women.cl.40.50 <- AD.MAR.c.cov.50[,11]
vector.med.MAR.c.cov.50.AD.men.cl.40.50 <- AD.MAR.c.cov.50[,12]

# Standard deviation

vector.sd.MAR.c.cov.50.AD.women.cl.15.25 <- AD.MAR.c.cov.50[,13]
vector.sd.MAR.c.cov.50.AD.men.cl.15.25 <- AD.MAR.c.cov.50[,14]

vector.sd.MAR.c.cov.50.AD.women.cl.25.40 <- AD.MAR.c.cov.50[,15]
vector.sd.MAR.c.cov.50.AD.men.cl.25.40 <- AD.MAR.c.cov.50[,16]

vector.sd.MAR.c.cov.50.AD.women.cl.40.50 <- AD.MAR.c.cov.50[,17]
vector.sd.MAR.c.cov.50.AD.men.cl.40.50 <- AD.MAR.c.cov.50[,18]



# 55



d.MAR.c.cov.55 <- d.MAR %>%
  select(contains("cov.MAR.c.55."))
AD.MAR.c.cov.55 <- d.MAR.c.cov.55 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.c.cov.55.AD.women.cl.15.25 <- AD.MAR.c.cov.55[,1]
vector.mean.MAR.c.cov.55.AD.men.cl.15.25 <- AD.MAR.c.cov.55[,2]

vector.mean.MAR.c.cov.55.AD.women.cl.25.40 <- AD.MAR.c.cov.55[,3]
vector.mean.MAR.c.cov.55.AD.men.cl.25.40 <- AD.MAR.c.cov.55[,4]

vector.mean.MAR.c.cov.55.AD.women.cl.40.50 <- AD.MAR.c.cov.55[,5]
vector.mean.MAR.c.cov.55.AD.men.cl.40.50 <- AD.MAR.c.cov.55[,6]

# Median

vector.med.MAR.c.cov.55.AD.women.cl.15.25 <- AD.MAR.c.cov.55[,7]
vector.med.MAR.c.cov.55.AD.men.cl.15.25 <- AD.MAR.c.cov.55[,8]

vector.med.MAR.c.cov.55.AD.women.cl.25.40 <- AD.MAR.c.cov.55[,9]
vector.med.MAR.c.cov.55.AD.men.cl.25.40 <- AD.MAR.c.cov.55[,10]

vector.med.MAR.c.cov.55.AD.women.cl.40.50 <- AD.MAR.c.cov.55[,11]
vector.med.MAR.c.cov.55.AD.men.cl.40.50 <- AD.MAR.c.cov.55[,12]

# Standard deviation

vector.sd.MAR.c.cov.55.AD.women.cl.15.25 <- AD.MAR.c.cov.55[,13]
vector.sd.MAR.c.cov.55.AD.men.cl.15.25 <- AD.MAR.c.cov.55[,14]

vector.sd.MAR.c.cov.55.AD.women.cl.25.40 <- AD.MAR.c.cov.55[,15]
vector.sd.MAR.c.cov.55.AD.men.cl.25.40 <- AD.MAR.c.cov.55[,16]

vector.sd.MAR.c.cov.55.AD.women.cl.40.50 <- AD.MAR.c.cov.55[,17]
vector.sd.MAR.c.cov.55.AD.men.cl.40.50 <- AD.MAR.c.cov.55[,18]



# 60



d.MAR.c.cov.60 <- d.MAR %>%
  select(contains("cov.MAR.c.60."))
AD.MAR.c.cov.60 <- d.MAR.c.cov.60 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.c.cov.60.AD.women.cl.15.25 <- AD.MAR.c.cov.60[,1]
vector.mean.MAR.c.cov.60.AD.men.cl.15.25 <- AD.MAR.c.cov.60[,2]

vector.mean.MAR.c.cov.60.AD.women.cl.25.40 <- AD.MAR.c.cov.60[,3]
vector.mean.MAR.c.cov.60.AD.men.cl.25.40 <- AD.MAR.c.cov.60[,4]

vector.mean.MAR.c.cov.60.AD.women.cl.40.50 <- AD.MAR.c.cov.60[,5]
vector.mean.MAR.c.cov.60.AD.men.cl.40.50 <- AD.MAR.c.cov.60[,6]

# Median

vector.med.MAR.c.cov.60.AD.women.cl.15.25 <- AD.MAR.c.cov.60[,7]
vector.med.MAR.c.cov.60.AD.men.cl.15.25 <- AD.MAR.c.cov.60[,8]

vector.med.MAR.c.cov.60.AD.women.cl.25.40 <- AD.MAR.c.cov.60[,9]
vector.med.MAR.c.cov.60.AD.men.cl.25.40 <- AD.MAR.c.cov.60[,10]

vector.med.MAR.c.cov.60.AD.women.cl.40.50 <- AD.MAR.c.cov.60[,11]
vector.med.MAR.c.cov.60.AD.men.cl.40.50 <- AD.MAR.c.cov.60[,12]

# Standard deviation

vector.sd.MAR.c.cov.60.AD.women.cl.15.25 <- AD.MAR.c.cov.60[,13]
vector.sd.MAR.c.cov.60.AD.men.cl.15.25 <- AD.MAR.c.cov.60[,14]

vector.sd.MAR.c.cov.60.AD.women.cl.25.40 <- AD.MAR.c.cov.60[,15]
vector.sd.MAR.c.cov.60.AD.men.cl.25.40 <- AD.MAR.c.cov.60[,16]

vector.sd.MAR.c.cov.60.AD.women.cl.40.50 <- AD.MAR.c.cov.60[,17]
vector.sd.MAR.c.cov.60.AD.men.cl.40.50 <- AD.MAR.c.cov.60[,18]



# 65



d.MAR.c.cov.65 <- d.MAR %>%
  select(contains("cov.MAR.c.65."))
AD.MAR.c.cov.65 <- d.MAR.c.cov.65 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.c.cov.65.AD.women.cl.15.25 <- AD.MAR.c.cov.65[,1]
vector.mean.MAR.c.cov.65.AD.men.cl.15.25 <- AD.MAR.c.cov.65[,2]

vector.mean.MAR.c.cov.65.AD.women.cl.25.40 <- AD.MAR.c.cov.65[,3]
vector.mean.MAR.c.cov.65.AD.men.cl.25.40 <- AD.MAR.c.cov.65[,4]

vector.mean.MAR.c.cov.65.AD.women.cl.40.50 <- AD.MAR.c.cov.65[,5]
vector.mean.MAR.c.cov.65.AD.men.cl.40.50 <- AD.MAR.c.cov.65[,6]

# Median

vector.med.MAR.c.cov.65.AD.women.cl.15.25 <- AD.MAR.c.cov.65[,7]
vector.med.MAR.c.cov.65.AD.men.cl.15.25 <- AD.MAR.c.cov.65[,8]

vector.med.MAR.c.cov.65.AD.women.cl.25.40 <- AD.MAR.c.cov.65[,9]
vector.med.MAR.c.cov.65.AD.men.cl.25.40 <- AD.MAR.c.cov.65[,10]

vector.med.MAR.c.cov.65.AD.women.cl.40.50 <- AD.MAR.c.cov.65[,11]
vector.med.MAR.c.cov.65.AD.men.cl.40.50 <- AD.MAR.c.cov.65[,12]

# Standard deviation

vector.sd.MAR.c.cov.65.AD.women.cl.15.25 <- AD.MAR.c.cov.65[,13]
vector.sd.MAR.c.cov.65.AD.men.cl.15.25 <- AD.MAR.c.cov.65[,14]

vector.sd.MAR.c.cov.65.AD.women.cl.25.40 <- AD.MAR.c.cov.65[,15]
vector.sd.MAR.c.cov.65.AD.men.cl.25.40 <- AD.MAR.c.cov.65[,16]

vector.sd.MAR.c.cov.65.AD.women.cl.40.50 <- AD.MAR.c.cov.65[,17]
vector.sd.MAR.c.cov.65.AD.men.cl.40.50 <- AD.MAR.c.cov.65[,18]



# 70



d.MAR.c.cov.70 <- d.MAR %>%
  select(contains("cov.MAR.c.70."))
AD.MAR.c.cov.70 <- d.MAR.c.cov.70 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.c.cov.70.AD.women.cl.15.25 <- AD.MAR.c.cov.70[,1]
vector.mean.MAR.c.cov.70.AD.men.cl.15.25 <- AD.MAR.c.cov.70[,2]

vector.mean.MAR.c.cov.70.AD.women.cl.25.40 <- AD.MAR.c.cov.70[,3]
vector.mean.MAR.c.cov.70.AD.men.cl.25.40 <- AD.MAR.c.cov.70[,4]

vector.mean.MAR.c.cov.70.AD.women.cl.40.50 <- AD.MAR.c.cov.70[,5]
vector.mean.MAR.c.cov.70.AD.men.cl.40.50 <- AD.MAR.c.cov.70[,6]

# Median

vector.med.MAR.c.cov.70.AD.women.cl.15.25 <- AD.MAR.c.cov.70[,7]
vector.med.MAR.c.cov.70.AD.men.cl.15.25 <- AD.MAR.c.cov.70[,8]

vector.med.MAR.c.cov.70.AD.women.cl.25.40 <- AD.MAR.c.cov.70[,9]
vector.med.MAR.c.cov.70.AD.men.cl.25.40 <- AD.MAR.c.cov.70[,10]

vector.med.MAR.c.cov.70.AD.women.cl.40.50 <- AD.MAR.c.cov.70[,11]
vector.med.MAR.c.cov.70.AD.men.cl.40.50 <- AD.MAR.c.cov.70[,12]

# Standard deviation

vector.sd.MAR.c.cov.70.AD.women.cl.15.25 <- AD.MAR.c.cov.70[,13]
vector.sd.MAR.c.cov.70.AD.men.cl.15.25 <- AD.MAR.c.cov.70[,14]

vector.sd.MAR.c.cov.70.AD.women.cl.25.40 <- AD.MAR.c.cov.70[,15]
vector.sd.MAR.c.cov.70.AD.men.cl.25.40 <- AD.MAR.c.cov.70[,16]

vector.sd.MAR.c.cov.70.AD.women.cl.40.50 <- AD.MAR.c.cov.70[,17]
vector.sd.MAR.c.cov.70.AD.men.cl.40.50 <- AD.MAR.c.cov.70[,18]



# 75




d.MAR.c.cov.75 <- d.MAR %>%
  select(contains("cov.MAR.c.75."))
AD.MAR.c.cov.75 <- d.MAR.c.cov.75 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.c.cov.75.AD.women.cl.15.25 <- AD.MAR.c.cov.75[,1]
vector.mean.MAR.c.cov.75.AD.men.cl.15.25 <- AD.MAR.c.cov.75[,2]

vector.mean.MAR.c.cov.75.AD.women.cl.25.40 <- AD.MAR.c.cov.75[,3]
vector.mean.MAR.c.cov.75.AD.men.cl.25.40 <- AD.MAR.c.cov.75[,4]

vector.mean.MAR.c.cov.75.AD.women.cl.40.50 <- AD.MAR.c.cov.75[,5]
vector.mean.MAR.c.cov.75.AD.men.cl.40.50 <- AD.MAR.c.cov.75[,6]

# Median

vector.med.MAR.c.cov.75.AD.women.cl.15.25 <- AD.MAR.c.cov.75[,7]
vector.med.MAR.c.cov.75.AD.men.cl.15.25 <- AD.MAR.c.cov.75[,8]

vector.med.MAR.c.cov.75.AD.women.cl.25.40 <- AD.MAR.c.cov.75[,9]
vector.med.MAR.c.cov.75.AD.men.cl.25.40 <- AD.MAR.c.cov.75[,10]

vector.med.MAR.c.cov.75.AD.women.cl.40.50 <- AD.MAR.c.cov.75[,11]
vector.med.MAR.c.cov.75.AD.men.cl.40.50 <- AD.MAR.c.cov.75[,12]

# Standard deviation

vector.sd.MAR.c.cov.75.AD.women.cl.15.25 <- AD.MAR.c.cov.75[,13]
vector.sd.MAR.c.cov.75.AD.men.cl.15.25 <- AD.MAR.c.cov.75[,14]

vector.sd.MAR.c.cov.75.AD.women.cl.25.40 <- AD.MAR.c.cov.75[,15]
vector.sd.MAR.c.cov.75.AD.men.cl.25.40 <- AD.MAR.c.cov.75[,16]

vector.sd.MAR.c.cov.75.AD.women.cl.40.50 <- AD.MAR.c.cov.75[,17]
vector.sd.MAR.c.cov.75.AD.men.cl.40.50 <- AD.MAR.c.cov.75[,18]



# 80



d.MAR.c.cov.80 <- d.MAR %>%
  select(contains("cov.MAR.c.80."))
AD.MAR.c.cov.80 <- d.MAR.c.cov.80 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.c.cov.80.AD.women.cl.15.25 <- AD.MAR.c.cov.80[,1]
vector.mean.MAR.c.cov.80.AD.men.cl.15.25 <- AD.MAR.c.cov.80[,2]

vector.mean.MAR.c.cov.80.AD.women.cl.25.40 <- AD.MAR.c.cov.80[,3]
vector.mean.MAR.c.cov.80.AD.men.cl.25.40 <- AD.MAR.c.cov.80[,4]

vector.mean.MAR.c.cov.80.AD.women.cl.40.50 <- AD.MAR.c.cov.80[,5]
vector.mean.MAR.c.cov.80.AD.men.cl.40.50 <- AD.MAR.c.cov.80[,6]

# Median

vector.med.MAR.c.cov.80.AD.women.cl.15.25 <- AD.MAR.c.cov.80[,7]
vector.med.MAR.c.cov.80.AD.men.cl.15.25 <- AD.MAR.c.cov.80[,8]

vector.med.MAR.c.cov.80.AD.women.cl.25.40 <- AD.MAR.c.cov.80[,9]
vector.med.MAR.c.cov.80.AD.men.cl.25.40 <- AD.MAR.c.cov.80[,10]

vector.med.MAR.c.cov.80.AD.women.cl.40.50 <- AD.MAR.c.cov.80[,11]
vector.med.MAR.c.cov.80.AD.men.cl.40.50 <- AD.MAR.c.cov.80[,12]

# Standard deviation

vector.sd.MAR.c.cov.80.AD.women.cl.15.25 <- AD.MAR.c.cov.80[,13]
vector.sd.MAR.c.cov.80.AD.men.cl.15.25 <- AD.MAR.c.cov.80[,14]

vector.sd.MAR.c.cov.80.AD.women.cl.25.40 <- AD.MAR.c.cov.80[,15]
vector.sd.MAR.c.cov.80.AD.men.cl.25.40 <- AD.MAR.c.cov.80[,16]

vector.sd.MAR.c.cov.80.AD.women.cl.40.50 <- AD.MAR.c.cov.80[,17]
vector.sd.MAR.c.cov.80.AD.men.cl.40.50 <- AD.MAR.c.cov.80[,18]



# 85



d.MAR.c.cov.85 <- d.MAR %>%
  select(contains("cov.MAR.c.85."))
AD.MAR.c.cov.85 <- d.MAR.c.cov.85 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.c.cov.85.AD.women.cl.15.25 <- AD.MAR.c.cov.85[,1]
vector.mean.MAR.c.cov.85.AD.men.cl.15.25 <- AD.MAR.c.cov.85[,2]

vector.mean.MAR.c.cov.85.AD.women.cl.25.40 <- AD.MAR.c.cov.85[,3]
vector.mean.MAR.c.cov.85.AD.men.cl.25.40 <- AD.MAR.c.cov.85[,4]

vector.mean.MAR.c.cov.85.AD.women.cl.40.50 <- AD.MAR.c.cov.85[,5]
vector.mean.MAR.c.cov.85.AD.men.cl.40.50 <- AD.MAR.c.cov.85[,6]

# Median

vector.med.MAR.c.cov.85.AD.women.cl.15.25 <- AD.MAR.c.cov.85[,7]
vector.med.MAR.c.cov.85.AD.men.cl.15.25 <- AD.MAR.c.cov.85[,8]

vector.med.MAR.c.cov.85.AD.women.cl.25.40 <- AD.MAR.c.cov.85[,9]
vector.med.MAR.c.cov.85.AD.men.cl.25.40 <- AD.MAR.c.cov.85[,10]

vector.med.MAR.c.cov.85.AD.women.cl.40.50 <- AD.MAR.c.cov.85[,11]
vector.med.MAR.c.cov.85.AD.men.cl.40.50 <- AD.MAR.c.cov.85[,12]

# Standard deviation

vector.sd.MAR.c.cov.85.AD.women.cl.15.25 <- AD.MAR.c.cov.85[,13]
vector.sd.MAR.c.cov.85.AD.men.cl.15.25 <- AD.MAR.c.cov.85[,14]

vector.sd.MAR.c.cov.85.AD.women.cl.25.40 <- AD.MAR.c.cov.85[,15]
vector.sd.MAR.c.cov.85.AD.men.cl.25.40 <- AD.MAR.c.cov.85[,16]

vector.sd.MAR.c.cov.85.AD.women.cl.40.50 <- AD.MAR.c.cov.85[,17]
vector.sd.MAR.c.cov.85.AD.men.cl.40.50 <- AD.MAR.c.cov.85[,18]



# 90



d.MAR.c.cov.90 <- d.MAR %>%
  select(contains("cov.MAR.c.90."))
AD.MAR.c.cov.90 <- d.MAR.c.cov.90 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.c.cov.90.AD.women.cl.15.25 <- AD.MAR.c.cov.90[,1]
vector.mean.MAR.c.cov.90.AD.men.cl.15.25 <- AD.MAR.c.cov.90[,2]

vector.mean.MAR.c.cov.90.AD.women.cl.25.40 <- AD.MAR.c.cov.90[,3]
vector.mean.MAR.c.cov.90.AD.men.cl.25.40 <- AD.MAR.c.cov.90[,4]

vector.mean.MAR.c.cov.90.AD.women.cl.40.50 <- AD.MAR.c.cov.90[,5]
vector.mean.MAR.c.cov.90.AD.men.cl.40.50 <- AD.MAR.c.cov.90[,6]

# Median

vector.med.MAR.c.cov.90.AD.women.cl.15.25 <- AD.MAR.c.cov.90[,7]
vector.med.MAR.c.cov.90.AD.men.cl.15.25 <- AD.MAR.c.cov.90[,8]

vector.med.MAR.c.cov.90.AD.women.cl.25.40 <- AD.MAR.c.cov.90[,9]
vector.med.MAR.c.cov.90.AD.men.cl.25.40 <- AD.MAR.c.cov.90[,10]

vector.med.MAR.c.cov.90.AD.women.cl.40.50 <- AD.MAR.c.cov.90[,11]
vector.med.MAR.c.cov.90.AD.men.cl.40.50 <- AD.MAR.c.cov.90[,12]

# Standard deviation

vector.sd.MAR.c.cov.90.AD.women.cl.15.25 <- AD.MAR.c.cov.90[,13]
vector.sd.MAR.c.cov.90.AD.men.cl.15.25 <- AD.MAR.c.cov.90[,14]

vector.sd.MAR.c.cov.90.AD.women.cl.25.40 <- AD.MAR.c.cov.90[,15]
vector.sd.MAR.c.cov.90.AD.men.cl.25.40 <- AD.MAR.c.cov.90[,16]

vector.sd.MAR.c.cov.90.AD.women.cl.40.50 <- AD.MAR.c.cov.90[,17]
vector.sd.MAR.c.cov.90.AD.men.cl.40.50 <- AD.MAR.c.cov.90[,18]



# 95



d.MAR.c.cov.95 <- d.MAR %>%
  select(contains("cov.MAR.c.95."))
AD.MAR.c.cov.95 <- d.MAR.c.cov.95 %>%
  select(contains(".AD.")) 



# Mean

vector.mean.MAR.c.cov.95.AD.women.cl.15.25 <- AD.MAR.c.cov.95[,1]
vector.mean.MAR.c.cov.95.AD.men.cl.15.25 <- AD.MAR.c.cov.95[,2]

vector.mean.MAR.c.cov.95.AD.women.cl.25.40 <- AD.MAR.c.cov.95[,3]
vector.mean.MAR.c.cov.95.AD.men.cl.25.40 <- AD.MAR.c.cov.95[,4]

vector.mean.MAR.c.cov.95.AD.women.cl.40.50 <- AD.MAR.c.cov.95[,5]
vector.mean.MAR.c.cov.95.AD.men.cl.40.50 <- AD.MAR.c.cov.95[,6]

# Median

vector.med.MAR.c.cov.95.AD.women.cl.15.25 <- AD.MAR.c.cov.95[,7]
vector.med.MAR.c.cov.95.AD.men.cl.15.25 <- AD.MAR.c.cov.95[,8]

vector.med.MAR.c.cov.95.AD.women.cl.25.40 <- AD.MAR.c.cov.95[,9]
vector.med.MAR.c.cov.95.AD.men.cl.25.40 <- AD.MAR.c.cov.95[,10]

vector.med.MAR.c.cov.95.AD.women.cl.40.50 <- AD.MAR.c.cov.95[,11]
vector.med.MAR.c.cov.95.AD.men.cl.40.50 <- AD.MAR.c.cov.95[,12]

# Standard deviation

vector.sd.MAR.c.cov.95.AD.women.cl.15.25 <- AD.MAR.c.cov.95[,13]
vector.sd.MAR.c.cov.95.AD.men.cl.15.25 <- AD.MAR.c.cov.95[,14]

vector.sd.MAR.c.cov.95.AD.women.cl.25.40 <- AD.MAR.c.cov.95[,15]
vector.sd.MAR.c.cov.95.AD.men.cl.25.40 <- AD.MAR.c.cov.95[,16]

vector.sd.MAR.c.cov.95.AD.women.cl.40.50 <- AD.MAR.c.cov.95[,17]
vector.sd.MAR.c.cov.95.AD.men.cl.40.50 <- AD.MAR.c.cov.95[,18]




## Computations ------


# Proportions MCAR and MAR (a) -------


mcar_mar_a_cov.35_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.a.cov.35.cl.prop.men15.25.F.15.25)

mcar_mar_a_cov.35_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.a.cov.35.cl.prop.women15.25.M.15.25)

mcar_mar_a_cov.35_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.a.cov.35.cl.prop.men25.40.F.15.25)

mcar_mar_a_cov.35_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.a.cov.35.cl.prop.women15.25.M.25.40)

mcar_mar_a_cov.35_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.a.cov.35.cl.prop.men25.40.F.25.40)

mcar_mar_a_cov.35_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.a.cov.35.cl.prop.women25.40.M.25.40)

mcar_mar_a_cov.35_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.a.cov.35.cl.prop.men40.50.F.15.25)


mcar_mar_a_cov.35_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.a.cov.35.cl.prop.women15.25.M.40.50)

mcar_mar_a_cov.35_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.a.cov.35.cl.prop.men40.50.F.25.40)

mcar_mar_a_cov.35_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.a.cov.35.cl.prop.women25.40.M.40.50)



# 40 - a


mcar_mar_a_cov.40_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.a.cov.40.cl.prop.men15.25.F.15.25)


mcar_mar_a_cov.40_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.a.cov.40.cl.prop.women15.25.M.15.25)


mcar_mar_a_cov.40_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.a.cov.40.cl.prop.men25.40.F.15.25)



mcar_mar_a_cov.40_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.a.cov.40.cl.prop.women15.25.M.25.40)

mcar_mar_a_cov.40_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.a.cov.40.cl.prop.men25.40.F.25.40)


mcar_mar_a_cov.40_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.a.cov.40.cl.prop.women25.40.M.25.40)


mcar_mar_a_cov.40_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.a.cov.40.cl.prop.men40.50.F.15.25)


mcar_mar_a_cov.40_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.a.cov.40.cl.prop.women15.25.M.40.50)


mcar_mar_a_cov.40_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.a.cov.40.cl.prop.men40.50.F.25.40)


mcar_mar_a_cov.40_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.a.cov.40.cl.prop.women25.40.M.40.50)



# 45 - a 


mcar_mar_a_cov.45_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.a.cov.45.cl.prop.men15.25.F.15.25)


mcar_mar_a_cov.45_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.a.cov.45.cl.prop.women15.25.M.15.25)


mcar_mar_a_cov.45_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.a.cov.45.cl.prop.men25.40.F.15.25)



mcar_mar_a_cov.45_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.a.cov.45.cl.prop.women15.25.M.25.40)

mcar_mar_a_cov.45_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.a.cov.45.cl.prop.men25.40.F.25.40)


mcar_mar_a_cov.45_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.a.cov.45.cl.prop.women25.40.M.25.40)


mcar_mar_a_cov.45_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.a.cov.45.cl.prop.men40.50.F.15.25)


mcar_mar_a_cov.45_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.a.cov.45.cl.prop.women15.25.M.40.50)


mcar_mar_a_cov.45_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.a.cov.45.cl.prop.men40.50.F.25.40)


mcar_mar_a_cov.45_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.a.cov.45.cl.prop.women25.40.M.40.50)



# 50 - a


mcar_mar_a_cov.50_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.a.cov.50.cl.prop.men15.25.F.15.25)


mcar_mar_a_cov.50_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.a.cov.50.cl.prop.women15.25.M.15.25)


mcar_mar_a_cov.50_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.a.cov.50.cl.prop.men25.40.F.15.25)



mcar_mar_a_cov.50_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.a.cov.50.cl.prop.women15.25.M.25.40)

mcar_mar_a_cov.50_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.a.cov.50.cl.prop.men25.40.F.25.40)


mcar_mar_a_cov.50_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.a.cov.50.cl.prop.women25.40.M.25.40)


mcar_mar_a_cov.50_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.a.cov.50.cl.prop.men40.50.F.15.25)


mcar_mar_a_cov.50_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.a.cov.50.cl.prop.women15.25.M.40.50)


mcar_mar_a_cov.50_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.a.cov.50.cl.prop.men40.50.F.25.40)


mcar_mar_a_cov.50_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.a.cov.50.cl.prop.women25.40.M.40.50)



# 55 - a


mcar_mar_a_cov.55_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.a.cov.55.cl.prop.men15.25.F.15.25)


mcar_mar_a_cov.55_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.a.cov.55.cl.prop.women15.25.M.15.25)


mcar_mar_a_cov.55_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.a.cov.55.cl.prop.men25.40.F.15.25)



mcar_mar_a_cov.55_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.a.cov.55.cl.prop.women15.25.M.25.40)

mcar_mar_a_cov.55_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.a.cov.55.cl.prop.men25.40.F.25.40)


mcar_mar_a_cov.55_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.a.cov.55.cl.prop.women25.40.M.25.40)


mcar_mar_a_cov.55_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.a.cov.55.cl.prop.men40.50.F.15.25)


mcar_mar_a_cov.55_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.a.cov.55.cl.prop.women15.25.M.40.50)


mcar_mar_a_cov.55_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.a.cov.55.cl.prop.men40.50.F.25.40)


mcar_mar_a_cov.55_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.a.cov.55.cl.prop.women25.40.M.40.50)



# 60 - a 


mcar_mar_a_cov.60_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.a.cov.60.cl.prop.men15.25.F.15.25)


mcar_mar_a_cov.60_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.a.cov.60.cl.prop.women15.25.M.15.25)


mcar_mar_a_cov.60_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.a.cov.60.cl.prop.men25.40.F.15.25)



mcar_mar_a_cov.60_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.a.cov.60.cl.prop.women15.25.M.25.40)

mcar_mar_a_cov.60_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.a.cov.60.cl.prop.men25.40.F.25.40)


mcar_mar_a_cov.60_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.a.cov.60.cl.prop.women25.40.M.25.40)


mcar_mar_a_cov.60_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.a.cov.60.cl.prop.men40.50.F.15.25)

mcar_mar_a_cov.60_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.a.cov.60.cl.prop.women15.25.M.40.50)


mcar_mar_a_cov.60_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.a.cov.60.cl.prop.men40.50.F.25.40)


mcar_mar_a_cov.60_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.a.cov.60.cl.prop.women25.40.M.40.50)



# 65 - a


mcar_mar_a_cov.65_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.a.cov.65.cl.prop.men15.25.F.15.25)


mcar_mar_a_cov.65_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.a.cov.65.cl.prop.women15.25.M.15.25)


mcar_mar_a_cov.65_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.a.cov.65.cl.prop.men25.40.F.15.25)



mcar_mar_a_cov.65_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.a.cov.65.cl.prop.women15.25.M.25.40)

mcar_mar_a_cov.65_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.a.cov.65.cl.prop.men25.40.F.25.40)


mcar_mar_a_cov.65_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.a.cov.65.cl.prop.women25.40.M.25.40)


mcar_mar_a_cov.65_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.a.cov.65.cl.prop.men40.50.F.15.25)

mcar_mar_a_cov.65_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.a.cov.65.cl.prop.women15.25.M.40.50)


mcar_mar_a_cov.65_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.a.cov.65.cl.prop.men40.50.F.25.40)


mcar_mar_a_cov.65_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.a.cov.65.cl.prop.women25.40.M.40.50)



# 65 - a


mcar_mar_a_cov.70_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.a.cov.70.cl.prop.men15.25.F.15.25)


mcar_mar_a_cov.70_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.a.cov.70.cl.prop.women15.25.M.15.25)


mcar_mar_a_cov.70_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.a.cov.70.cl.prop.men25.40.F.15.25)



mcar_mar_a_cov.70_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.a.cov.70.cl.prop.women15.25.M.25.40)

mcar_mar_a_cov.70_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.a.cov.70.cl.prop.men25.40.F.25.40)


mcar_mar_a_cov.70_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.a.cov.70.cl.prop.women25.40.M.25.40)


mcar_mar_a_cov.70_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.a.cov.70.cl.prop.men40.50.F.15.25)

mcar_mar_a_cov.70_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.a.cov.70.cl.prop.women15.25.M.40.50)

mcar_mar_a_cov.70_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.a.cov.70.cl.prop.men40.50.F.25.40)


mcar_mar_a_cov.70_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.a.cov.70.cl.prop.women25.40.M.40.50)



# 75 - a


mcar_mar_a_cov.75_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.a.cov.75.cl.prop.men15.25.F.15.25)


mcar_mar_a_cov.75_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.a.cov.75.cl.prop.women15.25.M.15.25)


mcar_mar_a_cov.75_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.a.cov.75.cl.prop.men25.40.F.15.25)



mcar_mar_a_cov.75_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.a.cov.75.cl.prop.women15.25.M.25.40)

mcar_mar_a_cov.75_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.a.cov.75.cl.prop.men25.40.F.25.40)


mcar_mar_a_cov.75_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.a.cov.75.cl.prop.women25.40.M.25.40)


mcar_mar_a_cov.75_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.a.cov.75.cl.prop.men40.50.F.15.25)


mcar_mar_a_cov.75_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.a.cov.75.cl.prop.women15.25.M.40.50)


mcar_mar_a_cov.75_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.a.cov.75.cl.prop.men40.50.F.25.40)


mcar_mar_a_cov.75_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.a.cov.75.cl.prop.women25.40.M.40.50)


# 80 - a


mcar_mar_a_cov.80_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.a.cov.80.cl.prop.men15.25.F.15.25)


mcar_mar_a_cov.80_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.a.cov.80.cl.prop.women15.25.M.15.25)


mcar_mar_a_cov.80_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.a.cov.80.cl.prop.men25.40.F.15.25)



mcar_mar_a_cov.80_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.a.cov.80.cl.prop.women15.25.M.25.40)

mcar_mar_a_cov.80_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.a.cov.80.cl.prop.men25.40.F.25.40)


mcar_mar_a_cov.80_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.a.cov.80.cl.prop.women25.40.M.25.40)


mcar_mar_a_cov.80_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.a.cov.80.cl.prop.men40.50.F.15.25)


mcar_mar_a_cov.80_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.a.cov.80.cl.prop.women15.25.M.40.50)


mcar_mar_a_cov.80_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.a.cov.80.cl.prop.men40.50.F.25.40)


mcar_mar_a_cov.80_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.a.cov.80.cl.prop.women25.40.M.40.50)



# 85 - a


mcar_mar_a_cov.85_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.a.cov.85.cl.prop.men15.25.F.15.25)


mcar_mar_a_cov.85_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.a.cov.85.cl.prop.women15.25.M.15.25)


mcar_mar_a_cov.85_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.a.cov.85.cl.prop.men25.40.F.15.25)



mcar_mar_a_cov.85_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.a.cov.85.cl.prop.women15.25.M.25.40)

mcar_mar_a_cov.85_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.a.cov.85.cl.prop.men25.40.F.25.40)


mcar_mar_a_cov.85_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.a.cov.85.cl.prop.women25.40.M.25.40)


mcar_mar_a_cov.85_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.a.cov.85.cl.prop.men40.50.F.15.25)


mcar_mar_a_cov.85_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.a.cov.85.cl.prop.women15.25.M.40.50)


mcar_mar_a_cov.85_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.a.cov.85.cl.prop.men40.50.F.25.40)


mcar_mar_a_cov.85_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.a.cov.85.cl.prop.women25.40.M.40.50)



# 90 - a


mcar_mar_a_cov.90_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.a.cov.90.cl.prop.men15.25.F.15.25)


mcar_mar_a_cov.90_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.a.cov.90.cl.prop.women15.25.M.15.25)


mcar_mar_a_cov.90_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.a.cov.90.cl.prop.men25.40.F.15.25)



mcar_mar_a_cov.90_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.a.cov.90.cl.prop.women15.25.M.25.40)

mcar_mar_a_cov.90_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.a.cov.90.cl.prop.men25.40.F.25.40)


mcar_mar_a_cov.90_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.a.cov.90.cl.prop.women25.40.M.25.40)


mcar_mar_a_cov.90_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.a.cov.90.cl.prop.men40.50.F.15.25)


mcar_mar_a_cov.90_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.a.cov.90.cl.prop.women15.25.M.40.50)


mcar_mar_a_cov.90_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.a.cov.90.cl.prop.men40.50.F.25.40)


mcar_mar_a_cov.90_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.a.cov.90.cl.prop.women25.40.M.40.50)



# 95 - a


mcar_mar_a_cov.95_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.a.cov.95.cl.prop.men15.25.F.15.25)


mcar_mar_a_cov.95_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.a.cov.95.cl.prop.women15.25.M.15.25)


mcar_mar_a_cov.95_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.a.cov.95.cl.prop.men25.40.F.15.25)



mcar_mar_a_cov.95_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.a.cov.95.cl.prop.women15.25.M.25.40)

mcar_mar_a_cov.95_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.a.cov.95.cl.prop.men25.40.F.25.40)


mcar_mar_a_cov.95_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.a.cov.95.cl.prop.women25.40.M.25.40)


mcar_mar_a_cov.95_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.a.cov.95.cl.prop.men40.50.F.15.25)


mcar_mar_a_cov.95_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.a.cov.95.cl.prop.women15.25.M.40.50)


mcar_mar_a_cov.95_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.a.cov.95.cl.prop.men40.50.F.25.40)


mcar_mar_a_cov.95_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.a.cov.95.cl.prop.women25.40.M.40.50)



# Proportions MCAR and MAR (b) -------


mcar_mar_b_cov.35_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.b.cov.35.cl.prop.men15.25.F.15.25)

mcar_mar_b_cov.35_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.b.cov.35.cl.prop.women15.25.M.15.25)

mcar_mar_b_cov.35_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.b.cov.35.cl.prop.men25.40.F.15.25)

mcar_mar_b_cov.35_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.b.cov.35.cl.prop.women15.25.M.25.40)

mcar_mar_b_cov.35_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.b.cov.35.cl.prop.men25.40.F.25.40)

mcar_mar_b_cov.35_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.b.cov.35.cl.prop.women25.40.M.25.40)

mcar_mar_b_cov.35_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.b.cov.35.cl.prop.men40.50.F.15.25)


mcar_mar_b_cov.35_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.b.cov.35.cl.prop.women15.25.M.40.50)

mcar_mar_b_cov.35_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.b.cov.35.cl.prop.men40.50.F.25.40)

mcar_mar_b_cov.35_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.b.cov.35.cl.prop.women25.40.M.40.50)



# 40 - b


mcar_mar_b_cov.40_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.b.cov.40.cl.prop.men15.25.F.15.25)


mcar_mar_b_cov.40_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.b.cov.40.cl.prop.women15.25.M.15.25)


mcar_mar_b_cov.40_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.b.cov.40.cl.prop.men25.40.F.15.25)



mcar_mar_b_cov.40_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.b.cov.40.cl.prop.women15.25.M.25.40)

mcar_mar_b_cov.40_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.b.cov.40.cl.prop.men25.40.F.25.40)


mcar_mar_b_cov.40_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.b.cov.40.cl.prop.women25.40.M.25.40)


mcar_mar_b_cov.40_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.b.cov.40.cl.prop.men40.50.F.15.25)


mcar_mar_b_cov.40_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.b.cov.40.cl.prop.women15.25.M.40.50)


mcar_mar_b_cov.40_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.b.cov.40.cl.prop.men40.50.F.25.40)


mcar_mar_b_cov.40_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.b.cov.40.cl.prop.women25.40.M.40.50)



# 45 - b 


mcar_mar_b_cov.45_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.b.cov.45.cl.prop.men15.25.F.15.25)


mcar_mar_b_cov.45_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.b.cov.45.cl.prop.women15.25.M.15.25)


mcar_mar_b_cov.45_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.b.cov.45.cl.prop.men25.40.F.15.25)



mcar_mar_b_cov.45_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.b.cov.45.cl.prop.women15.25.M.25.40)

mcar_mar_b_cov.45_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.b.cov.45.cl.prop.men25.40.F.25.40)


mcar_mar_b_cov.45_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.b.cov.45.cl.prop.women25.40.M.25.40)


mcar_mar_b_cov.45_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.b.cov.45.cl.prop.men40.50.F.15.25)


mcar_mar_b_cov.45_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.b.cov.45.cl.prop.women15.25.M.40.50)


mcar_mar_b_cov.45_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.b.cov.45.cl.prop.men40.50.F.25.40)


mcar_mar_b_cov.45_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.b.cov.45.cl.prop.women25.40.M.40.50)



# 50 - b


mcar_mar_b_cov.50_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.b.cov.50.cl.prop.men15.25.F.15.25)


mcar_mar_b_cov.50_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.b.cov.50.cl.prop.women15.25.M.15.25)


mcar_mar_b_cov.50_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.b.cov.50.cl.prop.men25.40.F.15.25)



mcar_mar_b_cov.50_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.b.cov.50.cl.prop.women15.25.M.25.40)

mcar_mar_b_cov.50_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.b.cov.50.cl.prop.men25.40.F.25.40)


mcar_mar_b_cov.50_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.b.cov.50.cl.prop.women25.40.M.25.40)


mcar_mar_b_cov.50_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.b.cov.50.cl.prop.men40.50.F.15.25)


mcar_mar_b_cov.50_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.b.cov.50.cl.prop.women15.25.M.40.50)


mcar_mar_b_cov.50_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.b.cov.50.cl.prop.men40.50.F.25.40)


mcar_mar_b_cov.50_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.b.cov.50.cl.prop.women25.40.M.40.50)



# 55 - b


mcar_mar_b_cov.55_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.b.cov.55.cl.prop.men15.25.F.15.25)


mcar_mar_b_cov.55_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.b.cov.55.cl.prop.women15.25.M.15.25)


mcar_mar_b_cov.55_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.b.cov.55.cl.prop.men25.40.F.15.25)



mcar_mar_b_cov.55_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.b.cov.55.cl.prop.women15.25.M.25.40)

mcar_mar_b_cov.55_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.b.cov.55.cl.prop.men25.40.F.25.40)


mcar_mar_b_cov.55_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.b.cov.55.cl.prop.women25.40.M.25.40)


mcar_mar_b_cov.55_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.b.cov.55.cl.prop.men40.50.F.15.25)


mcar_mar_b_cov.55_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.b.cov.55.cl.prop.women15.25.M.40.50)


mcar_mar_b_cov.55_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.b.cov.55.cl.prop.men40.50.F.25.40)


mcar_mar_b_cov.55_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.b.cov.55.cl.prop.women25.40.M.40.50)



# 60 - b 


mcar_mar_b_cov.60_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.b.cov.60.cl.prop.men15.25.F.15.25)


mcar_mar_b_cov.60_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.b.cov.60.cl.prop.women15.25.M.15.25)


mcar_mar_b_cov.60_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.b.cov.60.cl.prop.men25.40.F.15.25)



mcar_mar_b_cov.60_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.b.cov.60.cl.prop.women15.25.M.25.40)

mcar_mar_b_cov.60_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.b.cov.60.cl.prop.men25.40.F.25.40)


mcar_mar_b_cov.60_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.b.cov.60.cl.prop.women25.40.M.25.40)


mcar_mar_b_cov.60_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.b.cov.60.cl.prop.men40.50.F.15.25)

mcar_mar_b_cov.60_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.b.cov.60.cl.prop.women15.25.M.40.50)


mcar_mar_b_cov.60_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.b.cov.60.cl.prop.men40.50.F.25.40)


mcar_mar_b_cov.60_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.b.cov.60.cl.prop.women25.40.M.40.50)



# 65 - b


mcar_mar_b_cov.65_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.b.cov.65.cl.prop.men15.25.F.15.25)


mcar_mar_b_cov.65_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.b.cov.65.cl.prop.women15.25.M.15.25)


mcar_mar_b_cov.65_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.b.cov.65.cl.prop.men25.40.F.15.25)



mcar_mar_b_cov.65_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.b.cov.65.cl.prop.women15.25.M.25.40)

mcar_mar_b_cov.65_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.b.cov.65.cl.prop.men25.40.F.25.40)


mcar_mar_b_cov.65_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.b.cov.65.cl.prop.women25.40.M.25.40)


mcar_mar_b_cov.65_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.b.cov.65.cl.prop.men40.50.F.15.25)

mcar_mar_b_cov.65_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.b.cov.65.cl.prop.women15.25.M.40.50)


mcar_mar_b_cov.65_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.b.cov.65.cl.prop.men40.50.F.25.40)


mcar_mar_b_cov.65_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.b.cov.65.cl.prop.women25.40.M.40.50)



# 65 - b


mcar_mar_b_cov.70_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.b.cov.70.cl.prop.men15.25.F.15.25)


mcar_mar_b_cov.70_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.b.cov.70.cl.prop.women15.25.M.15.25)


mcar_mar_b_cov.70_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.b.cov.70.cl.prop.men25.40.F.15.25)



mcar_mar_b_cov.70_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.b.cov.70.cl.prop.women15.25.M.25.40)

mcar_mar_b_cov.70_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.b.cov.70.cl.prop.men25.40.F.25.40)


mcar_mar_b_cov.70_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.b.cov.70.cl.prop.women25.40.M.25.40)


mcar_mar_b_cov.70_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.b.cov.70.cl.prop.men40.50.F.15.25)

mcar_mar_b_cov.70_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.b.cov.70.cl.prop.women15.25.M.40.50)

mcar_mar_b_cov.70_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.b.cov.70.cl.prop.men40.50.F.25.40)


mcar_mar_b_cov.70_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.b.cov.70.cl.prop.women25.40.M.40.50)



# 75 - b


mcar_mar_b_cov.75_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.b.cov.75.cl.prop.men15.25.F.15.25)


mcar_mar_b_cov.75_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.b.cov.75.cl.prop.women15.25.M.15.25)


mcar_mar_b_cov.75_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.b.cov.75.cl.prop.men25.40.F.15.25)



mcar_mar_b_cov.75_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.b.cov.75.cl.prop.women15.25.M.25.40)

mcar_mar_b_cov.75_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.b.cov.75.cl.prop.men25.40.F.25.40)


mcar_mar_b_cov.75_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.b.cov.75.cl.prop.women25.40.M.25.40)


mcar_mar_b_cov.75_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.b.cov.75.cl.prop.men40.50.F.15.25)


mcar_mar_b_cov.75_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.b.cov.75.cl.prop.women15.25.M.40.50)


mcar_mar_b_cov.75_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.b.cov.75.cl.prop.men40.50.F.25.40)


mcar_mar_b_cov.75_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.b.cov.75.cl.prop.women25.40.M.40.50)


# 80 - b


mcar_mar_b_cov.80_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.b.cov.80.cl.prop.men15.25.F.15.25)


mcar_mar_b_cov.80_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.b.cov.80.cl.prop.women15.25.M.15.25)


mcar_mar_b_cov.80_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.b.cov.80.cl.prop.men25.40.F.15.25)



mcar_mar_b_cov.80_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.b.cov.80.cl.prop.women15.25.M.25.40)

mcar_mar_b_cov.80_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.b.cov.80.cl.prop.men25.40.F.25.40)


mcar_mar_b_cov.80_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.b.cov.80.cl.prop.women25.40.M.25.40)


mcar_mar_b_cov.80_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.b.cov.80.cl.prop.men40.50.F.15.25)


mcar_mar_b_cov.80_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.b.cov.80.cl.prop.women15.25.M.40.50)


mcar_mar_b_cov.80_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.b.cov.80.cl.prop.men40.50.F.25.40)


mcar_mar_b_cov.80_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.b.cov.80.cl.prop.women25.40.M.40.50)



# 85 - b


mcar_mar_b_cov.85_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.b.cov.85.cl.prop.men15.25.F.15.25)


mcar_mar_b_cov.85_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.b.cov.85.cl.prop.women15.25.M.15.25)


mcar_mar_b_cov.85_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.b.cov.85.cl.prop.men25.40.F.15.25)



mcar_mar_b_cov.85_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.b.cov.85.cl.prop.women15.25.M.25.40)

mcar_mar_b_cov.85_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.b.cov.85.cl.prop.men25.40.F.25.40)


mcar_mar_b_cov.85_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.b.cov.85.cl.prop.women25.40.M.25.40)


mcar_mar_b_cov.85_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.b.cov.85.cl.prop.men40.50.F.15.25)


mcar_mar_b_cov.85_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.b.cov.85.cl.prop.women15.25.M.40.50)


mcar_mar_b_cov.85_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.b.cov.85.cl.prop.men40.50.F.25.40)


mcar_mar_b_cov.85_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.b.cov.85.cl.prop.women25.40.M.40.50)



# 90 - b


mcar_mar_b_cov.90_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.b.cov.90.cl.prop.men15.25.F.15.25)


mcar_mar_b_cov.90_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.b.cov.90.cl.prop.women15.25.M.15.25)


mcar_mar_b_cov.90_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.b.cov.90.cl.prop.men25.40.F.15.25)



mcar_mar_b_cov.90_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.b.cov.90.cl.prop.women15.25.M.25.40)

mcar_mar_b_cov.90_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.b.cov.90.cl.prop.men25.40.F.25.40)


mcar_mar_b_cov.90_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.b.cov.90.cl.prop.women25.40.M.25.40)


mcar_mar_b_cov.90_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.b.cov.90.cl.prop.men40.50.F.15.25)


mcar_mar_b_cov.90_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.b.cov.90.cl.prop.women15.25.M.40.50)


mcar_mar_b_cov.90_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.b.cov.90.cl.prop.men40.50.F.25.40)


mcar_mar_b_cov.90_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.b.cov.90.cl.prop.women25.40.M.40.50)



# 95 - b


mcar_mar_b_cov.95_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.b.cov.95.cl.prop.men15.25.F.15.25)


mcar_mar_b_cov.95_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.b.cov.95.cl.prop.women15.25.M.15.25)


mcar_mar_b_cov.95_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.b.cov.95.cl.prop.men25.40.F.15.25)



mcar_mar_b_cov.95_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.b.cov.95.cl.prop.women15.25.M.25.40)

mcar_mar_b_cov.95_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.b.cov.95.cl.prop.men25.40.F.25.40)


mcar_mar_b_cov.95_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.b.cov.95.cl.prop.women25.40.M.25.40)


mcar_mar_b_cov.95_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.b.cov.95.cl.prop.men40.50.F.15.25)


mcar_mar_b_cov.95_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.b.cov.95.cl.prop.women15.25.M.40.50)


mcar_mar_b_cov.95_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.b.cov.95.cl.prop.men40.50.F.25.40)


mcar_mar_b_cov.95_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.b.cov.95.cl.prop.women25.40.M.40.50)



# Proportions MCAR and MAR (c) -------


mcar_mar_c_cov.35_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.c.cov.35.cl.prop.men15.25.F.15.25)

mcar_mar_c_cov.35_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.c.cov.35.cl.prop.women15.25.M.15.25)

mcar_mar_c_cov.35_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.c.cov.35.cl.prop.men25.40.F.15.25)

mcar_mar_c_cov.35_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.c.cov.35.cl.prop.women15.25.M.25.40)

mcar_mar_c_cov.35_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.c.cov.35.cl.prop.men25.40.F.25.40)

mcar_mar_c_cov.35_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.c.cov.35.cl.prop.women25.40.M.25.40)

mcar_mar_c_cov.35_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.c.cov.35.cl.prop.men40.50.F.15.25)


mcar_mar_c_cov.35_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.c.cov.35.cl.prop.women15.25.M.40.50)

mcar_mar_c_cov.35_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.c.cov.35.cl.prop.men40.50.F.25.40)

mcar_mar_c_cov.35_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.35.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.c.cov.35.cl.prop.women25.40.M.40.50)



# 40 - c


mcar_mar_c_cov.40_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.c.cov.40.cl.prop.men15.25.F.15.25)


mcar_mar_c_cov.40_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.c.cov.40.cl.prop.women15.25.M.15.25)


mcar_mar_c_cov.40_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.c.cov.40.cl.prop.men25.40.F.15.25)



mcar_mar_c_cov.40_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.c.cov.40.cl.prop.women15.25.M.25.40)

mcar_mar_c_cov.40_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.c.cov.40.cl.prop.men25.40.F.25.40)


mcar_mar_c_cov.40_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.c.cov.40.cl.prop.women25.40.M.25.40)


mcar_mar_c_cov.40_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.c.cov.40.cl.prop.men40.50.F.15.25)


mcar_mar_c_cov.40_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.c.cov.40.cl.prop.women15.25.M.40.50)


mcar_mar_c_cov.40_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.c.cov.40.cl.prop.men40.50.F.25.40)


mcar_mar_c_cov.40_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.40.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.c.cov.40.cl.prop.women25.40.M.40.50)



# 45 - c 


mcar_mar_c_cov.45_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.c.cov.45.cl.prop.men15.25.F.15.25)


mcar_mar_c_cov.45_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.c.cov.45.cl.prop.women15.25.M.15.25)


mcar_mar_c_cov.45_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.c.cov.45.cl.prop.men25.40.F.15.25)



mcar_mar_c_cov.45_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.c.cov.45.cl.prop.women15.25.M.25.40)

mcar_mar_c_cov.45_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.c.cov.45.cl.prop.men25.40.F.25.40)


mcar_mar_c_cov.45_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.c.cov.45.cl.prop.women25.40.M.25.40)


mcar_mar_c_cov.45_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.c.cov.45.cl.prop.men40.50.F.15.25)


mcar_mar_c_cov.45_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.c.cov.45.cl.prop.women15.25.M.40.50)


mcar_mar_c_cov.45_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.c.cov.45.cl.prop.men40.50.F.25.40)


mcar_mar_c_cov.45_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.45.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.c.cov.45.cl.prop.women25.40.M.40.50)



# 50 - c


mcar_mar_c_cov.50_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.c.cov.50.cl.prop.men15.25.F.15.25)


mcar_mar_c_cov.50_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.c.cov.50.cl.prop.women15.25.M.15.25)


mcar_mar_c_cov.50_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.c.cov.50.cl.prop.men25.40.F.15.25)



mcar_mar_c_cov.50_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.c.cov.50.cl.prop.women15.25.M.25.40)

mcar_mar_c_cov.50_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.c.cov.50.cl.prop.men25.40.F.25.40)


mcar_mar_c_cov.50_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.c.cov.50.cl.prop.women25.40.M.25.40)


mcar_mar_c_cov.50_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.c.cov.50.cl.prop.men40.50.F.15.25)


mcar_mar_c_cov.50_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.c.cov.50.cl.prop.women15.25.M.40.50)


mcar_mar_c_cov.50_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.c.cov.50.cl.prop.men40.50.F.25.40)


mcar_mar_c_cov.50_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.50.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.c.cov.50.cl.prop.women25.40.M.40.50)



# 55 - c


mcar_mar_c_cov.55_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.c.cov.55.cl.prop.men15.25.F.15.25)


mcar_mar_c_cov.55_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.c.cov.55.cl.prop.women15.25.M.15.25)


mcar_mar_c_cov.55_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.c.cov.55.cl.prop.men25.40.F.15.25)



mcar_mar_c_cov.55_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.c.cov.55.cl.prop.women15.25.M.25.40)

mcar_mar_c_cov.55_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.c.cov.55.cl.prop.men25.40.F.25.40)


mcar_mar_c_cov.55_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.c.cov.55.cl.prop.women25.40.M.25.40)


mcar_mar_c_cov.55_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.c.cov.55.cl.prop.men40.50.F.15.25)


mcar_mar_c_cov.55_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.c.cov.55.cl.prop.women15.25.M.40.50)


mcar_mar_c_cov.55_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.c.cov.55.cl.prop.men40.50.F.25.40)


mcar_mar_c_cov.55_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.55.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.c.cov.55.cl.prop.women25.40.M.40.50)



# 60 - c 


mcar_mar_c_cov.60_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.c.cov.60.cl.prop.men15.25.F.15.25)


mcar_mar_c_cov.60_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.c.cov.60.cl.prop.women15.25.M.15.25)


mcar_mar_c_cov.60_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.c.cov.60.cl.prop.men25.40.F.15.25)



mcar_mar_c_cov.60_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.c.cov.60.cl.prop.women15.25.M.25.40)

mcar_mar_c_cov.60_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.c.cov.60.cl.prop.men25.40.F.25.40)


mcar_mar_c_cov.60_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.c.cov.60.cl.prop.women25.40.M.25.40)


mcar_mar_c_cov.60_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.c.cov.60.cl.prop.men40.50.F.15.25)

mcar_mar_c_cov.60_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.c.cov.60.cl.prop.women15.25.M.40.50)


mcar_mar_c_cov.60_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.c.cov.60.cl.prop.men40.50.F.25.40)


mcar_mar_c_cov.60_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.60.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.c.cov.60.cl.prop.women25.40.M.40.50)



# 65 - c


mcar_mar_c_cov.65_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.c.cov.65.cl.prop.men15.25.F.15.25)


mcar_mar_c_cov.65_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.c.cov.65.cl.prop.women15.25.M.15.25)


mcar_mar_c_cov.65_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.c.cov.65.cl.prop.men25.40.F.15.25)



mcar_mar_c_cov.65_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.c.cov.65.cl.prop.women15.25.M.25.40)

mcar_mar_c_cov.65_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.c.cov.65.cl.prop.men25.40.F.25.40)


mcar_mar_c_cov.65_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.c.cov.65.cl.prop.women25.40.M.25.40)


mcar_mar_c_cov.65_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.c.cov.65.cl.prop.men40.50.F.15.25)

mcar_mar_c_cov.65_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.c.cov.65.cl.prop.women15.25.M.40.50)


mcar_mar_c_cov.65_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.c.cov.65.cl.prop.men40.50.F.25.40)


mcar_mar_c_cov.65_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.65.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.c.cov.65.cl.prop.women25.40.M.40.50)



# 65 - c


mcar_mar_c_cov.70_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.c.cov.70.cl.prop.men15.25.F.15.25)


mcar_mar_c_cov.70_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.c.cov.70.cl.prop.women15.25.M.15.25)


mcar_mar_c_cov.70_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.c.cov.70.cl.prop.men25.40.F.15.25)



mcar_mar_c_cov.70_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.c.cov.70.cl.prop.women15.25.M.25.40)

mcar_mar_c_cov.70_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.c.cov.70.cl.prop.men25.40.F.25.40)


mcar_mar_c_cov.70_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.c.cov.70.cl.prop.women25.40.M.25.40)


mcar_mar_c_cov.70_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.c.cov.70.cl.prop.men40.50.F.15.25)

mcar_mar_c_cov.70_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.c.cov.70.cl.prop.women15.25.M.40.50)

mcar_mar_c_cov.70_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.c.cov.70.cl.prop.men40.50.F.25.40)


mcar_mar_c_cov.70_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.70.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.c.cov.70.cl.prop.women25.40.M.40.50)



# 75 - c


mcar_mar_c_cov.75_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.c.cov.75.cl.prop.men15.25.F.15.25)


mcar_mar_c_cov.75_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.c.cov.75.cl.prop.women15.25.M.15.25)


mcar_mar_c_cov.75_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.c.cov.75.cl.prop.men25.40.F.15.25)



mcar_mar_c_cov.75_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.c.cov.75.cl.prop.women15.25.M.25.40)

mcar_mar_c_cov.75_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.c.cov.75.cl.prop.men25.40.F.25.40)


mcar_mar_c_cov.75_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.c.cov.75.cl.prop.women25.40.M.25.40)


mcar_mar_c_cov.75_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.c.cov.75.cl.prop.men40.50.F.15.25)


mcar_mar_c_cov.75_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.c.cov.75.cl.prop.women15.25.M.40.50)


mcar_mar_c_cov.75_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.c.cov.75.cl.prop.men40.50.F.25.40)


mcar_mar_c_cov.75_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.75.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.c.cov.75.cl.prop.women25.40.M.40.50)


# 80 - c


mcar_mar_c_cov.80_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.c.cov.80.cl.prop.men15.25.F.15.25)


mcar_mar_c_cov.80_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.c.cov.80.cl.prop.women15.25.M.15.25)


mcar_mar_c_cov.80_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.c.cov.80.cl.prop.men25.40.F.15.25)



mcar_mar_c_cov.80_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.c.cov.80.cl.prop.women15.25.M.25.40)

mcar_mar_c_cov.80_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.c.cov.80.cl.prop.men25.40.F.25.40)


mcar_mar_c_cov.80_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.c.cov.80.cl.prop.women25.40.M.25.40)


mcar_mar_c_cov.80_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.c.cov.80.cl.prop.men40.50.F.15.25)


mcar_mar_c_cov.80_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.c.cov.80.cl.prop.women15.25.M.40.50)


mcar_mar_c_cov.80_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.c.cov.80.cl.prop.men40.50.F.25.40)


mcar_mar_c_cov.80_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.80.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.c.cov.80.cl.prop.women25.40.M.40.50)



# 85 - c


mcar_mar_c_cov.85_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.c.cov.85.cl.prop.men15.25.F.15.25)


mcar_mar_c_cov.85_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.c.cov.85.cl.prop.women15.25.M.15.25)


mcar_mar_c_cov.85_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.c.cov.85.cl.prop.men25.40.F.15.25)



mcar_mar_c_cov.85_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.c.cov.85.cl.prop.women15.25.M.25.40)

mcar_mar_c_cov.85_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.c.cov.85.cl.prop.men25.40.F.25.40)


mcar_mar_c_cov.85_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.c.cov.85.cl.prop.women25.40.M.25.40)


mcar_mar_c_cov.85_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.c.cov.85.cl.prop.men40.50.F.15.25)


mcar_mar_c_cov.85_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.c.cov.85.cl.prop.women15.25.M.40.50)


mcar_mar_c_cov.85_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.c.cov.85.cl.prop.men40.50.F.25.40)


mcar_mar_c_cov.85_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.85.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.c.cov.85.cl.prop.women25.40.M.40.50)



# 90 - c


mcar_mar_c_cov.90_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.c.cov.90.cl.prop.men15.25.F.15.25)


mcar_mar_c_cov.90_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.c.cov.90.cl.prop.women15.25.M.15.25)


mcar_mar_c_cov.90_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.c.cov.90.cl.prop.men25.40.F.15.25)



mcar_mar_c_cov.90_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.c.cov.90.cl.prop.women15.25.M.25.40)

mcar_mar_c_cov.90_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.c.cov.90.cl.prop.men25.40.F.25.40)


mcar_mar_c_cov.90_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.c.cov.90.cl.prop.women25.40.M.25.40)


mcar_mar_c_cov.90_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.c.cov.90.cl.prop.men40.50.F.15.25)


mcar_mar_c_cov.90_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.c.cov.90.cl.prop.women15.25.M.40.50)


mcar_mar_c_cov.90_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.c.cov.90.cl.prop.men40.50.F.25.40)


mcar_mar_c_cov.90_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.90.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.c.cov.90.cl.prop.women25.40.M.40.50)



# 95 - c


mcar_mar_c_cov.95_prop.men15.25.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.men15.25.F.15.25, 
                                                           mar = vector.MAR.c.cov.95.cl.prop.men15.25.F.15.25)


mcar_mar_c_cov.95_prop.women15.25.M.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.women15.25.M.15.25, 
                                                             mar = vector.MAR.c.cov.95.cl.prop.women15.25.M.15.25)


mcar_mar_c_cov.95_prop.men25.40.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.men25.40.F.15.25, 
                                                           mar = vector.MAR.c.cov.95.cl.prop.men25.40.F.15.25)



mcar_mar_c_cov.95_prop.women15.25.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.women15.25.M.25.40, 
                                                             mar = vector.MAR.c.cov.95.cl.prop.women15.25.M.25.40)

mcar_mar_c_cov.95_prop.men25.40.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.men25.40.F.25.40, 
                                                           mar = vector.MAR.c.cov.95.cl.prop.men25.40.F.25.40)


mcar_mar_c_cov.95_prop.women25.40.M.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.women25.40.M.25.40, 
                                                             mar = vector.MAR.c.cov.95.cl.prop.women25.40.M.25.40)


mcar_mar_c_cov.95_prop.men40.50.F.15.25 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.men40.50.F.15.25, 
                                                           mar = vector.MAR.c.cov.95.cl.prop.men40.50.F.15.25)


mcar_mar_c_cov.95_prop.women15.25.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.women15.25.M.40.50, 
                                                             mar = vector.MAR.c.cov.95.cl.prop.women15.25.M.40.50)


mcar_mar_c_cov.95_prop.men40.50.F.25.40 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.men40.50.F.25.40, 
                                                           mar = vector.MAR.c.cov.95.cl.prop.men40.50.F.25.40)


mcar_mar_c_cov.95_prop.women25.40.M.40.50 <- wilcox.test.A.B(mcar = vector.MCAR.cov.95.cl.prop.women25.40.M.40.50, 
                                                             mar = vector.MAR.c.cov.95.cl.prop.women25.40.M.40.50)




# Age difference MCAR and MAR (a) -------


# mean

mcar_mar_a_cov.35_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.35.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.a.cov.35.AD.men.cl.15.25)


mcar_mar_a_cov.35_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.35.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.a.cov.35.AD.women.cl.15.25)


mcar_mar_a_cov.35_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.35.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.a.cov.35.AD.men.cl.25.40)


mcar_mar_a_cov.35_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.35.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.a.cov.35.AD.women.cl.25.40)


mcar_mar_a_cov.35_mean.men.cl.40.50<- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.35.AD.men.cl.40.50, 
                                                      mar = vector.mean.MAR.a.cov.35.AD.men.cl.40.50)


mcar_mar_a_cov.35_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.35.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.a.cov.35.AD.women.cl.40.50)


# median


mcar_mar_a_cov.35_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.35.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.a.cov.35.AD.men.cl.15.25)


mcar_mar_a_cov.35_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.35.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.a.cov.35.AD.women.cl.15.25)


mcar_mar_a_cov.35_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.35.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.a.cov.35.AD.men.cl.25.40)


mcar_mar_a_cov.35_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.35.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.a.cov.35.AD.women.cl.25.40)


mcar_mar_a_cov.35_med.men.cl.40.50<- wilcox.test.A.B(mcar = vector.med.MCAR.cov.35.AD.men.cl.40.50, 
                                                     mar = vector.med.MAR.a.cov.35.AD.men.cl.40.50)


mcar_mar_a_cov.35_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.35.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.a.cov.35.AD.women.cl.40.50)


# standard deviation

mcar_mar_a_cov.35_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.35.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.a.cov.35.AD.men.cl.15.25)


mcar_mar_a_cov.35_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.35.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.a.cov.35.AD.women.cl.15.25)


mcar_mar_a_cov.35_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.35.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.a.cov.35.AD.men.cl.25.40)


mcar_mar_a_cov.35_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.35.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.a.cov.35.AD.women.cl.25.40)


mcar_mar_a_cov.35_sd.men.cl.40.50<- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.35.AD.men.cl.40.50, 
                                                    mar = vector.sd.MAR.a.cov.35.AD.men.cl.40.50)


mcar_mar_a_cov.35_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.35.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.a.cov.35.AD.women.cl.40.50)


# 40 - a





# mean

mcar_mar_a_cov.40_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.40.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.a.cov.40.AD.men.cl.15.25)


mcar_mar_a_cov.40_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.40.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.a.cov.40.AD.women.cl.15.25)


mcar_mar_a_cov.40_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.40.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.a.cov.40.AD.men.cl.25.40)


mcar_mar_a_cov.40_mean.women.cl.25.40<- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.40.AD.women.cl.25.40, 
                                                        mar = vector.mean.MAR.a.cov.40.AD.women.cl.25.40)


mcar_mar_a_cov.40_mean.men.cl.40.50<- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.40.AD.men.cl.40.50, 
                                                      mar = vector.mean.MAR.a.cov.40.AD.men.cl.40.50)


mcar_mar_a_cov.40_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.40.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.a.cov.40.AD.women.cl.40.50)


# median


mcar_mar_a_cov.40_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.40.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.a.cov.40.AD.men.cl.15.25)


mcar_mar_a_cov.40_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.40.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.a.cov.40.AD.women.cl.15.25)


mcar_mar_a_cov.40_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.40.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.a.cov.40.AD.men.cl.25.40)


mcar_mar_a_cov.40_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.40.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.a.cov.40.AD.women.cl.25.40)


mcar_mar_a_cov.40_med.men.cl.40.50<- wilcox.test.A.B(mcar = vector.med.MCAR.cov.40.AD.men.cl.40.50, 
                                                     mar = vector.med.MAR.a.cov.40.AD.men.cl.40.50)


mcar_mar_a_cov.40_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.40.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.a.cov.40.AD.women.cl.40.50)


# standard deviation

mcar_mar_a_cov.40_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.40.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.a.cov.40.AD.men.cl.15.25)


mcar_mar_a_cov.40_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.40.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.a.cov.40.AD.women.cl.15.25)


mcar_mar_a_cov.40_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.40.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.a.cov.40.AD.men.cl.25.40)


mcar_mar_a_cov.40_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.40.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.a.cov.40.AD.women.cl.25.40)


mcar_mar_a_cov.40_sd.men.cl.40.50<- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.40.AD.men.cl.40.50, 
                                                    mar = vector.sd.MAR.a.cov.40.AD.men.cl.40.50)


mcar_mar_a_cov.40_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.40.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.a.cov.40.AD.women.cl.40.50)



# 45 - a





# mean

mcar_mar_a_cov.45_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.45.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.a.cov.45.AD.men.cl.15.25)


mcar_mar_a_cov.45_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.45.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.a.cov.45.AD.women.cl.15.25)


mcar_mar_a_cov.45_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.45.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.a.cov.45.AD.men.cl.25.40)


mcar_mar_a_cov.45_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.45.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.a.cov.45.AD.women.cl.25.40)


mcar_mar_a_cov.45_mean.men.cl.40.50<- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.45.AD.men.cl.40.50, 
                                                      mar = vector.mean.MAR.a.cov.45.AD.men.cl.40.50)


mcar_mar_a_cov.45_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.45.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.a.cov.45.AD.women.cl.40.50)


# median


mcar_mar_a_cov.45_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.45.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.a.cov.45.AD.men.cl.15.25)


mcar_mar_a_cov.45_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.45.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.a.cov.45.AD.women.cl.15.25)


mcar_mar_a_cov.45_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.45.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.a.cov.45.AD.men.cl.25.40)


mcar_mar_a_cov.45_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.45.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.a.cov.45.AD.women.cl.25.40)


mcar_mar_a_cov.45_med.men.cl.40.50<- wilcox.test.A.B(mcar = vector.med.MCAR.cov.45.AD.men.cl.40.50, 
                                                     mar = vector.med.MAR.a.cov.45.AD.men.cl.40.50)


mcar_mar_a_cov.45_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.45.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.a.cov.45.AD.women.cl.40.50)


# standard deviation

mcar_mar_a_cov.45_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.45.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.a.cov.45.AD.men.cl.15.25)


mcar_mar_a_cov.45_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.45.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.a.cov.45.AD.women.cl.15.25)


mcar_mar_a_cov.45_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.45.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.a.cov.45.AD.men.cl.25.40)


mcar_mar_a_cov.45_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.45.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.a.cov.45.AD.women.cl.25.40)


mcar_mar_a_cov.45_sd.men.cl.40.50<- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.45.AD.men.cl.40.50, 
                                                    mar = vector.sd.MAR.a.cov.45.AD.men.cl.40.50)


mcar_mar_a_cov.45_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.45.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.a.cov.45.AD.women.cl.40.50)



# 50 - a




# mean

mcar_mar_a_cov.50_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.50.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.a.cov.50.AD.men.cl.15.25)


mcar_mar_a_cov.50_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.50.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.a.cov.50.AD.women.cl.15.25)


mcar_mar_a_cov.50_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.50.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.a.cov.50.AD.men.cl.25.40)


mcar_mar_a_cov.50_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.50.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.a.cov.50.AD.women.cl.25.40)


mcar_mar_a_cov.50_mean.men.cl.40.50<- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.50.AD.men.cl.40.50, 
                                                      mar = vector.mean.MAR.a.cov.50.AD.men.cl.40.50)


mcar_mar_a_cov.50_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.50.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.a.cov.50.AD.women.cl.40.50)


# median


mcar_mar_a_cov.50_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.50.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.a.cov.50.AD.men.cl.15.25)


mcar_mar_a_cov.50_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.50.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.a.cov.50.AD.women.cl.15.25)


mcar_mar_a_cov.50_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.50.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.a.cov.50.AD.men.cl.25.40)


mcar_mar_a_cov.50_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.50.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.a.cov.50.AD.women.cl.25.40)


mcar_mar_a_cov.50_med.men.cl.40.50<- wilcox.test.A.B(mcar = vector.med.MCAR.cov.50.AD.men.cl.40.50, 
                                                     mar = vector.med.MAR.a.cov.50.AD.men.cl.40.50)


mcar_mar_a_cov.50_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.50.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.a.cov.50.AD.women.cl.40.50)


# standard deviation

mcar_mar_a_cov.50_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.50.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.a.cov.50.AD.men.cl.15.25)


mcar_mar_a_cov.50_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.50.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.a.cov.50.AD.women.cl.15.25)


mcar_mar_a_cov.50_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.50.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.a.cov.50.AD.men.cl.25.40)


mcar_mar_a_cov.50_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.50.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.a.cov.50.AD.women.cl.25.40)


mcar_mar_a_cov.50_sd.men.cl.40.50<- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.50.AD.men.cl.40.50, 
                                                    mar = vector.sd.MAR.a.cov.50.AD.men.cl.40.50)


mcar_mar_a_cov.50_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.50.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.a.cov.50.AD.women.cl.40.50)



# 55 - a




# mean

mcar_mar_a_cov.55_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.55.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.a.cov.55.AD.men.cl.15.25)


mcar_mar_a_cov.55_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.55.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.a.cov.55.AD.women.cl.15.25)


mcar_mar_a_cov.55_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.55.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.a.cov.55.AD.men.cl.25.40)


mcar_mar_a_cov.55_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.55.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.a.cov.55.AD.women.cl.25.40)


mcar_mar_a_cov.55_mean.men.cl.40.50<- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.55.AD.men.cl.40.50, 
                                                      mar = vector.mean.MAR.a.cov.55.AD.men.cl.40.50)


mcar_mar_a_cov.55_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.55.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.a.cov.55.AD.women.cl.40.50)


# median


mcar_mar_a_cov.55_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.55.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.a.cov.55.AD.men.cl.15.25)


mcar_mar_a_cov.55_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.55.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.a.cov.55.AD.women.cl.15.25)


mcar_mar_a_cov.55_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.55.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.a.cov.55.AD.men.cl.25.40)


mcar_mar_a_cov.55_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.55.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.a.cov.55.AD.women.cl.25.40)


mcar_mar_a_cov.55_med.men.cl.40.50<- wilcox.test.A.B(mcar = vector.med.MCAR.cov.55.AD.men.cl.40.50, 
                                                     mar = vector.med.MAR.a.cov.55.AD.men.cl.40.50)


mcar_mar_a_cov.55_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.55.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.a.cov.55.AD.women.cl.40.50)


# standard deviation

mcar_mar_a_cov.55_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.55.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.a.cov.55.AD.men.cl.15.25)


mcar_mar_a_cov.55_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.55.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.a.cov.55.AD.women.cl.15.25)


mcar_mar_a_cov.55_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.55.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.a.cov.55.AD.men.cl.25.40)


mcar_mar_a_cov.55_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.55.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.a.cov.55.AD.women.cl.25.40)


mcar_mar_a_cov.55_sd.men.cl.40.50<- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.55.AD.men.cl.40.50, 
                                                    mar = vector.sd.MAR.a.cov.55.AD.men.cl.40.50)


mcar_mar_a_cov.55_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.55.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.a.cov.55.AD.women.cl.40.50)



# 60 - a





# mean

mcar_mar_a_cov.60_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.60.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.a.cov.60.AD.men.cl.15.25)


mcar_mar_a_cov.60_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.60.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.a.cov.60.AD.women.cl.15.25)


mcar_mar_a_cov.60_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.60.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.a.cov.60.AD.men.cl.25.40)


mcar_mar_a_cov.60_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.60.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.a.cov.60.AD.women.cl.25.40)


mcar_mar_a_cov.60_mean.men.cl.40.50<- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.60.AD.men.cl.40.50, 
                                                      mar = vector.mean.MAR.a.cov.60.AD.men.cl.40.50)


mcar_mar_a_cov.60_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.60.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.a.cov.60.AD.women.cl.40.50)


# median


mcar_mar_a_cov.60_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.60.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.a.cov.60.AD.men.cl.15.25)


mcar_mar_a_cov.60_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.60.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.a.cov.60.AD.women.cl.15.25)


mcar_mar_a_cov.60_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.60.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.a.cov.60.AD.men.cl.25.40)


mcar_mar_a_cov.60_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.60.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.a.cov.60.AD.women.cl.25.40)


mcar_mar_a_cov.60_med.men.cl.40.50<- wilcox.test.A.B(mcar = vector.med.MCAR.cov.60.AD.men.cl.40.50, 
                                                     mar = vector.med.MAR.a.cov.60.AD.men.cl.40.50)


mcar_mar_a_cov.60_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.60.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.a.cov.60.AD.women.cl.40.50)


# standard deviation

mcar_mar_a_cov.60_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.60.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.a.cov.60.AD.men.cl.15.25)


mcar_mar_a_cov.60_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.60.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.a.cov.60.AD.women.cl.15.25)


mcar_mar_a_cov.60_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.60.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.a.cov.60.AD.men.cl.25.40)


mcar_mar_a_cov.60_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.60.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.a.cov.60.AD.women.cl.25.40)


mcar_mar_a_cov.60_sd.men.cl.40.50<- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.60.AD.men.cl.40.50, 
                                                    mar = vector.sd.MAR.a.cov.60.AD.men.cl.40.50)


mcar_mar_a_cov.60_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.60.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.a.cov.60.AD.women.cl.40.50)



# 65 - a




# mean

mcar_mar_a_cov.65_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.65.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.a.cov.65.AD.men.cl.15.25)


mcar_mar_a_cov.65_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.65.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.a.cov.65.AD.women.cl.15.25)


mcar_mar_a_cov.65_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.65.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.a.cov.65.AD.men.cl.25.40)


mcar_mar_a_cov.65_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.65.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.a.cov.65.AD.women.cl.25.40)


mcar_mar_a_cov.65_mean.men.cl.40.50<- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.65.AD.men.cl.40.50, 
                                                      mar = vector.mean.MAR.a.cov.65.AD.men.cl.40.50)


mcar_mar_a_cov.65_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.65.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.a.cov.65.AD.women.cl.40.50)


# median


mcar_mar_a_cov.65_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.65.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.a.cov.65.AD.men.cl.15.25)


mcar_mar_a_cov.65_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.65.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.a.cov.65.AD.women.cl.15.25)


mcar_mar_a_cov.65_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.65.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.a.cov.65.AD.men.cl.25.40)


mcar_mar_a_cov.65_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.65.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.a.cov.65.AD.women.cl.25.40)


mcar_mar_a_cov.65_med.men.cl.40.50<- wilcox.test.A.B(mcar = vector.med.MCAR.cov.65.AD.men.cl.40.50, 
                                                     mar = vector.med.MAR.a.cov.65.AD.men.cl.40.50)


mcar_mar_a_cov.65_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.65.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.a.cov.65.AD.women.cl.40.50)


# standard deviation

mcar_mar_a_cov.65_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.65.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.a.cov.65.AD.men.cl.15.25)


mcar_mar_a_cov.65_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.65.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.a.cov.65.AD.women.cl.15.25)


mcar_mar_a_cov.65_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.65.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.a.cov.65.AD.men.cl.25.40)


mcar_mar_a_cov.65_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.65.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.a.cov.65.AD.women.cl.25.40)


mcar_mar_a_cov.65_sd.men.cl.40.50<- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.65.AD.men.cl.40.50, 
                                                    mar = vector.sd.MAR.a.cov.65.AD.men.cl.40.50)


mcar_mar_a_cov.65_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.65.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.a.cov.65.AD.women.cl.40.50)



# 70 - a





# mean

mcar_mar_a_cov.70_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.70.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.a.cov.70.AD.men.cl.15.25)


mcar_mar_a_cov.70_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.70.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.a.cov.70.AD.women.cl.15.25)


mcar_mar_a_cov.70_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.70.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.a.cov.70.AD.men.cl.25.40)


mcar_mar_a_cov.70_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.70.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.a.cov.70.AD.women.cl.25.40)


mcar_mar_a_cov.70_mean.men.cl.40.50<- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.70.AD.men.cl.40.50, 
                                                      mar = vector.mean.MAR.a.cov.70.AD.men.cl.40.50)


mcar_mar_a_cov.70_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.70.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.a.cov.70.AD.women.cl.40.50)


# median


mcar_mar_a_cov.70_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.70.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.a.cov.70.AD.men.cl.15.25)


mcar_mar_a_cov.70_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.70.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.a.cov.70.AD.women.cl.15.25)


mcar_mar_a_cov.70_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.70.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.a.cov.70.AD.men.cl.25.40)


mcar_mar_a_cov.70_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.70.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.a.cov.70.AD.women.cl.25.40)


mcar_mar_a_cov.70_med.men.cl.40.50<- wilcox.test.A.B(mcar = vector.med.MCAR.cov.70.AD.men.cl.40.50, 
                                                     mar = vector.med.MAR.a.cov.70.AD.men.cl.40.50)


mcar_mar_a_cov.70_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.70.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.a.cov.70.AD.women.cl.40.50)


# standard deviation

mcar_mar_a_cov.70_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.70.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.a.cov.70.AD.men.cl.15.25)


mcar_mar_a_cov.70_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.70.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.a.cov.70.AD.women.cl.15.25)


mcar_mar_a_cov.70_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.70.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.a.cov.70.AD.men.cl.25.40)


mcar_mar_a_cov.70_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.70.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.a.cov.70.AD.women.cl.25.40)


mcar_mar_a_cov.70_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.70.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.a.cov.70.AD.men.cl.40.50)


mcar_mar_a_cov.70_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.70.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.a.cov.70.AD.women.cl.40.50)



# 75 - a




# mean

mcar_mar_a_cov.75_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.75.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.a.cov.75.AD.men.cl.15.25)


mcar_mar_a_cov.75_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.75.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.a.cov.75.AD.women.cl.15.25)


mcar_mar_a_cov.75_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.75.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.a.cov.75.AD.men.cl.25.40)


mcar_mar_a_cov.75_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.75.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.a.cov.75.AD.women.cl.25.40)


mcar_mar_a_cov.75_mean.men.cl.40.50<- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.75.AD.men.cl.40.50, 
                                                      mar = vector.mean.MAR.a.cov.75.AD.men.cl.40.50)


mcar_mar_a_cov.75_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.75.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.a.cov.75.AD.women.cl.40.50)


# median


mcar_mar_a_cov.75_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.75.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.a.cov.75.AD.men.cl.15.25)


mcar_mar_a_cov.75_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.75.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.a.cov.75.AD.women.cl.15.25)


mcar_mar_a_cov.75_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.75.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.a.cov.75.AD.men.cl.25.40)


mcar_mar_a_cov.75_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.75.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.a.cov.75.AD.women.cl.25.40)


mcar_mar_a_cov.75_med.men.cl.40.50<- wilcox.test.A.B(mcar = vector.med.MCAR.cov.75.AD.men.cl.40.50, 
                                                     mar = vector.med.MAR.a.cov.75.AD.men.cl.40.50)


mcar_mar_a_cov.75_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.75.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.a.cov.75.AD.women.cl.40.50)


# standard deviation

mcar_mar_a_cov.75_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.75.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.a.cov.75.AD.men.cl.15.25)


mcar_mar_a_cov.75_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.75.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.a.cov.75.AD.women.cl.15.25)


mcar_mar_a_cov.75_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.75.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.a.cov.75.AD.men.cl.25.40)


mcar_mar_a_cov.75_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.75.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.a.cov.75.AD.women.cl.25.40)


mcar_mar_a_cov.75_sd.men.cl.40.50<- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.75.AD.men.cl.40.50, 
                                                    mar = vector.sd.MAR.a.cov.75.AD.men.cl.40.50)


mcar_mar_a_cov.75_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.75.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.a.cov.75.AD.women.cl.40.50)



# 80 - a




# mean

mcar_mar_a_cov.80_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.80.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.a.cov.80.AD.men.cl.15.25)


mcar_mar_a_cov.80_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.80.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.a.cov.80.AD.women.cl.15.25)


mcar_mar_a_cov.80_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.80.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.a.cov.80.AD.men.cl.25.40)


mcar_mar_a_cov.80_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.80.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.a.cov.80.AD.women.cl.25.40)


mcar_mar_a_cov.80_mean.men.cl.40.50<- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.80.AD.men.cl.40.50, 
                                                      mar = vector.mean.MAR.a.cov.80.AD.men.cl.40.50)


mcar_mar_a_cov.80_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.80.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.a.cov.80.AD.women.cl.40.50)


# median


mcar_mar_a_cov.80_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.80.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.a.cov.80.AD.men.cl.15.25)


mcar_mar_a_cov.80_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.80.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.a.cov.80.AD.women.cl.15.25)


mcar_mar_a_cov.80_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.80.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.a.cov.80.AD.men.cl.25.40)


mcar_mar_a_cov.80_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.80.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.a.cov.80.AD.women.cl.25.40)


mcar_mar_a_cov.80_med.men.cl.40.50<- wilcox.test.A.B(mcar = vector.med.MCAR.cov.80.AD.men.cl.40.50, 
                                                     mar = vector.med.MAR.a.cov.80.AD.men.cl.40.50)


mcar_mar_a_cov.80_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.80.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.a.cov.80.AD.women.cl.40.50)


# standard deviation

mcar_mar_a_cov.80_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.80.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.a.cov.80.AD.men.cl.15.25)


mcar_mar_a_cov.80_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.80.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.a.cov.80.AD.women.cl.15.25)


mcar_mar_a_cov.80_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.80.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.a.cov.80.AD.men.cl.25.40)


mcar_mar_a_cov.80_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.80.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.a.cov.80.AD.women.cl.25.40)


mcar_mar_a_cov.80_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.80.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.a.cov.80.AD.men.cl.40.50)


mcar_mar_a_cov.80_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.80.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.a.cov.80.AD.women.cl.40.50)



# 85 - a




# mean

mcar_mar_a_cov.85_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.85.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.a.cov.85.AD.men.cl.15.25)


mcar_mar_a_cov.85_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.85.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.a.cov.85.AD.women.cl.15.25)


mcar_mar_a_cov.85_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.85.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.a.cov.85.AD.men.cl.25.40)


mcar_mar_a_cov.85_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.85.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.a.cov.85.AD.women.cl.25.40)


mcar_mar_a_cov.85_mean.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.85.AD.men.cl.40.50, 
                                                       mar = vector.mean.MAR.a.cov.85.AD.men.cl.40.50)


mcar_mar_a_cov.85_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.85.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.a.cov.85.AD.women.cl.40.50)


# median


mcar_mar_a_cov.85_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.85.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.a.cov.85.AD.men.cl.15.25)


mcar_mar_a_cov.85_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.85.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.a.cov.85.AD.women.cl.15.25)


mcar_mar_a_cov.85_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.85.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.a.cov.85.AD.men.cl.25.40)


mcar_mar_a_cov.85_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.85.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.a.cov.85.AD.women.cl.25.40)


mcar_mar_a_cov.85_med.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.85.AD.men.cl.40.50, 
                                                      mar = vector.med.MAR.a.cov.85.AD.men.cl.40.50)


mcar_mar_a_cov.85_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.85.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.a.cov.85.AD.women.cl.40.50)


# standard deviation

mcar_mar_a_cov.85_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.85.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.a.cov.85.AD.men.cl.15.25)


mcar_mar_a_cov.85_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.85.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.a.cov.85.AD.women.cl.15.25)


mcar_mar_a_cov.85_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.85.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.a.cov.85.AD.men.cl.25.40)


mcar_mar_a_cov.85_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.85.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.a.cov.85.AD.women.cl.25.40)


mcar_mar_a_cov.85_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.85.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.a.cov.85.AD.men.cl.40.50)


mcar_mar_a_cov.85_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.85.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.a.cov.85.AD.women.cl.40.50)



# 90 - a




# mean

mcar_mar_a_cov.90_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.90.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.a.cov.90.AD.men.cl.15.25)


mcar_mar_a_cov.90_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.90.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.a.cov.90.AD.women.cl.15.25)


mcar_mar_a_cov.90_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.90.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.a.cov.90.AD.men.cl.25.40)


mcar_mar_a_cov.90_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.90.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.a.cov.90.AD.women.cl.25.40)


mcar_mar_a_cov.90_mean.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.90.AD.men.cl.40.50, 
                                                       mar = vector.mean.MAR.a.cov.90.AD.men.cl.40.50)


mcar_mar_a_cov.90_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.90.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.a.cov.90.AD.women.cl.40.50)


# median


mcar_mar_a_cov.90_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.90.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.a.cov.90.AD.men.cl.15.25)


mcar_mar_a_cov.90_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.90.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.a.cov.90.AD.women.cl.15.25)


mcar_mar_a_cov.90_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.90.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.a.cov.90.AD.men.cl.25.40)


mcar_mar_a_cov.90_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.90.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.a.cov.90.AD.women.cl.25.40)


mcar_mar_a_cov.90_med.men.cl.40.50<- wilcox.test.A.B(mcar = vector.med.MCAR.cov.90.AD.men.cl.40.50, 
                                                     mar = vector.med.MAR.a.cov.90.AD.men.cl.40.50)


mcar_mar_a_cov.90_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.90.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.a.cov.90.AD.women.cl.40.50)


# standard deviation

mcar_mar_a_cov.90_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.90.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.a.cov.90.AD.men.cl.15.25)


mcar_mar_a_cov.90_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.90.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.a.cov.90.AD.women.cl.15.25)


mcar_mar_a_cov.90_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.90.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.a.cov.90.AD.men.cl.25.40)


mcar_mar_a_cov.90_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.90.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.a.cov.90.AD.women.cl.25.40)


mcar_mar_a_cov.90_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.90.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.a.cov.90.AD.men.cl.40.50)


mcar_mar_a_cov.90_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.90.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.a.cov.90.AD.women.cl.40.50)



# 95 - a




# mean

mcar_mar_a_cov.95_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.95.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.a.cov.95.AD.men.cl.15.25)


mcar_mar_a_cov.95_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.95.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.a.cov.95.AD.women.cl.15.25)


mcar_mar_a_cov.95_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.95.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.a.cov.95.AD.men.cl.25.40)


mcar_mar_a_cov.95_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.95.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.a.cov.95.AD.women.cl.25.40)


mcar_mar_a_cov.95_mean.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.95.AD.men.cl.40.50, 
                                                       mar = vector.mean.MAR.a.cov.95.AD.men.cl.40.50)


mcar_mar_a_cov.95_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.95.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.a.cov.95.AD.women.cl.40.50)


# median


mcar_mar_a_cov.95_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.95.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.a.cov.95.AD.men.cl.15.25)


mcar_mar_a_cov.95_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.95.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.a.cov.95.AD.women.cl.15.25)


mcar_mar_a_cov.95_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.95.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.a.cov.95.AD.men.cl.25.40)


mcar_mar_a_cov.95_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.95.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.a.cov.95.AD.women.cl.25.40)


mcar_mar_a_cov.95_med.men.cl.40.50<- wilcox.test.A.B(mcar = vector.med.MCAR.cov.95.AD.men.cl.40.50, 
                                                     mar = vector.med.MAR.a.cov.95.AD.men.cl.40.50)


mcar_mar_a_cov.95_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.95.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.a.cov.95.AD.women.cl.40.50)


# standard deviation

mcar_mar_a_cov.95_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.95.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.a.cov.95.AD.men.cl.15.25)


mcar_mar_a_cov.95_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.95.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.a.cov.95.AD.women.cl.15.25)


mcar_mar_a_cov.95_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.95.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.a.cov.95.AD.men.cl.25.40)


mcar_mar_a_cov.95_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.95.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.a.cov.95.AD.women.cl.25.40)


mcar_mar_a_cov.95_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.95.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.a.cov.95.AD.men.cl.40.50)


mcar_mar_a_cov.95_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.95.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.a.cov.95.AD.women.cl.40.50)




# Age difference MCAR and MAR (b) -------


# mean

mcar_mar_b_cov.35_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.35.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.b.cov.35.AD.men.cl.15.25)


mcar_mar_b_cov.35_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.35.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.b.cov.35.AD.women.cl.15.25)


mcar_mar_b_cov.35_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.35.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.b.cov.35.AD.men.cl.25.40)


mcar_mar_b_cov.35_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.35.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.b.cov.35.AD.women.cl.25.40)


mcar_mar_b_cov.35_mean.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.35.AD.men.cl.40.50, 
                                                       mar = vector.mean.MAR.b.cov.35.AD.men.cl.40.50)


mcar_mar_b_cov.35_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.35.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.b.cov.35.AD.women.cl.40.50)


# median


mcar_mar_b_cov.35_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.35.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.b.cov.35.AD.men.cl.15.25)


mcar_mar_b_cov.35_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.35.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.b.cov.35.AD.women.cl.15.25)


mcar_mar_b_cov.35_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.35.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.b.cov.35.AD.men.cl.25.40)


mcar_mar_b_cov.35_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.35.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.b.cov.35.AD.women.cl.25.40)


mcar_mar_b_cov.35_med.men.cl.40.50<- wilcox.test.A.B(mcar = vector.med.MCAR.cov.35.AD.men.cl.40.50, 
                                                     mar = vector.med.MAR.b.cov.35.AD.men.cl.40.50)


mcar_mar_b_cov.35_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.35.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.b.cov.35.AD.women.cl.40.50)


# standard deviation

mcar_mar_b_cov.35_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.35.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.b.cov.35.AD.men.cl.15.25)


mcar_mar_b_cov.35_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.35.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.b.cov.35.AD.women.cl.15.25)


mcar_mar_b_cov.35_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.35.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.b.cov.35.AD.men.cl.25.40)


mcar_mar_b_cov.35_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.35.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.b.cov.35.AD.women.cl.25.40)


mcar_mar_b_cov.35_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.35.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.b.cov.35.AD.men.cl.40.50)


mcar_mar_b_cov.35_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.35.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.b.cov.35.AD.women.cl.40.50)


# 40 - b





# mean

mcar_mar_b_cov.40_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.40.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.b.cov.40.AD.men.cl.15.25)


mcar_mar_b_cov.40_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.40.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.b.cov.40.AD.women.cl.15.25)


mcar_mar_b_cov.40_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.40.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.b.cov.40.AD.men.cl.25.40)


mcar_mar_b_cov.40_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.40.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.b.cov.40.AD.women.cl.25.40)


mcar_mar_b_cov.40_mean.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.40.AD.men.cl.40.50, 
                                                       mar = vector.mean.MAR.b.cov.40.AD.men.cl.40.50)


mcar_mar_b_cov.40_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.40.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.b.cov.40.AD.women.cl.40.50)


# median


mcar_mar_b_cov.40_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.40.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.b.cov.40.AD.men.cl.15.25)


mcar_mar_b_cov.40_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.40.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.b.cov.40.AD.women.cl.15.25)


mcar_mar_b_cov.40_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.40.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.b.cov.40.AD.men.cl.25.40)


mcar_mar_b_cov.40_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.40.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.b.cov.40.AD.women.cl.25.40)


mcar_mar_b_cov.40_med.men.cl.40.50<- wilcox.test.A.B(mcar = vector.med.MCAR.cov.40.AD.men.cl.40.50, 
                                                     mar = vector.med.MAR.b.cov.40.AD.men.cl.40.50)


mcar_mar_b_cov.40_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.40.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.b.cov.40.AD.women.cl.40.50)


# standard deviation

mcar_mar_b_cov.40_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.40.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.b.cov.40.AD.men.cl.15.25)


mcar_mar_b_cov.40_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.40.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.b.cov.40.AD.women.cl.15.25)


mcar_mar_b_cov.40_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.40.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.b.cov.40.AD.men.cl.25.40)


mcar_mar_b_cov.40_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.40.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.b.cov.40.AD.women.cl.25.40)


mcar_mar_b_cov.40_sd.men.cl.40.50<- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.40.AD.men.cl.40.50, 
                                                    mar = vector.sd.MAR.b.cov.40.AD.men.cl.40.50)


mcar_mar_b_cov.40_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.40.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.b.cov.40.AD.women.cl.40.50)



# 45 - b





# mean

mcar_mar_b_cov.45_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.45.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.b.cov.45.AD.men.cl.15.25)


mcar_mar_b_cov.45_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.45.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.b.cov.45.AD.women.cl.15.25)


mcar_mar_b_cov.45_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.45.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.b.cov.45.AD.men.cl.25.40)


mcar_mar_b_cov.45_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.45.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.b.cov.45.AD.women.cl.25.40)


mcar_mar_b_cov.45_mean.men.cl.40.50<- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.45.AD.men.cl.40.50, 
                                                      mar = vector.mean.MAR.b.cov.45.AD.men.cl.40.50)


mcar_mar_b_cov.45_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.45.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.b.cov.45.AD.women.cl.40.50)


# median


mcar_mar_b_cov.45_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.45.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.b.cov.45.AD.men.cl.15.25)


mcar_mar_b_cov.45_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.45.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.b.cov.45.AD.women.cl.15.25)


mcar_mar_b_cov.45_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.45.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.b.cov.45.AD.men.cl.25.40)


mcar_mar_b_cov.45_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.45.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.b.cov.45.AD.women.cl.25.40)


mcar_mar_b_cov.45_med.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.45.AD.men.cl.40.50, 
                                                      mar = vector.med.MAR.b.cov.45.AD.men.cl.40.50)


mcar_mar_b_cov.45_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.45.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.b.cov.45.AD.women.cl.40.50)


# standard deviation

mcar_mar_b_cov.45_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.45.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.b.cov.45.AD.men.cl.15.25)


mcar_mar_b_cov.45_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.45.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.b.cov.45.AD.women.cl.15.25)


mcar_mar_b_cov.45_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.45.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.b.cov.45.AD.men.cl.25.40)


mcar_mar_b_cov.45_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.45.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.b.cov.45.AD.women.cl.25.40)


mcar_mar_b_cov.45_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.45.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.b.cov.45.AD.men.cl.40.50)


mcar_mar_b_cov.45_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.45.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.b.cov.45.AD.women.cl.40.50)



# 50 - b





# mean

mcar_mar_b_cov.50_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.50.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.b.cov.50.AD.men.cl.15.25)


mcar_mar_b_cov.50_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.50.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.b.cov.50.AD.women.cl.15.25)


mcar_mar_b_cov.50_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.50.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.b.cov.50.AD.men.cl.25.40)


mcar_mar_b_cov.50_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.50.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.b.cov.50.AD.women.cl.25.40)


mcar_mar_b_cov.50_mean.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.50.AD.men.cl.40.50, 
                                                       mar = vector.mean.MAR.b.cov.50.AD.men.cl.40.50)


mcar_mar_b_cov.50_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.50.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.b.cov.50.AD.women.cl.40.50)


# median


mcar_mar_b_cov.50_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.50.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.b.cov.50.AD.men.cl.15.25)


mcar_mar_b_cov.50_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.50.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.b.cov.50.AD.women.cl.15.25)


mcar_mar_b_cov.50_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.50.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.b.cov.50.AD.men.cl.25.40)


mcar_mar_b_cov.50_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.50.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.b.cov.50.AD.women.cl.25.40)


mcar_mar_b_cov.50_med.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.50.AD.men.cl.40.50, 
                                                      mar = vector.med.MAR.b.cov.50.AD.men.cl.40.50)


mcar_mar_b_cov.50_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.50.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.b.cov.50.AD.women.cl.40.50)


# standard deviation

mcar_mar_b_cov.50_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.50.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.b.cov.50.AD.men.cl.15.25)


mcar_mar_b_cov.50_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.50.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.b.cov.50.AD.women.cl.15.25)


mcar_mar_b_cov.50_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.50.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.b.cov.50.AD.men.cl.25.40)


mcar_mar_b_cov.50_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.50.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.b.cov.50.AD.women.cl.25.40)


mcar_mar_b_cov.50_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.50.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.b.cov.50.AD.men.cl.40.50)


mcar_mar_b_cov.50_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.50.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.b.cov.50.AD.women.cl.40.50)



# 55 - b




# mean

mcar_mar_b_cov.55_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.55.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.b.cov.55.AD.men.cl.15.25)


mcar_mar_b_cov.55_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.55.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.b.cov.55.AD.women.cl.15.25)


mcar_mar_b_cov.55_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.55.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.b.cov.55.AD.men.cl.25.40)


mcar_mar_b_cov.55_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.55.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.b.cov.55.AD.women.cl.25.40)


mcar_mar_b_cov.55_mean.men.cl.40.50<- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.55.AD.men.cl.40.50, 
                                                      mar = vector.mean.MAR.b.cov.55.AD.men.cl.40.50)


mcar_mar_b_cov.55_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.55.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.b.cov.55.AD.women.cl.40.50)


# median


mcar_mar_b_cov.55_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.55.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.b.cov.55.AD.men.cl.15.25)


mcar_mar_b_cov.55_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.55.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.b.cov.55.AD.women.cl.15.25)


mcar_mar_b_cov.55_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.55.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.b.cov.55.AD.men.cl.25.40)


mcar_mar_b_cov.55_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.55.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.b.cov.55.AD.women.cl.25.40)


mcar_mar_b_cov.55_med.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.55.AD.men.cl.40.50, 
                                                      mar = vector.med.MAR.b.cov.55.AD.men.cl.40.50)


mcar_mar_b_cov.55_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.55.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.b.cov.55.AD.women.cl.40.50)


# standard deviation

mcar_mar_b_cov.55_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.55.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.b.cov.55.AD.men.cl.15.25)


mcar_mar_b_cov.55_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.55.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.b.cov.55.AD.women.cl.15.25)


mcar_mar_b_cov.55_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.55.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.b.cov.55.AD.men.cl.25.40)


mcar_mar_b_cov.55_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.55.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.b.cov.55.AD.women.cl.25.40)


mcar_mar_b_cov.55_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.55.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.b.cov.55.AD.men.cl.40.50)


mcar_mar_b_cov.55_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.55.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.b.cov.55.AD.women.cl.40.50)



# 60 - b





# mean

mcar_mar_b_cov.60_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.60.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.b.cov.60.AD.men.cl.15.25)


mcar_mar_b_cov.60_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.60.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.b.cov.60.AD.women.cl.15.25)


mcar_mar_b_cov.60_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.60.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.b.cov.60.AD.men.cl.25.40)


mcar_mar_b_cov.60_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.60.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.b.cov.60.AD.women.cl.25.40)


mcar_mar_b_cov.60_mean.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.60.AD.men.cl.40.50, 
                                                       mar = vector.mean.MAR.b.cov.60.AD.men.cl.40.50)


mcar_mar_b_cov.60_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.60.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.b.cov.60.AD.women.cl.40.50)


# median


mcar_mar_b_cov.60_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.60.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.b.cov.60.AD.men.cl.15.25)


mcar_mar_b_cov.60_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.60.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.b.cov.60.AD.women.cl.15.25)


mcar_mar_b_cov.60_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.60.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.b.cov.60.AD.men.cl.25.40)


mcar_mar_b_cov.60_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.60.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.b.cov.60.AD.women.cl.25.40)


mcar_mar_b_cov.60_med.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.60.AD.men.cl.40.50, 
                                                      mar = vector.med.MAR.b.cov.60.AD.men.cl.40.50)


mcar_mar_b_cov.60_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.60.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.b.cov.60.AD.women.cl.40.50)


# standard deviation

mcar_mar_b_cov.60_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.60.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.b.cov.60.AD.men.cl.15.25)


mcar_mar_b_cov.60_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.60.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.b.cov.60.AD.women.cl.15.25)


mcar_mar_b_cov.60_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.60.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.b.cov.60.AD.men.cl.25.40)


mcar_mar_b_cov.60_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.60.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.b.cov.60.AD.women.cl.25.40)


mcar_mar_b_cov.60_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.60.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.b.cov.60.AD.men.cl.40.50)


mcar_mar_b_cov.60_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.60.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.b.cov.60.AD.women.cl.40.50)



# 65 - b




# mean

mcar_mar_b_cov.65_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.65.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.b.cov.65.AD.men.cl.15.25)


mcar_mar_b_cov.65_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.65.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.b.cov.65.AD.women.cl.15.25)


mcar_mar_b_cov.65_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.65.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.b.cov.65.AD.men.cl.25.40)


mcar_mar_b_cov.65_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.65.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.b.cov.65.AD.women.cl.25.40)


mcar_mar_b_cov.65_mean.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.65.AD.men.cl.40.50, 
                                                       mar = vector.mean.MAR.b.cov.65.AD.men.cl.40.50)


mcar_mar_b_cov.65_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.65.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.b.cov.65.AD.women.cl.40.50)


# median


mcar_mar_b_cov.65_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.65.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.b.cov.65.AD.men.cl.15.25)


mcar_mar_b_cov.65_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.65.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.b.cov.65.AD.women.cl.15.25)


mcar_mar_b_cov.65_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.65.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.b.cov.65.AD.men.cl.25.40)


mcar_mar_b_cov.65_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.65.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.b.cov.65.AD.women.cl.25.40)


mcar_mar_b_cov.65_med.men.cl.40.50<- wilcox.test.A.B(mcar = vector.med.MCAR.cov.65.AD.men.cl.40.50, 
                                                     mar = vector.med.MAR.b.cov.65.AD.men.cl.40.50)


mcar_mar_b_cov.65_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.65.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.b.cov.65.AD.women.cl.40.50)


# standard deviation

mcar_mar_b_cov.65_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.65.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.b.cov.65.AD.men.cl.15.25)


mcar_mar_b_cov.65_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.65.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.b.cov.65.AD.women.cl.15.25)


mcar_mar_b_cov.65_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.65.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.b.cov.65.AD.men.cl.25.40)


mcar_mar_b_cov.65_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.65.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.b.cov.65.AD.women.cl.25.40)


mcar_mar_b_cov.65_sd.men.cl.40.50<- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.65.AD.men.cl.40.50, 
                                                    mar = vector.sd.MAR.b.cov.65.AD.men.cl.40.50)


mcar_mar_b_cov.65_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.65.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.b.cov.65.AD.women.cl.40.50)



# 70 - b





# mean

mcar_mar_b_cov.70_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.70.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.b.cov.70.AD.men.cl.15.25)


mcar_mar_b_cov.70_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.70.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.b.cov.70.AD.women.cl.15.25)


mcar_mar_b_cov.70_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.70.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.b.cov.70.AD.men.cl.25.40)


mcar_mar_b_cov.70_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.70.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.b.cov.70.AD.women.cl.25.40)


mcar_mar_b_cov.70_mean.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.70.AD.men.cl.40.50, 
                                                       mar = vector.mean.MAR.b.cov.70.AD.men.cl.40.50)


mcar_mar_b_cov.70_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.70.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.b.cov.70.AD.women.cl.40.50)


# median


mcar_mar_b_cov.70_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.70.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.b.cov.70.AD.men.cl.15.25)


mcar_mar_b_cov.70_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.70.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.b.cov.70.AD.women.cl.15.25)


mcar_mar_b_cov.70_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.70.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.b.cov.70.AD.men.cl.25.40)


mcar_mar_b_cov.70_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.70.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.b.cov.70.AD.women.cl.25.40)


mcar_mar_b_cov.70_med.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.70.AD.men.cl.40.50, 
                                                      mar = vector.med.MAR.b.cov.70.AD.men.cl.40.50)


mcar_mar_b_cov.70_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.70.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.b.cov.70.AD.women.cl.40.50)


# standard deviation

mcar_mar_b_cov.70_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.70.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.b.cov.70.AD.men.cl.15.25)


mcar_mar_b_cov.70_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.70.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.b.cov.70.AD.women.cl.15.25)


mcar_mar_b_cov.70_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.70.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.b.cov.70.AD.men.cl.25.40)


mcar_mar_b_cov.70_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.70.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.b.cov.70.AD.women.cl.25.40)


mcar_mar_b_cov.70_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.70.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.b.cov.70.AD.men.cl.40.50)


mcar_mar_b_cov.70_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.70.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.b.cov.70.AD.women.cl.40.50)



# 75 - b




# mean

mcar_mar_b_cov.75_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.75.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.b.cov.75.AD.men.cl.15.25)


mcar_mar_b_cov.75_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.75.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.b.cov.75.AD.women.cl.15.25)


mcar_mar_b_cov.75_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.75.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.b.cov.75.AD.men.cl.25.40)


mcar_mar_b_cov.75_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.75.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.b.cov.75.AD.women.cl.25.40)


mcar_mar_b_cov.75_mean.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.75.AD.men.cl.40.50, 
                                                       mar = vector.mean.MAR.b.cov.75.AD.men.cl.40.50)


mcar_mar_b_cov.75_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.75.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.b.cov.75.AD.women.cl.40.50)


# median


mcar_mar_b_cov.75_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.75.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.b.cov.75.AD.men.cl.15.25)


mcar_mar_b_cov.75_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.75.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.b.cov.75.AD.women.cl.15.25)


mcar_mar_b_cov.75_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.75.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.b.cov.75.AD.men.cl.25.40)


mcar_mar_b_cov.75_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.75.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.b.cov.75.AD.women.cl.25.40)


mcar_mar_b_cov.75_med.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.75.AD.men.cl.40.50, 
                                                      mar = vector.med.MAR.b.cov.75.AD.men.cl.40.50)


mcar_mar_b_cov.75_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.75.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.b.cov.75.AD.women.cl.40.50)


# standard deviation

mcar_mar_b_cov.75_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.75.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.b.cov.75.AD.men.cl.15.25)


mcar_mar_b_cov.75_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.75.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.b.cov.75.AD.women.cl.15.25)


mcar_mar_b_cov.75_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.75.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.b.cov.75.AD.men.cl.25.40)


mcar_mar_b_cov.75_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.75.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.b.cov.75.AD.women.cl.25.40)


mcar_mar_b_cov.75_sd.men.cl.40.50<- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.75.AD.men.cl.40.50, 
                                                    mar = vector.sd.MAR.b.cov.75.AD.men.cl.40.50)


mcar_mar_b_cov.75_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.75.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.b.cov.75.AD.women.cl.40.50)



# 80 - b




# mean

mcar_mar_b_cov.80_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.80.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.b.cov.80.AD.men.cl.15.25)


mcar_mar_b_cov.80_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.80.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.b.cov.80.AD.women.cl.15.25)


mcar_mar_b_cov.80_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.80.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.b.cov.80.AD.men.cl.25.40)


mcar_mar_b_cov.80_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.80.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.b.cov.80.AD.women.cl.25.40)


mcar_mar_b_cov.80_mean.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.80.AD.men.cl.40.50, 
                                                       mar = vector.mean.MAR.b.cov.80.AD.men.cl.40.50)


mcar_mar_b_cov.80_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.80.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.b.cov.80.AD.women.cl.40.50)


# median


mcar_mar_b_cov.80_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.80.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.b.cov.80.AD.men.cl.15.25)


mcar_mar_b_cov.80_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.80.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.b.cov.80.AD.women.cl.15.25)


mcar_mar_b_cov.80_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.80.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.b.cov.80.AD.men.cl.25.40)


mcar_mar_b_cov.80_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.80.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.b.cov.80.AD.women.cl.25.40)


mcar_mar_b_cov.80_med.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.80.AD.men.cl.40.50, 
                                                      mar = vector.med.MAR.b.cov.80.AD.men.cl.40.50)


mcar_mar_b_cov.80_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.80.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.b.cov.80.AD.women.cl.40.50)


# standard deviation

mcar_mar_b_cov.80_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.80.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.b.cov.80.AD.men.cl.15.25)


mcar_mar_b_cov.80_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.80.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.b.cov.80.AD.women.cl.15.25)


mcar_mar_b_cov.80_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.80.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.b.cov.80.AD.men.cl.25.40)


mcar_mar_b_cov.80_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.80.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.b.cov.80.AD.women.cl.25.40)


mcar_mar_b_cov.80_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.80.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.b.cov.80.AD.men.cl.40.50)


mcar_mar_b_cov.80_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.80.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.b.cov.80.AD.women.cl.40.50)



# 85 - b




# mean

mcar_mar_b_cov.85_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.85.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.b.cov.85.AD.men.cl.15.25)


mcar_mar_b_cov.85_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.85.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.b.cov.85.AD.women.cl.15.25)


mcar_mar_b_cov.85_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.85.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.b.cov.85.AD.men.cl.25.40)


mcar_mar_b_cov.85_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.85.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.b.cov.85.AD.women.cl.25.40)


mcar_mar_b_cov.85_mean.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.85.AD.men.cl.40.50, 
                                                       mar = vector.mean.MAR.b.cov.85.AD.men.cl.40.50)


mcar_mar_b_cov.85_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.85.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.b.cov.85.AD.women.cl.40.50)


# median


mcar_mar_b_cov.85_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.85.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.b.cov.85.AD.men.cl.15.25)


mcar_mar_b_cov.85_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.85.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.b.cov.85.AD.women.cl.15.25)


mcar_mar_b_cov.85_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.85.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.b.cov.85.AD.men.cl.25.40)


mcar_mar_b_cov.85_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.85.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.b.cov.85.AD.women.cl.25.40)


mcar_mar_b_cov.85_med.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.85.AD.men.cl.40.50, 
                                                      mar = vector.med.MAR.b.cov.85.AD.men.cl.40.50)


mcar_mar_b_cov.85_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.85.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.b.cov.85.AD.women.cl.40.50)


# standard deviation

mcar_mar_b_cov.85_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.85.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.b.cov.85.AD.men.cl.15.25)


mcar_mar_b_cov.85_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.85.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.b.cov.85.AD.women.cl.15.25)


mcar_mar_b_cov.85_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.85.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.b.cov.85.AD.men.cl.25.40)


mcar_mar_b_cov.85_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.85.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.b.cov.85.AD.women.cl.25.40)


mcar_mar_b_cov.85_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.85.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.b.cov.85.AD.men.cl.40.50)


mcar_mar_b_cov.85_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.85.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.b.cov.85.AD.women.cl.40.50)



# 90 - b




# mean

mcar_mar_b_cov.90_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.90.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.b.cov.90.AD.men.cl.15.25)


mcar_mar_b_cov.90_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.90.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.b.cov.90.AD.women.cl.15.25)


mcar_mar_b_cov.90_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.90.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.b.cov.90.AD.men.cl.25.40)


mcar_mar_b_cov.90_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.90.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.b.cov.90.AD.women.cl.25.40)


mcar_mar_b_cov.90_mean.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.90.AD.men.cl.40.50, 
                                                       mar = vector.mean.MAR.b.cov.90.AD.men.cl.40.50)


mcar_mar_b_cov.90_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.90.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.b.cov.90.AD.women.cl.40.50)


# median


mcar_mar_b_cov.90_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.90.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.b.cov.90.AD.men.cl.15.25)


mcar_mar_b_cov.90_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.90.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.b.cov.90.AD.women.cl.15.25)


mcar_mar_b_cov.90_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.90.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.b.cov.90.AD.men.cl.25.40)


mcar_mar_b_cov.90_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.90.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.b.cov.90.AD.women.cl.25.40)


mcar_mar_b_cov.90_med.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.90.AD.men.cl.40.50, 
                                                      mar = vector.med.MAR.b.cov.90.AD.men.cl.40.50)


mcar_mar_b_cov.90_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.90.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.b.cov.90.AD.women.cl.40.50)


# standard deviation

mcar_mar_b_cov.90_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.90.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.b.cov.90.AD.men.cl.15.25)


mcar_mar_b_cov.90_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.90.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.b.cov.90.AD.women.cl.15.25)


mcar_mar_b_cov.90_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.90.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.b.cov.90.AD.men.cl.25.40)


mcar_mar_b_cov.90_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.90.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.b.cov.90.AD.women.cl.25.40)


mcar_mar_b_cov.90_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.90.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.b.cov.90.AD.men.cl.40.50)


mcar_mar_b_cov.90_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.90.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.b.cov.90.AD.women.cl.40.50)



# 95 - b




# mean

mcar_mar_b_cov.95_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.95.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.b.cov.95.AD.men.cl.15.25)


mcar_mar_b_cov.95_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.95.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.b.cov.95.AD.women.cl.15.25)


mcar_mar_b_cov.95_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.95.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.b.cov.95.AD.men.cl.25.40)


mcar_mar_b_cov.95_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.95.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.b.cov.95.AD.women.cl.25.40)


mcar_mar_b_cov.95_mean.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.95.AD.men.cl.40.50, 
                                                       mar = vector.mean.MAR.b.cov.95.AD.men.cl.40.50)


mcar_mar_b_cov.95_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.95.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.b.cov.95.AD.women.cl.40.50)


# median


mcar_mar_b_cov.95_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.95.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.b.cov.95.AD.men.cl.15.25)


mcar_mar_b_cov.95_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.95.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.b.cov.95.AD.women.cl.15.25)


mcar_mar_b_cov.95_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.95.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.b.cov.95.AD.men.cl.25.40)


mcar_mar_b_cov.95_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.95.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.b.cov.95.AD.women.cl.25.40)


mcar_mar_b_cov.95_med.men.cl.40.50<- wilcox.test.A.B(mcar = vector.med.MCAR.cov.95.AD.men.cl.40.50, 
                                                     mar = vector.med.MAR.b.cov.95.AD.men.cl.40.50)


mcar_mar_b_cov.95_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.95.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.b.cov.95.AD.women.cl.40.50)


# standard deviation

mcar_mar_b_cov.95_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.95.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.b.cov.95.AD.men.cl.15.25)


mcar_mar_b_cov.95_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.95.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.b.cov.95.AD.women.cl.15.25)


mcar_mar_b_cov.95_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.95.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.b.cov.95.AD.men.cl.25.40)


mcar_mar_b_cov.95_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.95.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.b.cov.95.AD.women.cl.25.40)


mcar_mar_b_cov.95_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.95.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.b.cov.95.AD.men.cl.40.50)


mcar_mar_b_cov.95_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.95.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.b.cov.95.AD.women.cl.40.50)




# Age difference MCAR and MAR(c) -------


# mean

mcar_mar_c_cov.35_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.35.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.c.cov.35.AD.men.cl.15.25)


mcar_mar_c_cov.35_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.35.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.c.cov.35.AD.women.cl.15.25)


mcar_mar_c_cov.35_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.35.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.c.cov.35.AD.men.cl.25.40)


mcar_mar_c_cov.35_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.35.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.c.cov.35.AD.women.cl.25.40)


mcar_mar_c_cov.35_mean.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.35.AD.men.cl.40.50, 
                                                       mar = vector.mean.MAR.c.cov.35.AD.men.cl.40.50)


mcar_mar_c_cov.35_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.35.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.c.cov.35.AD.women.cl.40.50)


# median


mcar_mar_c_cov.35_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.35.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.c.cov.35.AD.men.cl.15.25)


mcar_mar_c_cov.35_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.35.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.c.cov.35.AD.women.cl.15.25)


mcar_mar_c_cov.35_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.35.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.c.cov.35.AD.men.cl.25.40)


mcar_mar_c_cov.35_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.35.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.c.cov.35.AD.women.cl.25.40)


mcar_mar_c_cov.35_med.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.35.AD.men.cl.40.50, 
                                                      mar = vector.med.MAR.c.cov.35.AD.men.cl.40.50)


mcar_mar_c_cov.35_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.35.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.c.cov.35.AD.women.cl.40.50)


# standard deviation

mcar_mar_c_cov.35_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.35.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.c.cov.35.AD.men.cl.15.25)


mcar_mar_c_cov.35_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.35.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.c.cov.35.AD.women.cl.15.25)


mcar_mar_c_cov.35_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.35.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.c.cov.35.AD.men.cl.25.40)


mcar_mar_c_cov.35_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.35.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.c.cov.35.AD.women.cl.25.40)


mcar_mar_c_cov.35_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.35.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.c.cov.35.AD.men.cl.40.50)


mcar_mar_c_cov.35_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.35.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.c.cov.35.AD.women.cl.40.50)


# 40 - c





# mean

mcar_mar_c_cov.40_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.40.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.c.cov.40.AD.men.cl.15.25)


mcar_mar_c_cov.40_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.40.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.c.cov.40.AD.women.cl.15.25)


mcar_mar_c_cov.40_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.40.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.c.cov.40.AD.men.cl.25.40)


mcar_mar_c_cov.40_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.40.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.c.cov.40.AD.women.cl.25.40)


mcar_mar_c_cov.40_mean.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.40.AD.men.cl.40.50, 
                                                       mar = vector.mean.MAR.c.cov.40.AD.men.cl.40.50)


mcar_mar_c_cov.40_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.40.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.c.cov.40.AD.women.cl.40.50)


# median


mcar_mar_c_cov.40_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.40.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.c.cov.40.AD.men.cl.15.25)


mcar_mar_c_cov.40_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.40.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.c.cov.40.AD.women.cl.15.25)


mcar_mar_c_cov.40_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.40.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.c.cov.40.AD.men.cl.25.40)


mcar_mar_c_cov.40_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.40.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.c.cov.40.AD.women.cl.25.40)


mcar_mar_c_cov.40_med.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.40.AD.men.cl.40.50, 
                                                      mar = vector.med.MAR.c.cov.40.AD.men.cl.40.50)


mcar_mar_c_cov.40_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.40.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.c.cov.40.AD.women.cl.40.50)


# standard deviation

mcar_mar_c_cov.40_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.40.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.c.cov.40.AD.men.cl.15.25)


mcar_mar_c_cov.40_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.40.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.c.cov.40.AD.women.cl.15.25)


mcar_mar_c_cov.40_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.40.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.c.cov.40.AD.men.cl.25.40)


mcar_mar_c_cov.40_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.40.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.c.cov.40.AD.women.cl.25.40)


mcar_mar_c_cov.40_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.40.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.c.cov.40.AD.men.cl.40.50)


mcar_mar_c_cov.40_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.40.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.c.cov.40.AD.women.cl.40.50)



# 45 - c





# mean

mcar_mar_c_cov.45_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.45.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.c.cov.45.AD.men.cl.15.25)


mcar_mar_c_cov.45_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.45.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.c.cov.45.AD.women.cl.15.25)


mcar_mar_c_cov.45_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.45.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.c.cov.45.AD.men.cl.25.40)


mcar_mar_c_cov.45_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.45.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.c.cov.45.AD.women.cl.25.40)


mcar_mar_c_cov.45_mean.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.45.AD.men.cl.40.50, 
                                                       mar = vector.mean.MAR.c.cov.45.AD.men.cl.40.50)


mcar_mar_c_cov.45_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.45.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.c.cov.45.AD.women.cl.40.50)


# median


mcar_mar_c_cov.45_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.45.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.c.cov.45.AD.men.cl.15.25)


mcar_mar_c_cov.45_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.45.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.c.cov.45.AD.women.cl.15.25)


mcar_mar_c_cov.45_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.45.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.c.cov.45.AD.men.cl.25.40)


mcar_mar_c_cov.45_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.45.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.c.cov.45.AD.women.cl.25.40)


mcar_mar_c_cov.45_med.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.45.AD.men.cl.40.50, 
                                                      mar = vector.med.MAR.c.cov.45.AD.men.cl.40.50)


mcar_mar_c_cov.45_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.45.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.c.cov.45.AD.women.cl.40.50)


# standard deviation

mcar_mar_c_cov.45_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.45.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.c.cov.45.AD.men.cl.15.25)


mcar_mar_c_cov.45_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.45.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.c.cov.45.AD.women.cl.15.25)


mcar_mar_c_cov.45_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.45.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.c.cov.45.AD.men.cl.25.40)


mcar_mar_c_cov.45_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.45.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.c.cov.45.AD.women.cl.25.40)


mcar_mar_c_cov.45_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.45.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.c.cov.45.AD.men.cl.40.50)


mcar_mar_c_cov.45_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.45.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.c.cov.45.AD.women.cl.40.50)



# 50 - c





# mean

mcar_mar_c_cov.50_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.50.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.c.cov.50.AD.men.cl.15.25)


mcar_mar_c_cov.50_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.50.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.c.cov.50.AD.women.cl.15.25)


mcar_mar_c_cov.50_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.50.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.c.cov.50.AD.men.cl.25.40)


mcar_mar_c_cov.50_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.50.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.c.cov.50.AD.women.cl.25.40)


mcar_mar_c_cov.50_mean.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.50.AD.men.cl.40.50, 
                                                       mar = vector.mean.MAR.c.cov.50.AD.men.cl.40.50)


mcar_mar_c_cov.50_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.50.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.c.cov.50.AD.women.cl.40.50)


# median


mcar_mar_c_cov.50_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.50.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.c.cov.50.AD.men.cl.15.25)


mcar_mar_c_cov.50_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.50.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.c.cov.50.AD.women.cl.15.25)


mcar_mar_c_cov.50_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.50.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.c.cov.50.AD.men.cl.25.40)


mcar_mar_c_cov.50_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.50.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.c.cov.50.AD.women.cl.25.40)


mcar_mar_c_cov.50_med.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.50.AD.men.cl.40.50, 
                                                      mar = vector.med.MAR.c.cov.50.AD.men.cl.40.50)


mcar_mar_c_cov.50_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.50.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.c.cov.50.AD.women.cl.40.50)


# standard deviation

mcar_mar_c_cov.50_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.50.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.c.cov.50.AD.men.cl.15.25)


mcar_mar_c_cov.50_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.50.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.c.cov.50.AD.women.cl.15.25)


mcar_mar_c_cov.50_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.50.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.c.cov.50.AD.men.cl.25.40)


mcar_mar_c_cov.50_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.50.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.c.cov.50.AD.women.cl.25.40)


mcar_mar_c_cov.50_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.50.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.c.cov.50.AD.men.cl.40.50)


mcar_mar_c_cov.50_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.50.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.c.cov.50.AD.women.cl.40.50)



# 55 - c




# mean

mcar_mar_c_cov.55_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.55.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.c.cov.55.AD.men.cl.15.25)


mcar_mar_c_cov.55_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.55.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.c.cov.55.AD.women.cl.15.25)


mcar_mar_c_cov.55_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.55.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.c.cov.55.AD.men.cl.25.40)


mcar_mar_c_cov.55_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.55.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.c.cov.55.AD.women.cl.25.40)


mcar_mar_c_cov.55_mean.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.55.AD.men.cl.40.50, 
                                                       mar = vector.mean.MAR.c.cov.55.AD.men.cl.40.50)


mcar_mar_c_cov.55_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.55.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.c.cov.55.AD.women.cl.40.50)


# median


mcar_mar_c_cov.55_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.55.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.c.cov.55.AD.men.cl.15.25)


mcar_mar_c_cov.55_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.55.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.c.cov.55.AD.women.cl.15.25)


mcar_mar_c_cov.55_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.55.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.c.cov.55.AD.men.cl.25.40)


mcar_mar_c_cov.55_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.55.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.c.cov.55.AD.women.cl.25.40)


mcar_mar_c_cov.55_med.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.55.AD.men.cl.40.50, 
                                                      mar = vector.med.MAR.c.cov.55.AD.men.cl.40.50)


mcar_mar_c_cov.55_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.55.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.c.cov.55.AD.women.cl.40.50)


# standard deviation

mcar_mar_c_cov.55_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.55.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.c.cov.55.AD.men.cl.15.25)


mcar_mar_c_cov.55_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.55.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.c.cov.55.AD.women.cl.15.25)


mcar_mar_c_cov.55_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.55.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.c.cov.55.AD.men.cl.25.40)


mcar_mar_c_cov.55_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.55.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.c.cov.55.AD.women.cl.25.40)


mcar_mar_c_cov.55_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.55.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.c.cov.55.AD.men.cl.40.50)


mcar_mar_c_cov.55_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.55.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.c.cov.55.AD.women.cl.40.50)



# 60 - c





# mean

mcar_mar_c_cov.60_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.60.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.c.cov.60.AD.men.cl.15.25)


mcar_mar_c_cov.60_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.60.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.c.cov.60.AD.women.cl.15.25)


mcar_mar_c_cov.60_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.60.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.c.cov.60.AD.men.cl.25.40)


mcar_mar_c_cov.60_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.60.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.c.cov.60.AD.women.cl.25.40)


mcar_mar_c_cov.60_mean.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.60.AD.men.cl.40.50, 
                                                       mar = vector.mean.MAR.c.cov.60.AD.men.cl.40.50)


mcar_mar_c_cov.60_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.60.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.c.cov.60.AD.women.cl.40.50)


# median


mcar_mar_c_cov.60_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.60.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.c.cov.60.AD.men.cl.15.25)


mcar_mar_c_cov.60_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.60.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.c.cov.60.AD.women.cl.15.25)


mcar_mar_c_cov.60_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.60.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.c.cov.60.AD.men.cl.25.40)


mcar_mar_c_cov.60_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.60.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.c.cov.60.AD.women.cl.25.40)


mcar_mar_c_cov.60_med.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.60.AD.men.cl.40.50, 
                                                      mar = vector.med.MAR.c.cov.60.AD.men.cl.40.50)


mcar_mar_c_cov.60_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.60.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.c.cov.60.AD.women.cl.40.50)


# standard deviation

mcar_mar_c_cov.60_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.60.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.c.cov.60.AD.men.cl.15.25)


mcar_mar_c_cov.60_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.60.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.c.cov.60.AD.women.cl.15.25)


mcar_mar_c_cov.60_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.60.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.c.cov.60.AD.men.cl.25.40)


mcar_mar_c_cov.60_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.60.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.c.cov.60.AD.women.cl.25.40)


mcar_mar_c_cov.60_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.60.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.c.cov.60.AD.men.cl.40.50)


mcar_mar_c_cov.60_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.60.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.c.cov.60.AD.women.cl.40.50)



# 65 - c




# mean

mcar_mar_c_cov.65_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.65.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.c.cov.65.AD.men.cl.15.25)


mcar_mar_c_cov.65_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.65.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.c.cov.65.AD.women.cl.15.25)


mcar_mar_c_cov.65_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.65.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.c.cov.65.AD.men.cl.25.40)


mcar_mar_c_cov.65_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.65.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.c.cov.65.AD.women.cl.25.40)


mcar_mar_c_cov.65_mean.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.65.AD.men.cl.40.50, 
                                                       mar = vector.mean.MAR.c.cov.65.AD.men.cl.40.50)


mcar_mar_c_cov.65_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.65.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.c.cov.65.AD.women.cl.40.50)


# median


mcar_mar_c_cov.65_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.65.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.c.cov.65.AD.men.cl.15.25)


mcar_mar_c_cov.65_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.65.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.c.cov.65.AD.women.cl.15.25)


mcar_mar_c_cov.65_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.65.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.c.cov.65.AD.men.cl.25.40)


mcar_mar_c_cov.65_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.65.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.c.cov.65.AD.women.cl.25.40)


mcar_mar_c_cov.65_med.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.65.AD.men.cl.40.50, 
                                                      mar = vector.med.MAR.c.cov.65.AD.men.cl.40.50)


mcar_mar_c_cov.65_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.65.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.c.cov.65.AD.women.cl.40.50)


# standard deviation

mcar_mar_c_cov.65_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.65.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.c.cov.65.AD.men.cl.15.25)


mcar_mar_c_cov.65_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.65.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.c.cov.65.AD.women.cl.15.25)


mcar_mar_c_cov.65_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.65.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.c.cov.65.AD.men.cl.25.40)


mcar_mar_c_cov.65_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.65.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.c.cov.65.AD.women.cl.25.40)


mcar_mar_c_cov.65_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.65.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.c.cov.65.AD.men.cl.40.50)


mcar_mar_c_cov.65_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.65.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.c.cov.65.AD.women.cl.40.50)



# 70 - c





# mean

mcar_mar_c_cov.70_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.70.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.c.cov.70.AD.men.cl.15.25)


mcar_mar_c_cov.70_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.70.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.c.cov.70.AD.women.cl.15.25)


mcar_mar_c_cov.70_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.70.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.c.cov.70.AD.men.cl.25.40)


mcar_mar_c_cov.70_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.70.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.c.cov.70.AD.women.cl.25.40)


mcar_mar_c_cov.70_mean.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.70.AD.men.cl.40.50, 
                                                       mar = vector.mean.MAR.c.cov.70.AD.men.cl.40.50)


mcar_mar_c_cov.70_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.70.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.c.cov.70.AD.women.cl.40.50)


# median


mcar_mar_c_cov.70_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.70.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.c.cov.70.AD.men.cl.15.25)


mcar_mar_c_cov.70_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.70.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.c.cov.70.AD.women.cl.15.25)


mcar_mar_c_cov.70_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.70.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.c.cov.70.AD.men.cl.25.40)


mcar_mar_c_cov.70_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.70.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.c.cov.70.AD.women.cl.25.40)


mcar_mar_c_cov.70_med.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.70.AD.men.cl.40.50, 
                                                      mar = vector.med.MAR.c.cov.70.AD.men.cl.40.50)


mcar_mar_c_cov.70_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.70.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.c.cov.70.AD.women.cl.40.50)


# standard deviation

mcar_mar_c_cov.70_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.70.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.c.cov.70.AD.men.cl.15.25)


mcar_mar_c_cov.70_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.70.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.c.cov.70.AD.women.cl.15.25)


mcar_mar_c_cov.70_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.70.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.c.cov.70.AD.men.cl.25.40)


mcar_mar_c_cov.70_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.70.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.c.cov.70.AD.women.cl.25.40)


mcar_mar_c_cov.70_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.70.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.c.cov.70.AD.men.cl.40.50)


mcar_mar_c_cov.70_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.70.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.c.cov.70.AD.women.cl.40.50)



# 75 - c




# mean

mcar_mar_c_cov.75_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.75.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.c.cov.75.AD.men.cl.15.25)


mcar_mar_c_cov.75_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.75.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.c.cov.75.AD.women.cl.15.25)


mcar_mar_c_cov.75_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.75.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.c.cov.75.AD.men.cl.25.40)


mcar_mar_c_cov.75_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.75.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.c.cov.75.AD.women.cl.25.40)


mcar_mar_c_cov.75_mean.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.75.AD.men.cl.40.50, 
                                                       mar = vector.mean.MAR.c.cov.75.AD.men.cl.40.50)


mcar_mar_c_cov.75_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.75.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.c.cov.75.AD.women.cl.40.50)


# median


mcar_mar_c_cov.75_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.75.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.c.cov.75.AD.men.cl.15.25)


mcar_mar_c_cov.75_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.75.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.c.cov.75.AD.women.cl.15.25)


mcar_mar_c_cov.75_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.75.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.c.cov.75.AD.men.cl.25.40)


mcar_mar_c_cov.75_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.75.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.c.cov.75.AD.women.cl.25.40)


mcar_mar_c_cov.75_med.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.75.AD.men.cl.40.50, 
                                                      mar = vector.med.MAR.c.cov.75.AD.men.cl.40.50)


mcar_mar_c_cov.75_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.75.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.c.cov.75.AD.women.cl.40.50)


# standard deviation

mcar_mar_c_cov.75_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.75.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.c.cov.75.AD.men.cl.15.25)


mcar_mar_c_cov.75_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.75.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.c.cov.75.AD.women.cl.15.25)


mcar_mar_c_cov.75_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.75.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.c.cov.75.AD.men.cl.25.40)


mcar_mar_c_cov.75_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.75.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.c.cov.75.AD.women.cl.25.40)


mcar_mar_c_cov.75_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.75.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.c.cov.75.AD.men.cl.40.50)


mcar_mar_c_cov.75_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.75.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.c.cov.75.AD.women.cl.40.50)



# 80 - c




# mean

mcar_mar_c_cov.80_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.80.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.c.cov.80.AD.men.cl.15.25)


mcar_mar_c_cov.80_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.80.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.c.cov.80.AD.women.cl.15.25)


mcar_mar_c_cov.80_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.80.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.c.cov.80.AD.men.cl.25.40)


mcar_mar_c_cov.80_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.80.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.c.cov.80.AD.women.cl.25.40)


mcar_mar_c_cov.80_mean.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.80.AD.men.cl.40.50, 
                                                       mar = vector.mean.MAR.c.cov.80.AD.men.cl.40.50)


mcar_mar_c_cov.80_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.80.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.c.cov.80.AD.women.cl.40.50)


# median


mcar_mar_c_cov.80_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.80.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.c.cov.80.AD.men.cl.15.25)


mcar_mar_c_cov.80_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.80.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.c.cov.80.AD.women.cl.15.25)


mcar_mar_c_cov.80_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.80.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.c.cov.80.AD.men.cl.25.40)


mcar_mar_c_cov.80_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.80.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.c.cov.80.AD.women.cl.25.40)


mcar_mar_c_cov.80_med.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.80.AD.men.cl.40.50, 
                                                      mar = vector.med.MAR.c.cov.80.AD.men.cl.40.50)


mcar_mar_c_cov.80_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.80.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.c.cov.80.AD.women.cl.40.50)


# standard deviation

mcar_mar_c_cov.80_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.80.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.c.cov.80.AD.men.cl.15.25)


mcar_mar_c_cov.80_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.80.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.c.cov.80.AD.women.cl.15.25)


mcar_mar_c_cov.80_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.80.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.c.cov.80.AD.men.cl.25.40)


mcar_mar_c_cov.80_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.80.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.c.cov.80.AD.women.cl.25.40)


mcar_mar_c_cov.80_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.80.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.c.cov.80.AD.men.cl.40.50)


mcar_mar_c_cov.80_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.80.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.c.cov.80.AD.women.cl.40.50)



# 85 - c




# mean

mcar_mar_c_cov.85_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.85.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.c.cov.85.AD.men.cl.15.25)


mcar_mar_c_cov.85_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.85.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.c.cov.85.AD.women.cl.15.25)


mcar_mar_c_cov.85_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.85.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.c.cov.85.AD.men.cl.25.40)


mcar_mar_c_cov.85_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.85.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.c.cov.85.AD.women.cl.25.40)


mcar_mar_c_cov.85_mean.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.85.AD.men.cl.40.50, 
                                                       mar = vector.mean.MAR.c.cov.85.AD.men.cl.40.50)


mcar_mar_c_cov.85_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.85.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.c.cov.85.AD.women.cl.40.50)


# median


mcar_mar_c_cov.85_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.85.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.c.cov.85.AD.men.cl.15.25)


mcar_mar_c_cov.85_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.85.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.c.cov.85.AD.women.cl.15.25)


mcar_mar_c_cov.85_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.85.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.c.cov.85.AD.men.cl.25.40)


mcar_mar_c_cov.85_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.85.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.c.cov.85.AD.women.cl.25.40)


mcar_mar_c_cov.85_med.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.85.AD.men.cl.40.50, 
                                                      mar = vector.med.MAR.c.cov.85.AD.men.cl.40.50)


mcar_mar_c_cov.85_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.85.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.c.cov.85.AD.women.cl.40.50)


# standard deviation

mcar_mar_c_cov.85_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.85.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.c.cov.85.AD.men.cl.15.25)


mcar_mar_c_cov.85_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.85.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.c.cov.85.AD.women.cl.15.25)


mcar_mar_c_cov.85_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.85.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.c.cov.85.AD.men.cl.25.40)


mcar_mar_c_cov.85_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.85.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.c.cov.85.AD.women.cl.25.40)


mcar_mar_c_cov.85_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.85.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.c.cov.85.AD.men.cl.40.50)


mcar_mar_c_cov.85_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.85.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.c.cov.85.AD.women.cl.40.50)



# 90 - c




# mean

mcar_mar_c_cov.90_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.90.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.c.cov.90.AD.men.cl.15.25)


mcar_mar_c_cov.90_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.90.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.c.cov.90.AD.women.cl.15.25)


mcar_mar_c_cov.90_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.90.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.c.cov.90.AD.men.cl.25.40)


mcar_mar_c_cov.90_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.90.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.c.cov.90.AD.women.cl.25.40)


mcar_mar_c_cov.90_mean.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.90.AD.men.cl.40.50, 
                                                       mar = vector.mean.MAR.c.cov.90.AD.men.cl.40.50)


mcar_mar_c_cov.90_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.90.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.c.cov.90.AD.women.cl.40.50)


# median


mcar_mar_c_cov.90_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.90.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.c.cov.90.AD.men.cl.15.25)


mcar_mar_c_cov.90_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.90.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.c.cov.90.AD.women.cl.15.25)


mcar_mar_c_cov.90_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.90.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.c.cov.90.AD.men.cl.25.40)


mcar_mar_c_cov.90_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.90.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.c.cov.90.AD.women.cl.25.40)


mcar_mar_c_cov.90_med.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.90.AD.men.cl.40.50, 
                                                      mar = vector.med.MAR.c.cov.90.AD.men.cl.40.50)


mcar_mar_c_cov.90_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.90.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.c.cov.90.AD.women.cl.40.50)


# standard deviation

mcar_mar_c_cov.90_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.90.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.c.cov.90.AD.men.cl.15.25)


mcar_mar_c_cov.90_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.90.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.c.cov.90.AD.women.cl.15.25)


mcar_mar_c_cov.90_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.90.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.c.cov.90.AD.men.cl.25.40)


mcar_mar_c_cov.90_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.90.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.c.cov.90.AD.women.cl.25.40)


mcar_mar_c_cov.90_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.90.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.c.cov.90.AD.men.cl.40.50)


mcar_mar_c_cov.90_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.90.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.c.cov.90.AD.women.cl.40.50)



# 95 - c




# mean

mcar_mar_c_cov.95_mean.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.95.AD.men.cl.15.25, 
                                                       mar = vector.mean.MAR.c.cov.95.AD.men.cl.15.25)


mcar_mar_c_cov.95_mean.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.95.AD.women.cl.15.25, 
                                                         mar = vector.mean.MAR.c.cov.95.AD.women.cl.15.25)


mcar_mar_c_cov.95_mean.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.95.AD.men.cl.25.40, 
                                                       mar = vector.mean.MAR.c.cov.95.AD.men.cl.25.40)


mcar_mar_c_cov.95_mean.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.95.AD.women.cl.25.40, 
                                                         mar = vector.mean.MAR.c.cov.95.AD.women.cl.25.40)


mcar_mar_c_cov.95_mean.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.95.AD.men.cl.40.50, 
                                                       mar = vector.mean.MAR.c.cov.95.AD.men.cl.40.50)


mcar_mar_c_cov.95_mean.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.mean.MCAR.cov.95.AD.women.cl.40.50, 
                                                         mar = vector.mean.MAR.c.cov.95.AD.women.cl.40.50)


# median


mcar_mar_c_cov.95_med.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.95.AD.men.cl.15.25, 
                                                      mar = vector.med.MAR.c.cov.95.AD.men.cl.15.25)


mcar_mar_c_cov.95_med.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.95.AD.women.cl.15.25, 
                                                        mar = vector.med.MAR.c.cov.95.AD.women.cl.15.25)


mcar_mar_c_cov.95_med.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.95.AD.men.cl.25.40, 
                                                      mar = vector.med.MAR.c.cov.95.AD.men.cl.25.40)


mcar_mar_c_cov.95_med.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.95.AD.women.cl.25.40, 
                                                        mar = vector.med.MAR.c.cov.95.AD.women.cl.25.40)


mcar_mar_c_cov.95_med.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.95.AD.men.cl.40.50, 
                                                      mar = vector.med.MAR.c.cov.95.AD.men.cl.40.50)


mcar_mar_c_cov.95_med.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.med.MCAR.cov.95.AD.women.cl.40.50, 
                                                        mar = vector.med.MAR.c.cov.95.AD.women.cl.40.50)


# standard deviation

mcar_mar_c_cov.95_sd.men.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.95.AD.men.cl.15.25, 
                                                     mar = vector.sd.MAR.c.cov.95.AD.men.cl.15.25)


mcar_mar_c_cov.95_sd.women.cl.15.25 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.95.AD.women.cl.15.25, 
                                                       mar = vector.sd.MAR.c.cov.95.AD.women.cl.15.25)


mcar_mar_c_cov.95_sd.men.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.95.AD.men.cl.25.40, 
                                                     mar = vector.sd.MAR.c.cov.95.AD.men.cl.25.40)


mcar_mar_c_cov.95_sd.women.cl.25.40 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.95.AD.women.cl.25.40, 
                                                       mar = vector.sd.MAR.c.cov.95.AD.women.cl.25.40)


mcar_mar_c_cov.95_sd.men.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.95.AD.men.cl.40.50, 
                                                     mar = vector.sd.MAR.c.cov.95.AD.men.cl.40.50)


mcar_mar_c_cov.95_sd.women.cl.40.50 <- wilcox.test.A.B(mcar = vector.sd.MCAR.cov.95.AD.women.cl.40.50, 
                                                       mar = vector.sd.MAR.c.cov.95.AD.women.cl.40.50)




# All together: summary table -------------



# MCAR - MAR - a proportions -------------- 



mcar_mar_a_prop.men15.25.F.15.25 <- c(mcar_mar_a_cov.35_prop.men15.25.F.15.25,
                                      mcar_mar_a_cov.40_prop.men15.25.F.15.25,
                                      mcar_mar_a_cov.45_prop.men15.25.F.15.25,
                                      mcar_mar_a_cov.50_prop.men15.25.F.15.25,
                                      mcar_mar_a_cov.55_prop.men15.25.F.15.25,
                                      mcar_mar_a_cov.60_prop.men15.25.F.15.25,
                                      mcar_mar_a_cov.65_prop.men15.25.F.15.25,
                                      mcar_mar_a_cov.70_prop.men15.25.F.15.25,
                                      mcar_mar_a_cov.75_prop.men15.25.F.15.25,
                                      mcar_mar_a_cov.80_prop.men15.25.F.15.25,
                                      mcar_mar_a_cov.85_prop.men15.25.F.15.25,
                                      mcar_mar_a_cov.90_prop.men15.25.F.15.25,
                                      mcar_mar_a_cov.95_prop.men15.25.F.15.25)



mcar_mar_a_prop.women15.25.M.15.25 <- c(mcar_mar_a_cov.35_prop.women15.25.M.15.25,
                                        mcar_mar_a_cov.40_prop.women15.25.M.15.25,
                                        mcar_mar_a_cov.45_prop.women15.25.M.15.25,
                                        mcar_mar_a_cov.50_prop.women15.25.M.15.25,
                                        mcar_mar_a_cov.55_prop.women15.25.M.15.25,
                                        mcar_mar_a_cov.60_prop.women15.25.M.15.25,
                                        mcar_mar_a_cov.65_prop.women15.25.M.15.25,
                                        mcar_mar_a_cov.70_prop.women15.25.M.15.25,
                                        mcar_mar_a_cov.75_prop.women15.25.M.15.25,
                                        mcar_mar_a_cov.80_prop.women15.25.M.15.25,
                                        mcar_mar_a_cov.85_prop.women15.25.M.15.25,
                                        mcar_mar_a_cov.90_prop.women15.25.M.15.25,
                                        mcar_mar_a_cov.95_prop.women15.25.M.15.25)

mcar_mar_a_prop.men25.40.F.15.25 <- c(mcar_mar_a_cov.35_prop.men25.40.F.15.25,
                                      mcar_mar_a_cov.40_prop.men25.40.F.15.25,
                                      mcar_mar_a_cov.45_prop.men25.40.F.15.25,
                                      mcar_mar_a_cov.50_prop.men25.40.F.15.25,
                                      mcar_mar_a_cov.55_prop.men25.40.F.15.25,
                                      mcar_mar_a_cov.60_prop.men25.40.F.15.25,
                                      mcar_mar_a_cov.65_prop.men25.40.F.15.25,
                                      mcar_mar_a_cov.70_prop.men25.40.F.15.25,
                                      mcar_mar_a_cov.75_prop.men25.40.F.15.25,
                                      mcar_mar_a_cov.80_prop.men25.40.F.15.25,
                                      mcar_mar_a_cov.85_prop.men25.40.F.15.25,
                                      mcar_mar_a_cov.90_prop.men25.40.F.15.25,
                                      mcar_mar_a_cov.95_prop.men25.40.F.15.25)


mcar_mar_a_prop.women15.25.M.25.40 <- c(mcar_mar_a_cov.35_prop.women15.25.M.25.40,
                                        mcar_mar_a_cov.40_prop.women15.25.M.25.40,
                                        mcar_mar_a_cov.45_prop.women15.25.M.25.40,
                                        mcar_mar_a_cov.50_prop.women15.25.M.25.40,
                                        mcar_mar_a_cov.55_prop.women15.25.M.25.40,
                                        mcar_mar_a_cov.60_prop.women15.25.M.25.40,
                                        mcar_mar_a_cov.65_prop.women15.25.M.25.40,
                                        mcar_mar_a_cov.70_prop.women15.25.M.25.40,
                                        mcar_mar_a_cov.75_prop.women15.25.M.25.40,
                                        mcar_mar_a_cov.80_prop.women15.25.M.25.40,
                                        mcar_mar_a_cov.85_prop.women15.25.M.25.40,
                                        mcar_mar_a_cov.90_prop.women15.25.M.25.40,
                                        mcar_mar_a_cov.95_prop.women15.25.M.25.40)

mcar_mar_a_prop.men25.40.F.25.40 <- c(mcar_mar_a_cov.35_prop.men25.40.F.25.40,
                                      mcar_mar_a_cov.40_prop.men25.40.F.25.40,
                                      mcar_mar_a_cov.45_prop.men25.40.F.25.40,
                                      mcar_mar_a_cov.50_prop.men25.40.F.25.40,
                                      mcar_mar_a_cov.55_prop.men25.40.F.25.40,
                                      mcar_mar_a_cov.60_prop.men25.40.F.25.40,
                                      mcar_mar_a_cov.65_prop.men25.40.F.25.40,
                                      mcar_mar_a_cov.70_prop.men25.40.F.25.40,
                                      mcar_mar_a_cov.75_prop.men25.40.F.25.40,
                                      mcar_mar_a_cov.80_prop.men25.40.F.25.40,
                                      mcar_mar_a_cov.85_prop.men25.40.F.25.40,
                                      mcar_mar_a_cov.90_prop.men25.40.F.25.40,
                                      mcar_mar_a_cov.95_prop.men25.40.F.25.40)


mcar_mar_a_prop.women25.40.M.25.40 <- c(mcar_mar_a_cov.35_prop.women25.40.M.25.40,
                                        mcar_mar_a_cov.40_prop.women25.40.M.25.40,
                                        mcar_mar_a_cov.45_prop.women25.40.M.25.40,
                                        mcar_mar_a_cov.50_prop.women25.40.M.25.40,
                                        mcar_mar_a_cov.55_prop.women25.40.M.25.40,
                                        mcar_mar_a_cov.60_prop.women25.40.M.25.40,
                                        mcar_mar_a_cov.65_prop.women25.40.M.25.40,
                                        mcar_mar_a_cov.70_prop.women25.40.M.25.40,
                                        mcar_mar_a_cov.75_prop.women25.40.M.25.40,
                                        mcar_mar_a_cov.80_prop.women25.40.M.25.40,
                                        mcar_mar_a_cov.85_prop.women25.40.M.25.40,
                                        mcar_mar_a_cov.90_prop.women25.40.M.25.40,
                                        mcar_mar_a_cov.95_prop.women25.40.M.25.40)


mcar_mar_a_prop.men40.50.F.15.25 <- c(mcar_mar_a_cov.35_prop.men40.50.F.15.25,
                                      mcar_mar_a_cov.40_prop.men40.50.F.15.25,
                                      mcar_mar_a_cov.45_prop.men40.50.F.15.25,
                                      mcar_mar_a_cov.50_prop.men40.50.F.15.25,
                                      mcar_mar_a_cov.55_prop.men40.50.F.15.25,
                                      mcar_mar_a_cov.60_prop.men40.50.F.15.25,
                                      mcar_mar_a_cov.65_prop.men40.50.F.15.25,
                                      mcar_mar_a_cov.70_prop.men40.50.F.15.25,
                                      mcar_mar_a_cov.75_prop.men40.50.F.15.25,
                                      mcar_mar_a_cov.80_prop.men40.50.F.15.25,
                                      mcar_mar_a_cov.85_prop.men40.50.F.15.25,
                                      mcar_mar_a_cov.90_prop.men40.50.F.15.25,
                                      mcar_mar_a_cov.95_prop.men40.50.F.15.25)


mcar_mar_a_prop.women15.25.M.40.50 <- c(mcar_mar_a_cov.35_prop.women15.25.M.40.50,
                                        mcar_mar_a_cov.40_prop.women15.25.M.40.50,
                                        mcar_mar_a_cov.45_prop.women15.25.M.40.50,
                                        mcar_mar_a_cov.50_prop.women15.25.M.40.50,
                                        mcar_mar_a_cov.55_prop.women15.25.M.40.50,
                                        mcar_mar_a_cov.60_prop.women15.25.M.40.50,
                                        mcar_mar_a_cov.65_prop.women15.25.M.40.50,
                                        mcar_mar_a_cov.70_prop.women15.25.M.40.50,
                                        mcar_mar_a_cov.75_prop.women15.25.M.40.50,
                                        mcar_mar_a_cov.80_prop.women15.25.M.40.50,
                                        mcar_mar_a_cov.85_prop.women15.25.M.40.50,
                                        mcar_mar_a_cov.90_prop.women15.25.M.40.50,
                                        mcar_mar_a_cov.95_prop.women15.25.M.40.50)


mcar_mar_a_prop.men40.50.F.25.40 <- c(mcar_mar_a_cov.35_prop.men40.50.F.25.40,
                                      mcar_mar_a_cov.40_prop.men40.50.F.25.40,
                                      mcar_mar_a_cov.45_prop.men40.50.F.25.40,
                                      mcar_mar_a_cov.50_prop.men40.50.F.25.40,
                                      mcar_mar_a_cov.55_prop.men40.50.F.25.40,
                                      mcar_mar_a_cov.60_prop.men40.50.F.25.40,
                                      mcar_mar_a_cov.65_prop.men40.50.F.25.40,
                                      mcar_mar_a_cov.70_prop.men40.50.F.25.40,
                                      mcar_mar_a_cov.75_prop.men40.50.F.25.40,
                                      mcar_mar_a_cov.80_prop.men40.50.F.25.40,
                                      mcar_mar_a_cov.85_prop.men40.50.F.25.40,
                                      mcar_mar_a_cov.90_prop.men40.50.F.25.40,
                                      mcar_mar_a_cov.95_prop.men40.50.F.25.40)


mcar_mar_a_prop.women25.40.M.40.50 <- c(mcar_mar_a_cov.35_prop.women25.40.M.40.50,
                                        mcar_mar_a_cov.40_prop.women25.40.M.40.50,
                                        mcar_mar_a_cov.45_prop.women25.40.M.40.50,
                                        mcar_mar_a_cov.50_prop.women25.40.M.40.50,
                                        mcar_mar_a_cov.55_prop.women25.40.M.40.50,
                                        mcar_mar_a_cov.60_prop.women25.40.M.40.50,
                                        mcar_mar_a_cov.65_prop.women25.40.M.40.50,
                                        mcar_mar_a_cov.70_prop.women25.40.M.40.50,
                                        mcar_mar_a_cov.75_prop.women25.40.M.40.50,
                                        mcar_mar_a_cov.80_prop.women25.40.M.40.50,
                                        mcar_mar_a_cov.85_prop.women25.40.M.40.50,
                                        mcar_mar_a_cov.90_prop.women25.40.M.40.50,
                                        mcar_mar_a_cov.95_prop.women25.40.M.40.50)







# MCAR - MAR - b proportions --------------

mcar_mar_b_prop.men15.25.F.15.25 <- c(mcar_mar_b_cov.35_prop.men15.25.F.15.25,
                                      mcar_mar_b_cov.40_prop.men15.25.F.15.25,
                                      mcar_mar_b_cov.45_prop.men15.25.F.15.25,
                                      mcar_mar_b_cov.50_prop.men15.25.F.15.25,
                                      mcar_mar_b_cov.55_prop.men15.25.F.15.25,
                                      mcar_mar_b_cov.60_prop.men15.25.F.15.25,
                                      mcar_mar_b_cov.65_prop.men15.25.F.15.25,
                                      mcar_mar_b_cov.70_prop.men15.25.F.15.25,
                                      mcar_mar_b_cov.75_prop.men15.25.F.15.25,
                                      mcar_mar_b_cov.80_prop.men15.25.F.15.25,
                                      mcar_mar_b_cov.85_prop.men15.25.F.15.25,
                                      mcar_mar_b_cov.90_prop.men15.25.F.15.25,
                                      mcar_mar_b_cov.95_prop.men15.25.F.15.25)



mcar_mar_b_prop.women15.25.M.15.25 <- c(mcar_mar_b_cov.35_prop.women15.25.M.15.25,
                                        mcar_mar_b_cov.40_prop.women15.25.M.15.25,
                                        mcar_mar_b_cov.45_prop.women15.25.M.15.25,
                                        mcar_mar_b_cov.50_prop.women15.25.M.15.25,
                                        mcar_mar_b_cov.55_prop.women15.25.M.15.25,
                                        mcar_mar_b_cov.60_prop.women15.25.M.15.25,
                                        mcar_mar_b_cov.65_prop.women15.25.M.15.25,
                                        mcar_mar_b_cov.70_prop.women15.25.M.15.25,
                                        mcar_mar_b_cov.75_prop.women15.25.M.15.25,
                                        mcar_mar_b_cov.80_prop.women15.25.M.15.25,
                                        mcar_mar_b_cov.85_prop.women15.25.M.15.25,
                                        mcar_mar_b_cov.90_prop.women15.25.M.15.25,
                                        mcar_mar_b_cov.95_prop.women15.25.M.15.25)

mcar_mar_b_prop.men25.40.F.15.25 <- c(mcar_mar_b_cov.35_prop.men25.40.F.15.25,
                                      mcar_mar_b_cov.40_prop.men25.40.F.15.25,
                                      mcar_mar_b_cov.45_prop.men25.40.F.15.25,
                                      mcar_mar_b_cov.50_prop.men25.40.F.15.25,
                                      mcar_mar_b_cov.55_prop.men25.40.F.15.25,
                                      mcar_mar_b_cov.60_prop.men25.40.F.15.25,
                                      mcar_mar_b_cov.65_prop.men25.40.F.15.25,
                                      mcar_mar_b_cov.70_prop.men25.40.F.15.25,
                                      mcar_mar_b_cov.75_prop.men25.40.F.15.25,
                                      mcar_mar_b_cov.80_prop.men25.40.F.15.25,
                                      mcar_mar_b_cov.85_prop.men25.40.F.15.25,
                                      mcar_mar_b_cov.90_prop.men25.40.F.15.25,
                                      mcar_mar_b_cov.95_prop.men25.40.F.15.25)


mcar_mar_b_prop.women15.25.M.25.40 <- c(mcar_mar_b_cov.35_prop.women15.25.M.25.40,
                                        mcar_mar_b_cov.40_prop.women15.25.M.25.40,
                                        mcar_mar_b_cov.45_prop.women15.25.M.25.40,
                                        mcar_mar_b_cov.50_prop.women15.25.M.25.40,
                                        mcar_mar_b_cov.55_prop.women15.25.M.25.40,
                                        mcar_mar_b_cov.60_prop.women15.25.M.25.40,
                                        mcar_mar_b_cov.65_prop.women15.25.M.25.40,
                                        mcar_mar_b_cov.70_prop.women15.25.M.25.40,
                                        mcar_mar_b_cov.75_prop.women15.25.M.25.40,
                                        mcar_mar_b_cov.80_prop.women15.25.M.25.40,
                                        mcar_mar_b_cov.85_prop.women15.25.M.25.40,
                                        mcar_mar_b_cov.90_prop.women15.25.M.25.40,
                                        mcar_mar_b_cov.95_prop.women15.25.M.25.40)

mcar_mar_b_prop.men25.40.F.25.40 <- c(mcar_mar_b_cov.35_prop.men25.40.F.25.40,
                                      mcar_mar_b_cov.40_prop.men25.40.F.25.40,
                                      mcar_mar_b_cov.45_prop.men25.40.F.25.40,
                                      mcar_mar_b_cov.50_prop.men25.40.F.25.40,
                                      mcar_mar_b_cov.55_prop.men25.40.F.25.40,
                                      mcar_mar_b_cov.60_prop.men25.40.F.25.40,
                                      mcar_mar_b_cov.65_prop.men25.40.F.25.40,
                                      mcar_mar_b_cov.70_prop.men25.40.F.25.40,
                                      mcar_mar_b_cov.75_prop.men25.40.F.25.40,
                                      mcar_mar_b_cov.80_prop.men25.40.F.25.40,
                                      mcar_mar_b_cov.85_prop.men25.40.F.25.40,
                                      mcar_mar_b_cov.90_prop.men25.40.F.25.40,
                                      mcar_mar_b_cov.95_prop.men25.40.F.25.40)


mcar_mar_b_prop.women25.40.M.25.40 <- c(mcar_mar_b_cov.35_prop.women25.40.M.25.40,
                                        mcar_mar_b_cov.40_prop.women25.40.M.25.40,
                                        mcar_mar_b_cov.45_prop.women25.40.M.25.40,
                                        mcar_mar_b_cov.50_prop.women25.40.M.25.40,
                                        mcar_mar_b_cov.55_prop.women25.40.M.25.40,
                                        mcar_mar_b_cov.60_prop.women25.40.M.25.40,
                                        mcar_mar_b_cov.65_prop.women25.40.M.25.40,
                                        mcar_mar_b_cov.70_prop.women25.40.M.25.40,
                                        mcar_mar_b_cov.75_prop.women25.40.M.25.40,
                                        mcar_mar_b_cov.80_prop.women25.40.M.25.40,
                                        mcar_mar_b_cov.85_prop.women25.40.M.25.40,
                                        mcar_mar_b_cov.90_prop.women25.40.M.25.40,
                                        mcar_mar_b_cov.95_prop.women25.40.M.25.40)


mcar_mar_b_prop.men40.50.F.15.25 <- c(mcar_mar_b_cov.35_prop.men40.50.F.15.25,
                                      mcar_mar_b_cov.40_prop.men40.50.F.15.25,
                                      mcar_mar_b_cov.45_prop.men40.50.F.15.25,
                                      mcar_mar_b_cov.50_prop.men40.50.F.15.25,
                                      mcar_mar_b_cov.55_prop.men40.50.F.15.25,
                                      mcar_mar_b_cov.60_prop.men40.50.F.15.25,
                                      mcar_mar_b_cov.65_prop.men40.50.F.15.25,
                                      mcar_mar_b_cov.70_prop.men40.50.F.15.25,
                                      mcar_mar_b_cov.75_prop.men40.50.F.15.25,
                                      mcar_mar_b_cov.80_prop.men40.50.F.15.25,
                                      mcar_mar_b_cov.85_prop.men40.50.F.15.25,
                                      mcar_mar_b_cov.90_prop.men40.50.F.15.25,
                                      mcar_mar_b_cov.95_prop.men40.50.F.15.25)


mcar_mar_b_prop.women15.25.M.40.50 <- c(mcar_mar_b_cov.35_prop.women15.25.M.40.50,
                                        mcar_mar_b_cov.40_prop.women15.25.M.40.50,
                                        mcar_mar_b_cov.45_prop.women15.25.M.40.50,
                                        mcar_mar_b_cov.50_prop.women15.25.M.40.50,
                                        mcar_mar_b_cov.55_prop.women15.25.M.40.50,
                                        mcar_mar_b_cov.60_prop.women15.25.M.40.50,
                                        mcar_mar_b_cov.65_prop.women15.25.M.40.50,
                                        mcar_mar_b_cov.70_prop.women15.25.M.40.50,
                                        mcar_mar_b_cov.75_prop.women15.25.M.40.50,
                                        mcar_mar_b_cov.80_prop.women15.25.M.40.50,
                                        mcar_mar_b_cov.85_prop.women15.25.M.40.50,
                                        mcar_mar_b_cov.90_prop.women15.25.M.40.50,
                                        mcar_mar_b_cov.95_prop.women15.25.M.40.50)


mcar_mar_b_prop.men40.50.F.25.40 <- c(mcar_mar_b_cov.35_prop.men40.50.F.25.40,
                                      mcar_mar_b_cov.40_prop.men40.50.F.25.40,
                                      mcar_mar_b_cov.45_prop.men40.50.F.25.40,
                                      mcar_mar_b_cov.50_prop.men40.50.F.25.40,
                                      mcar_mar_b_cov.55_prop.men40.50.F.25.40,
                                      mcar_mar_b_cov.60_prop.men40.50.F.25.40,
                                      mcar_mar_b_cov.65_prop.men40.50.F.25.40,
                                      mcar_mar_b_cov.70_prop.men40.50.F.25.40,
                                      mcar_mar_b_cov.75_prop.men40.50.F.25.40,
                                      mcar_mar_b_cov.80_prop.men40.50.F.25.40,
                                      mcar_mar_b_cov.85_prop.men40.50.F.25.40,
                                      mcar_mar_b_cov.90_prop.men40.50.F.25.40,
                                      mcar_mar_b_cov.95_prop.men40.50.F.25.40)


mcar_mar_b_prop.women25.40.M.40.50 <- c(mcar_mar_b_cov.35_prop.women25.40.M.40.50,
                                        mcar_mar_b_cov.40_prop.women25.40.M.40.50,
                                        mcar_mar_b_cov.45_prop.women25.40.M.40.50,
                                        mcar_mar_b_cov.50_prop.women25.40.M.40.50,
                                        mcar_mar_b_cov.55_prop.women25.40.M.40.50,
                                        mcar_mar_b_cov.60_prop.women25.40.M.40.50,
                                        mcar_mar_b_cov.65_prop.women25.40.M.40.50,
                                        mcar_mar_b_cov.70_prop.women25.40.M.40.50,
                                        mcar_mar_b_cov.75_prop.women25.40.M.40.50,
                                        mcar_mar_b_cov.80_prop.women25.40.M.40.50,
                                        mcar_mar_b_cov.85_prop.women25.40.M.40.50,
                                        mcar_mar_b_cov.90_prop.women25.40.M.40.50,
                                        mcar_mar_b_cov.95_prop.women25.40.M.40.50)









# MCAR - MAR - c proportions --------------

mcar_mar_c_prop.men15.25.F.15.25 <- c(mcar_mar_c_cov.35_prop.men15.25.F.15.25,
                                      mcar_mar_c_cov.40_prop.men15.25.F.15.25,
                                      mcar_mar_c_cov.45_prop.men15.25.F.15.25,
                                      mcar_mar_c_cov.50_prop.men15.25.F.15.25,
                                      mcar_mar_c_cov.55_prop.men15.25.F.15.25,
                                      mcar_mar_c_cov.60_prop.men15.25.F.15.25,
                                      mcar_mar_c_cov.65_prop.men15.25.F.15.25,
                                      mcar_mar_c_cov.70_prop.men15.25.F.15.25,
                                      mcar_mar_c_cov.75_prop.men15.25.F.15.25,
                                      mcar_mar_c_cov.80_prop.men15.25.F.15.25,
                                      mcar_mar_c_cov.85_prop.men15.25.F.15.25,
                                      mcar_mar_c_cov.90_prop.men15.25.F.15.25,
                                      mcar_mar_c_cov.95_prop.men15.25.F.15.25)



mcar_mar_c_prop.women15.25.M.15.25 <- c(mcar_mar_c_cov.35_prop.women15.25.M.15.25,
                                        mcar_mar_c_cov.40_prop.women15.25.M.15.25,
                                        mcar_mar_c_cov.45_prop.women15.25.M.15.25,
                                        mcar_mar_c_cov.50_prop.women15.25.M.15.25,
                                        mcar_mar_c_cov.55_prop.women15.25.M.15.25,
                                        mcar_mar_c_cov.60_prop.women15.25.M.15.25,
                                        mcar_mar_c_cov.65_prop.women15.25.M.15.25,
                                        mcar_mar_c_cov.70_prop.women15.25.M.15.25,
                                        mcar_mar_c_cov.75_prop.women15.25.M.15.25,
                                        mcar_mar_c_cov.80_prop.women15.25.M.15.25,
                                        mcar_mar_c_cov.85_prop.women15.25.M.15.25,
                                        mcar_mar_c_cov.90_prop.women15.25.M.15.25,
                                        mcar_mar_c_cov.95_prop.women15.25.M.15.25)

mcar_mar_c_prop.men25.40.F.15.25 <- c(mcar_mar_c_cov.35_prop.men25.40.F.15.25,
                                      mcar_mar_c_cov.40_prop.men25.40.F.15.25,
                                      mcar_mar_c_cov.45_prop.men25.40.F.15.25,
                                      mcar_mar_c_cov.50_prop.men25.40.F.15.25,
                                      mcar_mar_c_cov.55_prop.men25.40.F.15.25,
                                      mcar_mar_c_cov.60_prop.men25.40.F.15.25,
                                      mcar_mar_c_cov.65_prop.men25.40.F.15.25,
                                      mcar_mar_c_cov.70_prop.men25.40.F.15.25,
                                      mcar_mar_c_cov.75_prop.men25.40.F.15.25,
                                      mcar_mar_c_cov.80_prop.men25.40.F.15.25,
                                      mcar_mar_c_cov.85_prop.men25.40.F.15.25,
                                      mcar_mar_c_cov.90_prop.men25.40.F.15.25,
                                      mcar_mar_c_cov.95_prop.men25.40.F.15.25)


mcar_mar_c_prop.women15.25.M.25.40 <- c(mcar_mar_c_cov.35_prop.women15.25.M.25.40,
                                        mcar_mar_c_cov.40_prop.women15.25.M.25.40,
                                        mcar_mar_c_cov.45_prop.women15.25.M.25.40,
                                        mcar_mar_c_cov.50_prop.women15.25.M.25.40,
                                        mcar_mar_c_cov.55_prop.women15.25.M.25.40,
                                        mcar_mar_c_cov.60_prop.women15.25.M.25.40,
                                        mcar_mar_c_cov.65_prop.women15.25.M.25.40,
                                        mcar_mar_c_cov.70_prop.women15.25.M.25.40,
                                        mcar_mar_c_cov.75_prop.women15.25.M.25.40,
                                        mcar_mar_c_cov.80_prop.women15.25.M.25.40,
                                        mcar_mar_c_cov.85_prop.women15.25.M.25.40,
                                        mcar_mar_c_cov.90_prop.women15.25.M.25.40,
                                        mcar_mar_c_cov.95_prop.women15.25.M.25.40)

mcar_mar_c_prop.men25.40.F.25.40 <- c(mcar_mar_c_cov.35_prop.men25.40.F.25.40,
                                      mcar_mar_c_cov.40_prop.men25.40.F.25.40,
                                      mcar_mar_c_cov.45_prop.men25.40.F.25.40,
                                      mcar_mar_c_cov.50_prop.men25.40.F.25.40,
                                      mcar_mar_c_cov.55_prop.men25.40.F.25.40,
                                      mcar_mar_c_cov.60_prop.men25.40.F.25.40,
                                      mcar_mar_c_cov.65_prop.men25.40.F.25.40,
                                      mcar_mar_c_cov.70_prop.men25.40.F.25.40,
                                      mcar_mar_c_cov.75_prop.men25.40.F.25.40,
                                      mcar_mar_c_cov.80_prop.men25.40.F.25.40,
                                      mcar_mar_c_cov.85_prop.men25.40.F.25.40,
                                      mcar_mar_c_cov.90_prop.men25.40.F.25.40,
                                      mcar_mar_c_cov.95_prop.men25.40.F.25.40)


mcar_mar_c_prop.women25.40.M.25.40 <- c(mcar_mar_c_cov.35_prop.women25.40.M.25.40,
                                        mcar_mar_c_cov.40_prop.women25.40.M.25.40,
                                        mcar_mar_c_cov.45_prop.women25.40.M.25.40,
                                        mcar_mar_c_cov.50_prop.women25.40.M.25.40,
                                        mcar_mar_c_cov.55_prop.women25.40.M.25.40,
                                        mcar_mar_c_cov.60_prop.women25.40.M.25.40,
                                        mcar_mar_c_cov.65_prop.women25.40.M.25.40,
                                        mcar_mar_c_cov.70_prop.women25.40.M.25.40,
                                        mcar_mar_c_cov.75_prop.women25.40.M.25.40,
                                        mcar_mar_c_cov.80_prop.women25.40.M.25.40,
                                        mcar_mar_c_cov.85_prop.women25.40.M.25.40,
                                        mcar_mar_c_cov.90_prop.women25.40.M.25.40,
                                        mcar_mar_c_cov.95_prop.women25.40.M.25.40)


mcar_mar_c_prop.men40.50.F.15.25 <- c(mcar_mar_c_cov.35_prop.men40.50.F.15.25,
                                      mcar_mar_c_cov.40_prop.men40.50.F.15.25,
                                      mcar_mar_c_cov.45_prop.men40.50.F.15.25,
                                      mcar_mar_c_cov.50_prop.men40.50.F.15.25,
                                      mcar_mar_c_cov.55_prop.men40.50.F.15.25,
                                      mcar_mar_c_cov.60_prop.men40.50.F.15.25,
                                      mcar_mar_c_cov.65_prop.men40.50.F.15.25,
                                      mcar_mar_c_cov.70_prop.men40.50.F.15.25,
                                      mcar_mar_c_cov.75_prop.men40.50.F.15.25,
                                      mcar_mar_c_cov.80_prop.men40.50.F.15.25,
                                      mcar_mar_c_cov.85_prop.men40.50.F.15.25,
                                      mcar_mar_c_cov.90_prop.men40.50.F.15.25,
                                      mcar_mar_c_cov.95_prop.men40.50.F.15.25)


mcar_mar_c_prop.women15.25.M.40.50 <- c(mcar_mar_c_cov.35_prop.women15.25.M.40.50,
                                        mcar_mar_c_cov.40_prop.women15.25.M.40.50,
                                        mcar_mar_c_cov.45_prop.women15.25.M.40.50,
                                        mcar_mar_c_cov.50_prop.women15.25.M.40.50,
                                        mcar_mar_c_cov.55_prop.women15.25.M.40.50,
                                        mcar_mar_c_cov.60_prop.women15.25.M.40.50,
                                        mcar_mar_c_cov.65_prop.women15.25.M.40.50,
                                        mcar_mar_c_cov.70_prop.women15.25.M.40.50,
                                        mcar_mar_c_cov.75_prop.women15.25.M.40.50,
                                        mcar_mar_c_cov.80_prop.women15.25.M.40.50,
                                        mcar_mar_c_cov.85_prop.women15.25.M.40.50,
                                        mcar_mar_c_cov.90_prop.women15.25.M.40.50,
                                        mcar_mar_c_cov.95_prop.women15.25.M.40.50)


mcar_mar_c_prop.men40.50.F.25.40 <- c(mcar_mar_c_cov.35_prop.men40.50.F.25.40,
                                      mcar_mar_c_cov.40_prop.men40.50.F.25.40,
                                      mcar_mar_c_cov.45_prop.men40.50.F.25.40,
                                      mcar_mar_c_cov.50_prop.men40.50.F.25.40,
                                      mcar_mar_c_cov.55_prop.men40.50.F.25.40,
                                      mcar_mar_c_cov.60_prop.men40.50.F.25.40,
                                      mcar_mar_c_cov.65_prop.men40.50.F.25.40,
                                      mcar_mar_c_cov.70_prop.men40.50.F.25.40,
                                      mcar_mar_c_cov.75_prop.men40.50.F.25.40,
                                      mcar_mar_c_cov.80_prop.men40.50.F.25.40,
                                      mcar_mar_c_cov.85_prop.men40.50.F.25.40,
                                      mcar_mar_c_cov.90_prop.men40.50.F.25.40,
                                      mcar_mar_c_cov.95_prop.men40.50.F.25.40)


mcar_mar_c_prop.women25.40.M.40.50 <- c(mcar_mar_c_cov.35_prop.women25.40.M.40.50,
                                        mcar_mar_c_cov.40_prop.women25.40.M.40.50,
                                        mcar_mar_c_cov.45_prop.women25.40.M.40.50,
                                        mcar_mar_c_cov.50_prop.women25.40.M.40.50,
                                        mcar_mar_c_cov.55_prop.women25.40.M.40.50,
                                        mcar_mar_c_cov.60_prop.women25.40.M.40.50,
                                        mcar_mar_c_cov.65_prop.women25.40.M.40.50,
                                        mcar_mar_c_cov.70_prop.women25.40.M.40.50,
                                        mcar_mar_c_cov.75_prop.women25.40.M.40.50,
                                        mcar_mar_c_cov.80_prop.women25.40.M.40.50,
                                        mcar_mar_c_cov.85_prop.women25.40.M.40.50,
                                        mcar_mar_c_cov.90_prop.women25.40.M.40.50,
                                        mcar_mar_c_cov.95_prop.women25.40.M.40.50)



# MCAR - MAR - a AD -------------- 




# mean

mcar_mar_a_mean.men.cl.15.25 <- c(mcar_mar_a_cov.35_mean.men.cl.15.25,
                                  mcar_mar_a_cov.40_mean.men.cl.15.25,
                                  mcar_mar_a_cov.45_mean.men.cl.15.25,
                                  mcar_mar_a_cov.50_mean.men.cl.15.25,
                                  mcar_mar_a_cov.55_mean.men.cl.15.25,
                                  mcar_mar_a_cov.60_mean.men.cl.15.25,
                                  mcar_mar_a_cov.65_mean.men.cl.15.25,
                                  mcar_mar_a_cov.70_mean.men.cl.15.25,
                                  mcar_mar_a_cov.75_mean.men.cl.15.25,
                                  mcar_mar_a_cov.80_mean.men.cl.15.25,
                                  mcar_mar_a_cov.85_mean.men.cl.15.25,
                                  mcar_mar_a_cov.90_mean.men.cl.15.25,
                                  mcar_mar_a_cov.95_mean.men.cl.15.25)



mcar_mar_a_mean.women.cl.15.25 <- c(mcar_mar_a_cov.35_mean.women.cl.15.25,
                                    mcar_mar_a_cov.40_mean.women.cl.15.25,
                                    mcar_mar_a_cov.45_mean.women.cl.15.25,
                                    mcar_mar_a_cov.50_mean.women.cl.15.25,
                                    mcar_mar_a_cov.55_mean.women.cl.15.25,
                                    mcar_mar_a_cov.60_mean.women.cl.15.25,
                                    mcar_mar_a_cov.65_mean.women.cl.15.25,
                                    mcar_mar_a_cov.70_mean.women.cl.15.25,
                                    mcar_mar_a_cov.75_mean.women.cl.15.25,
                                    mcar_mar_a_cov.80_mean.women.cl.15.25,
                                    mcar_mar_a_cov.85_mean.women.cl.15.25,
                                    mcar_mar_a_cov.90_mean.women.cl.15.25,
                                    mcar_mar_a_cov.95_mean.women.cl.15.25)



mcar_mar_a_mean.men.cl.25.40 <- c(mcar_mar_a_cov.35_mean.men.cl.25.40,
                                  mcar_mar_a_cov.40_mean.men.cl.25.40,
                                  mcar_mar_a_cov.45_mean.men.cl.25.40,
                                  mcar_mar_a_cov.50_mean.men.cl.25.40,
                                  mcar_mar_a_cov.55_mean.men.cl.25.40,
                                  mcar_mar_a_cov.60_mean.men.cl.25.40,
                                  mcar_mar_a_cov.65_mean.men.cl.25.40,
                                  mcar_mar_a_cov.70_mean.men.cl.25.40,
                                  mcar_mar_a_cov.75_mean.men.cl.25.40,
                                  mcar_mar_a_cov.80_mean.men.cl.25.40,
                                  mcar_mar_a_cov.85_mean.men.cl.25.40,
                                  mcar_mar_a_cov.90_mean.men.cl.25.40,
                                  mcar_mar_a_cov.95_mean.men.cl.25.40)

mcar_mar_a_mean.women.cl.25.40 <- c(mcar_mar_a_cov.35_mean.women.cl.25.40,
                                    mcar_mar_a_cov.40_mean.women.cl.25.40,
                                    mcar_mar_a_cov.45_mean.women.cl.25.40,
                                    mcar_mar_a_cov.50_mean.women.cl.25.40,
                                    mcar_mar_a_cov.55_mean.women.cl.25.40,
                                    mcar_mar_a_cov.60_mean.women.cl.25.40,
                                    mcar_mar_a_cov.65_mean.women.cl.25.40,
                                    mcar_mar_a_cov.70_mean.women.cl.25.40,
                                    mcar_mar_a_cov.75_mean.women.cl.25.40,
                                    mcar_mar_a_cov.80_mean.women.cl.25.40,
                                    mcar_mar_a_cov.85_mean.women.cl.25.40,
                                    mcar_mar_a_cov.90_mean.women.cl.25.40,
                                    mcar_mar_a_cov.95_mean.women.cl.25.40)



mcar_mar_a_mean.men.cl.40.50 <- c(mcar_mar_a_cov.35_mean.men.cl.40.50,
                                  mcar_mar_a_cov.40_mean.men.cl.40.50,
                                  mcar_mar_a_cov.45_mean.men.cl.40.50,
                                  mcar_mar_a_cov.50_mean.men.cl.40.50,
                                  mcar_mar_a_cov.55_mean.men.cl.40.50,
                                  mcar_mar_a_cov.60_mean.men.cl.40.50,
                                  mcar_mar_a_cov.65_mean.men.cl.40.50,
                                  mcar_mar_a_cov.70_mean.men.cl.40.50,
                                  mcar_mar_a_cov.75_mean.men.cl.40.50,
                                  mcar_mar_a_cov.80_mean.men.cl.40.50,
                                  mcar_mar_a_cov.85_mean.men.cl.40.50,
                                  mcar_mar_a_cov.90_mean.men.cl.40.50,
                                  mcar_mar_a_cov.95_mean.men.cl.40.50)



mcar_mar_a_mean.women.cl.40.50 <- c(mcar_mar_a_cov.35_mean.women.cl.40.50,
                                    mcar_mar_a_cov.40_mean.women.cl.40.50,
                                    mcar_mar_a_cov.45_mean.women.cl.40.50,
                                    mcar_mar_a_cov.50_mean.women.cl.40.50,
                                    mcar_mar_a_cov.55_mean.women.cl.40.50,
                                    mcar_mar_a_cov.60_mean.women.cl.40.50,
                                    mcar_mar_a_cov.65_mean.women.cl.40.50,
                                    mcar_mar_a_cov.70_mean.women.cl.40.50,
                                    mcar_mar_a_cov.75_mean.women.cl.40.50,
                                    mcar_mar_a_cov.80_mean.women.cl.40.50,
                                    mcar_mar_a_cov.85_mean.women.cl.40.50,
                                    mcar_mar_a_cov.90_mean.women.cl.40.50,
                                    mcar_mar_a_cov.95_mean.women.cl.40.50)


# median


mcar_mar_a_med.men.cl.15.25 <- c(mcar_mar_a_cov.35_med.men.cl.15.25,
                                 mcar_mar_a_cov.40_med.men.cl.15.25,
                                 mcar_mar_a_cov.45_med.men.cl.15.25,
                                 mcar_mar_a_cov.50_med.men.cl.15.25,
                                 mcar_mar_a_cov.55_med.men.cl.15.25,
                                 mcar_mar_a_cov.60_med.men.cl.15.25,
                                 mcar_mar_a_cov.65_med.men.cl.15.25,
                                 mcar_mar_a_cov.70_med.men.cl.15.25,
                                 mcar_mar_a_cov.75_med.men.cl.15.25,
                                 mcar_mar_a_cov.80_med.men.cl.15.25,
                                 mcar_mar_a_cov.85_med.men.cl.15.25,
                                 mcar_mar_a_cov.90_med.men.cl.15.25,
                                 mcar_mar_a_cov.95_med.men.cl.15.25)



mcar_mar_a_med.women.cl.15.25 <- c(mcar_mar_a_cov.35_med.women.cl.15.25,
                                   mcar_mar_a_cov.40_med.women.cl.15.25,
                                   mcar_mar_a_cov.45_med.women.cl.15.25,
                                   mcar_mar_a_cov.50_med.women.cl.15.25,
                                   mcar_mar_a_cov.55_med.women.cl.15.25,
                                   mcar_mar_a_cov.60_med.women.cl.15.25,
                                   mcar_mar_a_cov.65_med.women.cl.15.25,
                                   mcar_mar_a_cov.70_med.women.cl.15.25,
                                   mcar_mar_a_cov.75_med.women.cl.15.25,
                                   mcar_mar_a_cov.80_med.women.cl.15.25,
                                   mcar_mar_a_cov.85_med.women.cl.15.25,
                                   mcar_mar_a_cov.90_med.women.cl.15.25,
                                   mcar_mar_a_cov.95_med.women.cl.15.25)



mcar_mar_a_med.men.cl.25.40 <- c(mcar_mar_a_cov.35_med.men.cl.25.40,
                                 mcar_mar_a_cov.40_med.men.cl.25.40,
                                 mcar_mar_a_cov.45_med.men.cl.25.40,
                                 mcar_mar_a_cov.50_med.men.cl.25.40,
                                 mcar_mar_a_cov.55_med.men.cl.25.40,
                                 mcar_mar_a_cov.60_med.men.cl.25.40,
                                 mcar_mar_a_cov.65_med.men.cl.25.40,
                                 mcar_mar_a_cov.70_med.men.cl.25.40,
                                 mcar_mar_a_cov.75_med.men.cl.25.40,
                                 mcar_mar_a_cov.80_med.men.cl.25.40,
                                 mcar_mar_a_cov.85_med.men.cl.25.40,
                                 mcar_mar_a_cov.90_med.men.cl.25.40,
                                 mcar_mar_a_cov.95_med.men.cl.25.40)


mcar_mar_a_med.women.cl.25.40 <- c(mcar_mar_a_cov.35_med.women.cl.25.40,
                                   mcar_mar_a_cov.40_med.women.cl.25.40,
                                   mcar_mar_a_cov.45_med.women.cl.25.40,
                                   mcar_mar_a_cov.50_med.women.cl.25.40,
                                   mcar_mar_a_cov.55_med.women.cl.25.40,
                                   mcar_mar_a_cov.60_med.women.cl.25.40,
                                   mcar_mar_a_cov.65_med.women.cl.25.40,
                                   mcar_mar_a_cov.70_med.women.cl.25.40,
                                   mcar_mar_a_cov.75_med.women.cl.25.40,
                                   mcar_mar_a_cov.80_med.women.cl.25.40,
                                   mcar_mar_a_cov.85_med.women.cl.25.40,
                                   mcar_mar_a_cov.90_med.women.cl.25.40,
                                   mcar_mar_a_cov.95_med.women.cl.25.40)




mcar_mar_a_med.men.cl.40.50 <- c(mcar_mar_a_cov.35_med.men.cl.40.50,
                                 mcar_mar_a_cov.40_med.men.cl.40.50,
                                 mcar_mar_a_cov.45_med.men.cl.40.50,
                                 mcar_mar_a_cov.50_med.men.cl.40.50,
                                 mcar_mar_a_cov.55_med.men.cl.40.50,
                                 mcar_mar_a_cov.60_med.men.cl.40.50,
                                 mcar_mar_a_cov.65_med.men.cl.40.50,
                                 mcar_mar_a_cov.70_med.men.cl.40.50,
                                 mcar_mar_a_cov.75_med.men.cl.40.50,
                                 mcar_mar_a_cov.80_med.men.cl.40.50,
                                 mcar_mar_a_cov.85_med.men.cl.40.50,
                                 mcar_mar_a_cov.90_med.men.cl.40.50,
                                 mcar_mar_a_cov.95_med.men.cl.40.50)



mcar_mar_a_med.women.cl.40.50 <- c(mcar_mar_a_cov.35_med.women.cl.40.50,
                                   mcar_mar_a_cov.40_med.women.cl.40.50,
                                   mcar_mar_a_cov.45_med.women.cl.40.50,
                                   mcar_mar_a_cov.50_med.women.cl.40.50,
                                   mcar_mar_a_cov.55_med.women.cl.40.50,
                                   mcar_mar_a_cov.60_med.women.cl.40.50,
                                   mcar_mar_a_cov.65_med.women.cl.40.50,
                                   mcar_mar_a_cov.70_med.women.cl.40.50,
                                   mcar_mar_a_cov.75_med.women.cl.40.50,
                                   mcar_mar_a_cov.80_med.women.cl.40.50,
                                   mcar_mar_a_cov.85_med.women.cl.40.50,
                                   mcar_mar_a_cov.90_med.women.cl.40.50,
                                   mcar_mar_a_cov.95_med.women.cl.40.50)



# standard deviation

mcar_mar_a_sd.men.cl.15.25 <- c(mcar_mar_a_cov.35_sd.men.cl.15.25,
                                mcar_mar_a_cov.40_sd.men.cl.15.25,
                                mcar_mar_a_cov.45_sd.men.cl.15.25,
                                mcar_mar_a_cov.50_sd.men.cl.15.25,
                                mcar_mar_a_cov.55_sd.men.cl.15.25,
                                mcar_mar_a_cov.60_sd.men.cl.15.25,
                                mcar_mar_a_cov.65_sd.men.cl.15.25,
                                mcar_mar_a_cov.70_sd.men.cl.15.25,
                                mcar_mar_a_cov.75_sd.men.cl.15.25,
                                mcar_mar_a_cov.80_sd.men.cl.15.25,
                                mcar_mar_a_cov.85_sd.men.cl.15.25,
                                mcar_mar_a_cov.90_sd.men.cl.15.25,
                                mcar_mar_a_cov.95_sd.men.cl.15.25)



mcar_mar_a_sd.women.cl.15.25 <- c(mcar_mar_a_cov.35_sd.women.cl.15.25,
                                  mcar_mar_a_cov.40_sd.women.cl.15.25,
                                  mcar_mar_a_cov.45_sd.women.cl.15.25,
                                  mcar_mar_a_cov.50_sd.women.cl.15.25,
                                  mcar_mar_a_cov.55_sd.women.cl.15.25,
                                  mcar_mar_a_cov.60_sd.women.cl.15.25,
                                  mcar_mar_a_cov.65_sd.women.cl.15.25,
                                  mcar_mar_a_cov.70_sd.women.cl.15.25,
                                  mcar_mar_a_cov.75_sd.women.cl.15.25,
                                  mcar_mar_a_cov.80_sd.women.cl.15.25,
                                  mcar_mar_a_cov.85_sd.women.cl.15.25,
                                  mcar_mar_a_cov.90_sd.women.cl.15.25,
                                  mcar_mar_a_cov.95_sd.women.cl.15.25)



mcar_mar_a_sd.men.cl.25.40 <- c(mcar_mar_a_cov.35_sd.men.cl.25.40,
                                mcar_mar_a_cov.40_sd.men.cl.25.40,
                                mcar_mar_a_cov.45_sd.men.cl.25.40,
                                mcar_mar_a_cov.50_sd.men.cl.25.40,
                                mcar_mar_a_cov.55_sd.men.cl.25.40,
                                mcar_mar_a_cov.60_sd.men.cl.25.40,
                                mcar_mar_a_cov.65_sd.men.cl.25.40,
                                mcar_mar_a_cov.70_sd.men.cl.25.40,
                                mcar_mar_a_cov.75_sd.men.cl.25.40,
                                mcar_mar_a_cov.80_sd.men.cl.25.40,
                                mcar_mar_a_cov.85_sd.men.cl.25.40,
                                mcar_mar_a_cov.90_sd.men.cl.25.40,
                                mcar_mar_a_cov.95_sd.men.cl.25.40)



mcar_mar_a_sd.women.cl.25.40 <- c(mcar_mar_a_cov.35_sd.women.cl.25.40,
                                  mcar_mar_a_cov.40_sd.women.cl.25.40,
                                  mcar_mar_a_cov.45_sd.women.cl.25.40,
                                  mcar_mar_a_cov.50_sd.women.cl.25.40,
                                  mcar_mar_a_cov.55_sd.women.cl.25.40,
                                  mcar_mar_a_cov.60_sd.women.cl.25.40,
                                  mcar_mar_a_cov.65_sd.women.cl.25.40,
                                  mcar_mar_a_cov.70_sd.women.cl.25.40,
                                  mcar_mar_a_cov.75_sd.women.cl.25.40,
                                  mcar_mar_a_cov.80_sd.women.cl.25.40,
                                  mcar_mar_a_cov.85_sd.women.cl.25.40,
                                  mcar_mar_a_cov.90_sd.women.cl.25.40,
                                  mcar_mar_a_cov.95_sd.women.cl.25.40)




mcar_mar_a_sd.men.cl.40.50 <- c(mcar_mar_a_cov.35_sd.men.cl.40.50,
                                mcar_mar_a_cov.40_sd.men.cl.40.50,
                                mcar_mar_a_cov.45_sd.men.cl.40.50,
                                mcar_mar_a_cov.50_sd.men.cl.40.50,
                                mcar_mar_a_cov.55_sd.men.cl.40.50,
                                mcar_mar_a_cov.60_sd.men.cl.40.50,
                                mcar_mar_a_cov.65_sd.men.cl.40.50,
                                mcar_mar_a_cov.70_sd.men.cl.40.50,
                                mcar_mar_a_cov.75_sd.men.cl.40.50,
                                mcar_mar_a_cov.80_sd.men.cl.40.50,
                                mcar_mar_a_cov.85_sd.men.cl.40.50,
                                mcar_mar_a_cov.90_sd.men.cl.40.50,
                                mcar_mar_a_cov.95_sd.men.cl.40.50)



mcar_mar_a_sd.women.cl.40.50 <- c(mcar_mar_a_cov.35_sd.women.cl.40.50,
                                  mcar_mar_a_cov.40_sd.women.cl.40.50,
                                  mcar_mar_a_cov.45_sd.women.cl.40.50,
                                  mcar_mar_a_cov.50_sd.women.cl.40.50,
                                  mcar_mar_a_cov.55_sd.women.cl.40.50,
                                  mcar_mar_a_cov.60_sd.women.cl.40.50,
                                  mcar_mar_a_cov.65_sd.women.cl.40.50,
                                  mcar_mar_a_cov.70_sd.women.cl.40.50,
                                  mcar_mar_a_cov.75_sd.women.cl.40.50,
                                  mcar_mar_a_cov.80_sd.women.cl.40.50,
                                  mcar_mar_a_cov.85_sd.women.cl.40.50,
                                  mcar_mar_a_cov.90_sd.women.cl.40.50,
                                  mcar_mar_a_cov.95_sd.women.cl.40.50)





# MCAR - MAR - b AD -------------- 





# mean

mcar_mar_b_mean.men.cl.15.25 <- c(mcar_mar_b_cov.35_mean.men.cl.15.25,
                                  mcar_mar_b_cov.40_mean.men.cl.15.25,
                                  mcar_mar_b_cov.45_mean.men.cl.15.25,
                                  mcar_mar_b_cov.50_mean.men.cl.15.25,
                                  mcar_mar_b_cov.55_mean.men.cl.15.25,
                                  mcar_mar_b_cov.60_mean.men.cl.15.25,
                                  mcar_mar_b_cov.65_mean.men.cl.15.25,
                                  mcar_mar_b_cov.70_mean.men.cl.15.25,
                                  mcar_mar_b_cov.75_mean.men.cl.15.25,
                                  mcar_mar_b_cov.80_mean.men.cl.15.25,
                                  mcar_mar_b_cov.85_mean.men.cl.15.25,
                                  mcar_mar_b_cov.90_mean.men.cl.15.25,
                                  mcar_mar_b_cov.95_mean.men.cl.15.25)



mcar_mar_b_mean.women.cl.15.25 <- c(mcar_mar_b_cov.35_mean.women.cl.15.25,
                                    mcar_mar_b_cov.40_mean.women.cl.15.25,
                                    mcar_mar_b_cov.45_mean.women.cl.15.25,
                                    mcar_mar_b_cov.50_mean.women.cl.15.25,
                                    mcar_mar_b_cov.55_mean.women.cl.15.25,
                                    mcar_mar_b_cov.60_mean.women.cl.15.25,
                                    mcar_mar_b_cov.65_mean.women.cl.15.25,
                                    mcar_mar_b_cov.70_mean.women.cl.15.25,
                                    mcar_mar_b_cov.75_mean.women.cl.15.25,
                                    mcar_mar_b_cov.80_mean.women.cl.15.25,
                                    mcar_mar_b_cov.85_mean.women.cl.15.25,
                                    mcar_mar_b_cov.90_mean.women.cl.15.25,
                                    mcar_mar_b_cov.95_mean.women.cl.15.25)



mcar_mar_b_mean.men.cl.25.40 <- c(mcar_mar_b_cov.35_mean.men.cl.25.40,
                                  mcar_mar_b_cov.40_mean.men.cl.25.40,
                                  mcar_mar_b_cov.45_mean.men.cl.25.40,
                                  mcar_mar_b_cov.50_mean.men.cl.25.40,
                                  mcar_mar_b_cov.55_mean.men.cl.25.40,
                                  mcar_mar_b_cov.60_mean.men.cl.25.40,
                                  mcar_mar_b_cov.65_mean.men.cl.25.40,
                                  mcar_mar_b_cov.70_mean.men.cl.25.40,
                                  mcar_mar_b_cov.75_mean.men.cl.25.40,
                                  mcar_mar_b_cov.80_mean.men.cl.25.40,
                                  mcar_mar_b_cov.85_mean.men.cl.25.40,
                                  mcar_mar_b_cov.90_mean.men.cl.25.40,
                                  mcar_mar_b_cov.95_mean.men.cl.25.40)


mcar_mar_b_mean.women.cl.25.40 <- c(mcar_mar_b_cov.35_mean.women.cl.25.40,
                                    mcar_mar_b_cov.40_mean.women.cl.25.40,
                                    mcar_mar_b_cov.45_mean.women.cl.25.40,
                                    mcar_mar_b_cov.50_mean.women.cl.25.40,
                                    mcar_mar_b_cov.55_mean.women.cl.25.40,
                                    mcar_mar_b_cov.60_mean.women.cl.25.40,
                                    mcar_mar_b_cov.65_mean.women.cl.25.40,
                                    mcar_mar_b_cov.70_mean.women.cl.25.40,
                                    mcar_mar_b_cov.75_mean.women.cl.25.40,
                                    mcar_mar_b_cov.80_mean.women.cl.25.40,
                                    mcar_mar_b_cov.85_mean.women.cl.25.40,
                                    mcar_mar_b_cov.90_mean.women.cl.25.40,
                                    mcar_mar_b_cov.95_mean.women.cl.25.40)


mcar_mar_b_mean.men.cl.40.50 <- c(mcar_mar_b_cov.35_mean.men.cl.40.50,
                                  mcar_mar_b_cov.40_mean.men.cl.40.50,
                                  mcar_mar_b_cov.45_mean.men.cl.40.50,
                                  mcar_mar_b_cov.50_mean.men.cl.40.50,
                                  mcar_mar_b_cov.55_mean.men.cl.40.50,
                                  mcar_mar_b_cov.60_mean.men.cl.40.50,
                                  mcar_mar_b_cov.65_mean.men.cl.40.50,
                                  mcar_mar_b_cov.70_mean.men.cl.40.50,
                                  mcar_mar_b_cov.75_mean.men.cl.40.50,
                                  mcar_mar_b_cov.80_mean.men.cl.40.50,
                                  mcar_mar_b_cov.85_mean.men.cl.40.50,
                                  mcar_mar_b_cov.90_mean.men.cl.40.50,
                                  mcar_mar_b_cov.95_mean.men.cl.40.50)



mcar_mar_b_mean.women.cl.40.50 <- c(mcar_mar_b_cov.35_mean.women.cl.40.50,
                                    mcar_mar_b_cov.40_mean.women.cl.40.50,
                                    mcar_mar_b_cov.45_mean.women.cl.40.50,
                                    mcar_mar_b_cov.50_mean.women.cl.40.50,
                                    mcar_mar_b_cov.55_mean.women.cl.40.50,
                                    mcar_mar_b_cov.60_mean.women.cl.40.50,
                                    mcar_mar_b_cov.65_mean.women.cl.40.50,
                                    mcar_mar_b_cov.70_mean.women.cl.40.50,
                                    mcar_mar_b_cov.75_mean.women.cl.40.50,
                                    mcar_mar_b_cov.80_mean.women.cl.40.50,
                                    mcar_mar_b_cov.85_mean.women.cl.40.50,
                                    mcar_mar_b_cov.90_mean.women.cl.40.50,
                                    mcar_mar_b_cov.95_mean.women.cl.40.50)


# median


mcar_mar_b_med.men.cl.15.25 <- c(mcar_mar_b_cov.35_med.men.cl.15.25,
                                 mcar_mar_b_cov.40_med.men.cl.15.25,
                                 mcar_mar_b_cov.45_med.men.cl.15.25,
                                 mcar_mar_b_cov.50_med.men.cl.15.25,
                                 mcar_mar_b_cov.55_med.men.cl.15.25,
                                 mcar_mar_b_cov.60_med.men.cl.15.25,
                                 mcar_mar_b_cov.65_med.men.cl.15.25,
                                 mcar_mar_b_cov.70_med.men.cl.15.25,
                                 mcar_mar_b_cov.75_med.men.cl.15.25,
                                 mcar_mar_b_cov.80_med.men.cl.15.25,
                                 mcar_mar_b_cov.85_med.men.cl.15.25,
                                 mcar_mar_b_cov.90_med.men.cl.15.25,
                                 mcar_mar_b_cov.95_med.men.cl.15.25)



mcar_mar_b_med.women.cl.15.25 <- c(mcar_mar_b_cov.35_med.women.cl.15.25,
                                   mcar_mar_b_cov.40_med.women.cl.15.25,
                                   mcar_mar_b_cov.45_med.women.cl.15.25,
                                   mcar_mar_b_cov.50_med.women.cl.15.25,
                                   mcar_mar_b_cov.55_med.women.cl.15.25,
                                   mcar_mar_b_cov.60_med.women.cl.15.25,
                                   mcar_mar_b_cov.65_med.women.cl.15.25,
                                   mcar_mar_b_cov.70_med.women.cl.15.25,
                                   mcar_mar_b_cov.75_med.women.cl.15.25,
                                   mcar_mar_b_cov.80_med.women.cl.15.25,
                                   mcar_mar_b_cov.85_med.women.cl.15.25,
                                   mcar_mar_b_cov.90_med.women.cl.15.25,
                                   mcar_mar_b_cov.95_med.women.cl.15.25)



mcar_mar_b_med.men.cl.25.40 <- c(mcar_mar_b_cov.35_med.men.cl.25.40,
                                 mcar_mar_b_cov.40_med.men.cl.25.40,
                                 mcar_mar_b_cov.45_med.men.cl.25.40,
                                 mcar_mar_b_cov.50_med.men.cl.25.40,
                                 mcar_mar_b_cov.55_med.men.cl.25.40,
                                 mcar_mar_b_cov.60_med.men.cl.25.40,
                                 mcar_mar_b_cov.65_med.men.cl.25.40,
                                 mcar_mar_b_cov.70_med.men.cl.25.40,
                                 mcar_mar_b_cov.75_med.men.cl.25.40,
                                 mcar_mar_b_cov.80_med.men.cl.25.40,
                                 mcar_mar_b_cov.85_med.men.cl.25.40,
                                 mcar_mar_b_cov.90_med.men.cl.25.40,
                                 mcar_mar_b_cov.95_med.men.cl.25.40)



mcar_mar_b_med.women.cl.25.40 <- c(mcar_mar_b_cov.35_med.women.cl.25.40,
                                   mcar_mar_b_cov.40_med.women.cl.25.40,
                                   mcar_mar_b_cov.45_med.women.cl.25.40,
                                   mcar_mar_b_cov.50_med.women.cl.25.40,
                                   mcar_mar_b_cov.55_med.women.cl.25.40,
                                   mcar_mar_b_cov.60_med.women.cl.25.40,
                                   mcar_mar_b_cov.65_med.women.cl.25.40,
                                   mcar_mar_b_cov.70_med.women.cl.25.40,
                                   mcar_mar_b_cov.75_med.women.cl.25.40,
                                   mcar_mar_b_cov.80_med.women.cl.25.40,
                                   mcar_mar_b_cov.85_med.women.cl.25.40,
                                   mcar_mar_b_cov.90_med.women.cl.25.40,
                                   mcar_mar_b_cov.95_med.women.cl.25.40)


mcar_mar_b_med.men.cl.40.50 <- c(mcar_mar_b_cov.35_med.men.cl.40.50,
                                 mcar_mar_b_cov.40_med.men.cl.40.50,
                                 mcar_mar_b_cov.45_med.men.cl.40.50,
                                 mcar_mar_b_cov.50_med.men.cl.40.50,
                                 mcar_mar_b_cov.55_med.men.cl.40.50,
                                 mcar_mar_b_cov.60_med.men.cl.40.50,
                                 mcar_mar_b_cov.65_med.men.cl.40.50,
                                 mcar_mar_b_cov.70_med.men.cl.40.50,
                                 mcar_mar_b_cov.75_med.men.cl.40.50,
                                 mcar_mar_b_cov.80_med.men.cl.40.50,
                                 mcar_mar_b_cov.85_med.men.cl.40.50,
                                 mcar_mar_b_cov.90_med.men.cl.40.50,
                                 mcar_mar_b_cov.95_med.men.cl.40.50)



mcar_mar_b_med.women.cl.40.50 <- c(mcar_mar_b_cov.35_med.women.cl.40.50,
                                   mcar_mar_b_cov.40_med.women.cl.40.50,
                                   mcar_mar_b_cov.45_med.women.cl.40.50,
                                   mcar_mar_b_cov.50_med.women.cl.40.50,
                                   mcar_mar_b_cov.55_med.women.cl.40.50,
                                   mcar_mar_b_cov.60_med.women.cl.40.50,
                                   mcar_mar_b_cov.65_med.women.cl.40.50,
                                   mcar_mar_b_cov.70_med.women.cl.40.50,
                                   mcar_mar_b_cov.75_med.women.cl.40.50,
                                   mcar_mar_b_cov.80_med.women.cl.40.50,
                                   mcar_mar_b_cov.85_med.women.cl.40.50,
                                   mcar_mar_b_cov.90_med.women.cl.40.50,
                                   mcar_mar_b_cov.95_med.women.cl.40.50)



# standard deviation

mcar_mar_b_sd.men.cl.15.25 <- c(mcar_mar_b_cov.35_sd.men.cl.15.25,
                                mcar_mar_b_cov.40_sd.men.cl.15.25,
                                mcar_mar_b_cov.45_sd.men.cl.15.25,
                                mcar_mar_b_cov.50_sd.men.cl.15.25,
                                mcar_mar_b_cov.55_sd.men.cl.15.25,
                                mcar_mar_b_cov.60_sd.men.cl.15.25,
                                mcar_mar_b_cov.65_sd.men.cl.15.25,
                                mcar_mar_b_cov.70_sd.men.cl.15.25,
                                mcar_mar_b_cov.75_sd.men.cl.15.25,
                                mcar_mar_b_cov.80_sd.men.cl.15.25,
                                mcar_mar_b_cov.85_sd.men.cl.15.25,
                                mcar_mar_b_cov.90_sd.men.cl.15.25,
                                mcar_mar_b_cov.95_sd.men.cl.15.25)



mcar_mar_b_sd.women.cl.15.25 <- c(mcar_mar_b_cov.35_sd.women.cl.15.25,
                                  mcar_mar_b_cov.40_sd.women.cl.15.25,
                                  mcar_mar_b_cov.45_sd.women.cl.15.25,
                                  mcar_mar_b_cov.50_sd.women.cl.15.25,
                                  mcar_mar_b_cov.55_sd.women.cl.15.25,
                                  mcar_mar_b_cov.60_sd.women.cl.15.25,
                                  mcar_mar_b_cov.65_sd.women.cl.15.25,
                                  mcar_mar_b_cov.70_sd.women.cl.15.25,
                                  mcar_mar_b_cov.75_sd.women.cl.15.25,
                                  mcar_mar_b_cov.80_sd.women.cl.15.25,
                                  mcar_mar_b_cov.85_sd.women.cl.15.25,
                                  mcar_mar_b_cov.90_sd.women.cl.15.25,
                                  mcar_mar_b_cov.95_sd.women.cl.15.25)



mcar_mar_b_sd.men.cl.25.40 <- c(mcar_mar_b_cov.35_sd.men.cl.25.40,
                                mcar_mar_b_cov.40_sd.men.cl.25.40,
                                mcar_mar_b_cov.45_sd.men.cl.25.40,
                                mcar_mar_b_cov.50_sd.men.cl.25.40,
                                mcar_mar_b_cov.55_sd.men.cl.25.40,
                                mcar_mar_b_cov.60_sd.men.cl.25.40,
                                mcar_mar_b_cov.65_sd.men.cl.25.40,
                                mcar_mar_b_cov.70_sd.men.cl.25.40,
                                mcar_mar_b_cov.75_sd.men.cl.25.40,
                                mcar_mar_b_cov.80_sd.men.cl.25.40,
                                mcar_mar_b_cov.85_sd.men.cl.25.40,
                                mcar_mar_b_cov.90_sd.men.cl.25.40,
                                mcar_mar_b_cov.95_sd.men.cl.25.40)


mcar_mar_b_sd.women.cl.25.40 <- c(mcar_mar_b_cov.35_sd.women.cl.25.40,
                                  mcar_mar_b_cov.40_sd.women.cl.25.40,
                                  mcar_mar_b_cov.45_sd.women.cl.25.40,
                                  mcar_mar_b_cov.50_sd.women.cl.25.40,
                                  mcar_mar_b_cov.55_sd.women.cl.25.40,
                                  mcar_mar_b_cov.60_sd.women.cl.25.40,
                                  mcar_mar_b_cov.65_sd.women.cl.25.40,
                                  mcar_mar_b_cov.70_sd.women.cl.25.40,
                                  mcar_mar_b_cov.75_sd.women.cl.25.40,
                                  mcar_mar_b_cov.80_sd.women.cl.25.40,
                                  mcar_mar_b_cov.85_sd.women.cl.25.40,
                                  mcar_mar_b_cov.90_sd.women.cl.25.40,
                                  mcar_mar_b_cov.95_sd.women.cl.25.40)



mcar_mar_b_sd.men.cl.40.50 <- c(mcar_mar_b_cov.35_sd.men.cl.40.50,
                                mcar_mar_b_cov.40_sd.men.cl.40.50,
                                mcar_mar_b_cov.45_sd.men.cl.40.50,
                                mcar_mar_b_cov.50_sd.men.cl.40.50,
                                mcar_mar_b_cov.55_sd.men.cl.40.50,
                                mcar_mar_b_cov.60_sd.men.cl.40.50,
                                mcar_mar_b_cov.65_sd.men.cl.40.50,
                                mcar_mar_b_cov.70_sd.men.cl.40.50,
                                mcar_mar_b_cov.75_sd.men.cl.40.50,
                                mcar_mar_b_cov.80_sd.men.cl.40.50,
                                mcar_mar_b_cov.85_sd.men.cl.40.50,
                                mcar_mar_b_cov.90_sd.men.cl.40.50,
                                mcar_mar_b_cov.95_sd.men.cl.40.50)



mcar_mar_b_sd.women.cl.40.50 <- c(mcar_mar_b_cov.35_sd.women.cl.40.50,
                                  mcar_mar_b_cov.40_sd.women.cl.40.50,
                                  mcar_mar_b_cov.45_sd.women.cl.40.50,
                                  mcar_mar_b_cov.50_sd.women.cl.40.50,
                                  mcar_mar_b_cov.55_sd.women.cl.40.50,
                                  mcar_mar_b_cov.60_sd.women.cl.40.50,
                                  mcar_mar_b_cov.65_sd.women.cl.40.50,
                                  mcar_mar_b_cov.70_sd.women.cl.40.50,
                                  mcar_mar_b_cov.75_sd.women.cl.40.50,
                                  mcar_mar_b_cov.80_sd.women.cl.40.50,
                                  mcar_mar_b_cov.85_sd.women.cl.40.50,
                                  mcar_mar_b_cov.90_sd.women.cl.40.50,
                                  mcar_mar_b_cov.95_sd.women.cl.40.50)





# MCAR - MAR - c AD -------------- 






# mean

mcar_mar_c_mean.men.cl.15.25 <- c(mcar_mar_c_cov.35_mean.men.cl.15.25,
                                  mcar_mar_c_cov.40_mean.men.cl.15.25,
                                  mcar_mar_c_cov.45_mean.men.cl.15.25,
                                  mcar_mar_c_cov.50_mean.men.cl.15.25,
                                  mcar_mar_c_cov.55_mean.men.cl.15.25,
                                  mcar_mar_c_cov.60_mean.men.cl.15.25,
                                  mcar_mar_c_cov.65_mean.men.cl.15.25,
                                  mcar_mar_c_cov.70_mean.men.cl.15.25,
                                  mcar_mar_c_cov.75_mean.men.cl.15.25,
                                  mcar_mar_c_cov.80_mean.men.cl.15.25,
                                  mcar_mar_c_cov.85_mean.men.cl.15.25,
                                  mcar_mar_c_cov.90_mean.men.cl.15.25,
                                  mcar_mar_c_cov.95_mean.men.cl.15.25)



mcar_mar_c_mean.women.cl.15.25 <- c(mcar_mar_c_cov.35_mean.women.cl.15.25,
                                    mcar_mar_c_cov.40_mean.women.cl.15.25,
                                    mcar_mar_c_cov.45_mean.women.cl.15.25,
                                    mcar_mar_c_cov.50_mean.women.cl.15.25,
                                    mcar_mar_c_cov.55_mean.women.cl.15.25,
                                    mcar_mar_c_cov.60_mean.women.cl.15.25,
                                    mcar_mar_c_cov.65_mean.women.cl.15.25,
                                    mcar_mar_c_cov.70_mean.women.cl.15.25,
                                    mcar_mar_c_cov.75_mean.women.cl.15.25,
                                    mcar_mar_c_cov.80_mean.women.cl.15.25,
                                    mcar_mar_c_cov.85_mean.women.cl.15.25,
                                    mcar_mar_c_cov.90_mean.women.cl.15.25,
                                    mcar_mar_c_cov.95_mean.women.cl.15.25)



mcar_mar_c_mean.men.cl.25.40 <- c(mcar_mar_c_cov.35_mean.men.cl.25.40,
                                  mcar_mar_c_cov.40_mean.men.cl.25.40,
                                  mcar_mar_c_cov.45_mean.men.cl.25.40,
                                  mcar_mar_c_cov.50_mean.men.cl.25.40,
                                  mcar_mar_c_cov.55_mean.men.cl.25.40,
                                  mcar_mar_c_cov.60_mean.men.cl.25.40,
                                  mcar_mar_c_cov.65_mean.men.cl.25.40,
                                  mcar_mar_c_cov.70_mean.men.cl.25.40,
                                  mcar_mar_c_cov.75_mean.men.cl.25.40,
                                  mcar_mar_c_cov.80_mean.men.cl.25.40,
                                  mcar_mar_c_cov.85_mean.men.cl.25.40,
                                  mcar_mar_c_cov.90_mean.men.cl.25.40,
                                  mcar_mar_c_cov.95_mean.men.cl.25.40)


mcar_mar_c_mean.women.cl.25.40 <- c(mcar_mar_c_cov.35_mean.women.cl.25.40,
                                    mcar_mar_c_cov.40_mean.women.cl.25.40,
                                    mcar_mar_c_cov.45_mean.women.cl.25.40,
                                    mcar_mar_c_cov.50_mean.women.cl.25.40,
                                    mcar_mar_c_cov.55_mean.women.cl.25.40,
                                    mcar_mar_c_cov.60_mean.women.cl.25.40,
                                    mcar_mar_c_cov.65_mean.women.cl.25.40,
                                    mcar_mar_c_cov.70_mean.women.cl.25.40,
                                    mcar_mar_c_cov.75_mean.women.cl.25.40,
                                    mcar_mar_c_cov.80_mean.women.cl.25.40,
                                    mcar_mar_c_cov.85_mean.women.cl.25.40,
                                    mcar_mar_c_cov.90_mean.women.cl.25.40,
                                    mcar_mar_c_cov.95_mean.women.cl.25.40)




mcar_mar_c_mean.men.cl.40.50 <- c(mcar_mar_c_cov.35_mean.men.cl.40.50,
                                  mcar_mar_c_cov.40_mean.men.cl.40.50,
                                  mcar_mar_c_cov.45_mean.men.cl.40.50,
                                  mcar_mar_c_cov.50_mean.men.cl.40.50,
                                  mcar_mar_c_cov.55_mean.men.cl.40.50,
                                  mcar_mar_c_cov.60_mean.men.cl.40.50,
                                  mcar_mar_c_cov.65_mean.men.cl.40.50,
                                  mcar_mar_c_cov.70_mean.men.cl.40.50,
                                  mcar_mar_c_cov.75_mean.men.cl.40.50,
                                  mcar_mar_c_cov.80_mean.men.cl.40.50,
                                  mcar_mar_c_cov.85_mean.men.cl.40.50,
                                  mcar_mar_c_cov.90_mean.men.cl.40.50,
                                  mcar_mar_c_cov.95_mean.men.cl.40.50)



mcar_mar_c_mean.women.cl.40.50 <- c(mcar_mar_c_cov.35_mean.women.cl.40.50,
                                    mcar_mar_c_cov.40_mean.women.cl.40.50,
                                    mcar_mar_c_cov.45_mean.women.cl.40.50,
                                    mcar_mar_c_cov.50_mean.women.cl.40.50,
                                    mcar_mar_c_cov.55_mean.women.cl.40.50,
                                    mcar_mar_c_cov.60_mean.women.cl.40.50,
                                    mcar_mar_c_cov.65_mean.women.cl.40.50,
                                    mcar_mar_c_cov.70_mean.women.cl.40.50,
                                    mcar_mar_c_cov.75_mean.women.cl.40.50,
                                    mcar_mar_c_cov.80_mean.women.cl.40.50,
                                    mcar_mar_c_cov.85_mean.women.cl.40.50,
                                    mcar_mar_c_cov.90_mean.women.cl.40.50,
                                    mcar_mar_c_cov.95_mean.women.cl.40.50)


# median


mcar_mar_c_med.men.cl.15.25 <- c(mcar_mar_c_cov.35_med.men.cl.15.25,
                                 mcar_mar_c_cov.40_med.men.cl.15.25,
                                 mcar_mar_c_cov.45_med.men.cl.15.25,
                                 mcar_mar_c_cov.50_med.men.cl.15.25,
                                 mcar_mar_c_cov.55_med.men.cl.15.25,
                                 mcar_mar_c_cov.60_med.men.cl.15.25,
                                 mcar_mar_c_cov.65_med.men.cl.15.25,
                                 mcar_mar_c_cov.70_med.men.cl.15.25,
                                 mcar_mar_c_cov.75_med.men.cl.15.25,
                                 mcar_mar_c_cov.80_med.men.cl.15.25,
                                 mcar_mar_c_cov.85_med.men.cl.15.25,
                                 mcar_mar_c_cov.90_med.men.cl.15.25,
                                 mcar_mar_c_cov.95_med.men.cl.15.25)



mcar_mar_c_med.women.cl.15.25 <- c(mcar_mar_c_cov.35_med.women.cl.15.25,
                                   mcar_mar_c_cov.40_med.women.cl.15.25,
                                   mcar_mar_c_cov.45_med.women.cl.15.25,
                                   mcar_mar_c_cov.50_med.women.cl.15.25,
                                   mcar_mar_c_cov.55_med.women.cl.15.25,
                                   mcar_mar_c_cov.60_med.women.cl.15.25,
                                   mcar_mar_c_cov.65_med.women.cl.15.25,
                                   mcar_mar_c_cov.70_med.women.cl.15.25,
                                   mcar_mar_c_cov.75_med.women.cl.15.25,
                                   mcar_mar_c_cov.80_med.women.cl.15.25,
                                   mcar_mar_c_cov.85_med.women.cl.15.25,
                                   mcar_mar_c_cov.90_med.women.cl.15.25,
                                   mcar_mar_c_cov.95_med.women.cl.15.25)



mcar_mar_c_med.men.cl.25.40 <- c(mcar_mar_c_cov.35_med.men.cl.25.40,
                                 mcar_mar_c_cov.40_med.men.cl.25.40,
                                 mcar_mar_c_cov.45_med.men.cl.25.40,
                                 mcar_mar_c_cov.50_med.men.cl.25.40,
                                 mcar_mar_c_cov.55_med.men.cl.25.40,
                                 mcar_mar_c_cov.60_med.men.cl.25.40,
                                 mcar_mar_c_cov.65_med.men.cl.25.40,
                                 mcar_mar_c_cov.70_med.men.cl.25.40,
                                 mcar_mar_c_cov.75_med.men.cl.25.40,
                                 mcar_mar_c_cov.80_med.men.cl.25.40,
                                 mcar_mar_c_cov.85_med.men.cl.25.40,
                                 mcar_mar_c_cov.90_med.men.cl.25.40,
                                 mcar_mar_c_cov.95_med.men.cl.25.40)




mcar_mar_c_med.women.cl.25.40 <- c(mcar_mar_c_cov.35_med.women.cl.25.40,
                                   mcar_mar_c_cov.40_med.women.cl.25.40,
                                   mcar_mar_c_cov.45_med.women.cl.25.40,
                                   mcar_mar_c_cov.50_med.women.cl.25.40,
                                   mcar_mar_c_cov.55_med.women.cl.25.40,
                                   mcar_mar_c_cov.60_med.women.cl.25.40,
                                   mcar_mar_c_cov.65_med.women.cl.25.40,
                                   mcar_mar_c_cov.70_med.women.cl.25.40,
                                   mcar_mar_c_cov.75_med.women.cl.25.40,
                                   mcar_mar_c_cov.80_med.women.cl.25.40,
                                   mcar_mar_c_cov.85_med.women.cl.25.40,
                                   mcar_mar_c_cov.90_med.women.cl.25.40,
                                   mcar_mar_c_cov.95_med.women.cl.25.40)


mcar_mar_c_med.men.cl.40.50 <- c(mcar_mar_c_cov.35_med.men.cl.40.50,
                                 mcar_mar_c_cov.40_med.men.cl.40.50,
                                 mcar_mar_c_cov.45_med.men.cl.40.50,
                                 mcar_mar_c_cov.50_med.men.cl.40.50,
                                 mcar_mar_c_cov.55_med.men.cl.40.50,
                                 mcar_mar_c_cov.60_med.men.cl.40.50,
                                 mcar_mar_c_cov.65_med.men.cl.40.50,
                                 mcar_mar_c_cov.70_med.men.cl.40.50,
                                 mcar_mar_c_cov.75_med.men.cl.40.50,
                                 mcar_mar_c_cov.80_med.men.cl.40.50,
                                 mcar_mar_c_cov.85_med.men.cl.40.50,
                                 mcar_mar_c_cov.90_med.men.cl.40.50,
                                 mcar_mar_c_cov.95_med.men.cl.40.50)



mcar_mar_c_med.women.cl.40.50 <- c(mcar_mar_c_cov.35_med.women.cl.40.50,
                                   mcar_mar_c_cov.40_med.women.cl.40.50,
                                   mcar_mar_c_cov.45_med.women.cl.40.50,
                                   mcar_mar_c_cov.50_med.women.cl.40.50,
                                   mcar_mar_c_cov.55_med.women.cl.40.50,
                                   mcar_mar_c_cov.60_med.women.cl.40.50,
                                   mcar_mar_c_cov.65_med.women.cl.40.50,
                                   mcar_mar_c_cov.70_med.women.cl.40.50,
                                   mcar_mar_c_cov.75_med.women.cl.40.50,
                                   mcar_mar_c_cov.80_med.women.cl.40.50,
                                   mcar_mar_c_cov.85_med.women.cl.40.50,
                                   mcar_mar_c_cov.90_med.women.cl.40.50,
                                   mcar_mar_c_cov.95_med.women.cl.40.50)



# standard deviation

mcar_mar_c_sd.men.cl.15.25 <- c(mcar_mar_c_cov.35_sd.men.cl.15.25,
                                mcar_mar_c_cov.40_sd.men.cl.15.25,
                                mcar_mar_c_cov.45_sd.men.cl.15.25,
                                mcar_mar_c_cov.50_sd.men.cl.15.25,
                                mcar_mar_c_cov.55_sd.men.cl.15.25,
                                mcar_mar_c_cov.60_sd.men.cl.15.25,
                                mcar_mar_c_cov.65_sd.men.cl.15.25,
                                mcar_mar_c_cov.70_sd.men.cl.15.25,
                                mcar_mar_c_cov.75_sd.men.cl.15.25,
                                mcar_mar_c_cov.80_sd.men.cl.15.25,
                                mcar_mar_c_cov.85_sd.men.cl.15.25,
                                mcar_mar_c_cov.90_sd.men.cl.15.25,
                                mcar_mar_c_cov.95_sd.men.cl.15.25)



mcar_mar_c_sd.women.cl.15.25 <- c(mcar_mar_c_cov.35_sd.women.cl.15.25,
                                  mcar_mar_c_cov.40_sd.women.cl.15.25,
                                  mcar_mar_c_cov.45_sd.women.cl.15.25,
                                  mcar_mar_c_cov.50_sd.women.cl.15.25,
                                  mcar_mar_c_cov.55_sd.women.cl.15.25,
                                  mcar_mar_c_cov.60_sd.women.cl.15.25,
                                  mcar_mar_c_cov.65_sd.women.cl.15.25,
                                  mcar_mar_c_cov.70_sd.women.cl.15.25,
                                  mcar_mar_c_cov.75_sd.women.cl.15.25,
                                  mcar_mar_c_cov.80_sd.women.cl.15.25,
                                  mcar_mar_c_cov.85_sd.women.cl.15.25,
                                  mcar_mar_c_cov.90_sd.women.cl.15.25,
                                  mcar_mar_c_cov.95_sd.women.cl.15.25)



mcar_mar_c_sd.men.cl.25.40 <- c(mcar_mar_c_cov.35_sd.men.cl.25.40,
                                mcar_mar_c_cov.40_sd.men.cl.25.40,
                                mcar_mar_c_cov.45_sd.men.cl.25.40,
                                mcar_mar_c_cov.50_sd.men.cl.25.40,
                                mcar_mar_c_cov.55_sd.men.cl.25.40,
                                mcar_mar_c_cov.60_sd.men.cl.25.40,
                                mcar_mar_c_cov.65_sd.men.cl.25.40,
                                mcar_mar_c_cov.70_sd.men.cl.25.40,
                                mcar_mar_c_cov.75_sd.men.cl.25.40,
                                mcar_mar_c_cov.80_sd.men.cl.25.40,
                                mcar_mar_c_cov.85_sd.men.cl.25.40,
                                mcar_mar_c_cov.90_sd.men.cl.25.40,
                                mcar_mar_c_cov.95_sd.men.cl.25.40)



mcar_mar_c_sd.women.cl.25.40 <- c(mcar_mar_c_cov.35_sd.women.cl.25.40,
                                  mcar_mar_c_cov.40_sd.women.cl.25.40,
                                  mcar_mar_c_cov.45_sd.women.cl.25.40,
                                  mcar_mar_c_cov.50_sd.women.cl.25.40,
                                  mcar_mar_c_cov.55_sd.women.cl.25.40,
                                  mcar_mar_c_cov.60_sd.women.cl.25.40,
                                  mcar_mar_c_cov.65_sd.women.cl.25.40,
                                  mcar_mar_c_cov.70_sd.women.cl.25.40,
                                  mcar_mar_c_cov.75_sd.women.cl.25.40,
                                  mcar_mar_c_cov.80_sd.women.cl.25.40,
                                  mcar_mar_c_cov.85_sd.women.cl.25.40,
                                  mcar_mar_c_cov.90_sd.women.cl.25.40,
                                  mcar_mar_c_cov.95_sd.women.cl.25.40)



mcar_mar_c_sd.men.cl.40.50 <- c(mcar_mar_c_cov.35_sd.men.cl.40.50,
                                mcar_mar_c_cov.40_sd.men.cl.40.50,
                                mcar_mar_c_cov.45_sd.men.cl.40.50,
                                mcar_mar_c_cov.50_sd.men.cl.40.50,
                                mcar_mar_c_cov.55_sd.men.cl.40.50,
                                mcar_mar_c_cov.60_sd.men.cl.40.50,
                                mcar_mar_c_cov.65_sd.men.cl.40.50,
                                mcar_mar_c_cov.70_sd.men.cl.40.50,
                                mcar_mar_c_cov.75_sd.men.cl.40.50,
                                mcar_mar_c_cov.80_sd.men.cl.40.50,
                                mcar_mar_c_cov.85_sd.men.cl.40.50,
                                mcar_mar_c_cov.90_sd.men.cl.40.50,
                                mcar_mar_c_cov.95_sd.men.cl.40.50)



mcar_mar_c_sd.women.cl.40.50 <- c(mcar_mar_c_cov.35_sd.women.cl.40.50,
                                  mcar_mar_c_cov.40_sd.women.cl.40.50,
                                  mcar_mar_c_cov.45_sd.women.cl.40.50,
                                  mcar_mar_c_cov.50_sd.women.cl.40.50,
                                  mcar_mar_c_cov.55_sd.women.cl.40.50,
                                  mcar_mar_c_cov.60_sd.women.cl.40.50,
                                  mcar_mar_c_cov.65_sd.women.cl.40.50,
                                  mcar_mar_c_cov.70_sd.women.cl.40.50,
                                  mcar_mar_c_cov.75_sd.women.cl.40.50,
                                  mcar_mar_c_cov.80_sd.women.cl.40.50,
                                  mcar_mar_c_cov.85_sd.women.cl.40.50,
                                  mcar_mar_c_cov.90_sd.women.cl.40.50,
                                  mcar_mar_c_cov.95_sd.women.cl.40.50)



# Rows for the final table --------------


# 
# # Proportions -------------
# 
# # a
# 
# mcar_mar_a_prop.men15.25.F.15.25
# mcar_mar_a_prop.women15.25.M.15.25
# mcar_mar_a_prop.men25.40.F.15.25
# mcar_mar_a_prop.women15.25.M.25.40
# mcar_mar_a_prop.men25.40.F.25.40
# mcar_mar_a_prop.women25.40.M.25.40
# mcar_mar_a_prop.men40.50.F.15.25
# mcar_mar_a_prop.women15.25.M.40.50
# mcar_mar_a_prop.men40.50.F.25.40
# mcar_mar_a_prop.women25.40.M.40.50
# 
# # b
# 
# mcar_mar_b_prop.men15.25.F.15.25
# mcar_mar_b_prop.women15.25.M.15.25
# mcar_mar_b_prop.men25.40.F.15.25
# mcar_mar_b_prop.women15.25.M.25.40
# mcar_mar_b_prop.men25.40.F.25.40
# mcar_mar_b_prop.women25.40.M.25.40
# mcar_mar_b_prop.men40.50.F.15.25
# mcar_mar_b_prop.women15.25.M.40.50
# mcar_mar_b_prop.men40.50.F.25.40
# mcar_mar_b_prop.women25.40.M.40.50
# 
# 
# # c
# 
# mcar_mar_c_prop.men15.25.F.15.25
# mcar_mar_c_prop.women15.25.M.15.25
# mcar_mar_c_prop.men25.40.F.15.25
# mcar_mar_c_prop.women15.25.M.25.40
# mcar_mar_c_prop.men25.40.F.25.40
# mcar_mar_c_prop.women25.40.M.25.40
# mcar_mar_c_prop.men40.50.F.15.25
# mcar_mar_c_prop.women15.25.M.40.50
# mcar_mar_c_prop.men40.50.F.25.40
# mcar_mar_c_prop.women25.40.M.40.50
# 
# 
# 



# # Age difference --------------


# 
# # a
# 
# # mean
# 
# mcar_mar_a_mean.men.cl.15.25
# mcar_mar_a_mean.women.cl.15.25
# mcar_mar_a_mean.men.cl.25.40
# mcar_mar_a_mean.women.cl.15.25
# mcar_mar_a_mean.men.cl.40.50
# mcar_mar_a_mean.women.cl.40.50
# 
# # median
# 
# mcar_mar_a_med.men.cl.15.25
# mcar_mar_a_med.women.cl.15.25
# mcar_mar_a_med.men.cl.25.40
# mcar_mar_a_med.women.cl.15.25
# mcar_mar_a_med.men.cl.40.50
# mcar_mar_a_med.women.cl.40.50
# 
# 
# # standard deviation
# 
# mcar_mar_a_sd.men.cl.15.25
# mcar_mar_a_sd.women.cl.15.25
# mcar_mar_a_sd.men.cl.25.40
# mcar_mar_a_sd.women.cl.15.25
# mcar_mar_a_sd.men.cl.40.50
# mcar_mar_a_sd.women.cl.40.50
# 
# 
# # b
# 
# 
# # mean
# 
# mcar_mar_b_mean.men.cl.15.25
# mcar_mar_b_mean.women.cl.15.25
# mcar_mar_b_mean.men.cl.25.40
# mcar_mar_b_mean.women.cl.15.25
# mcar_mar_b_mean.men.cl.40.50
# mcar_mar_b_mean.women.cl.40.50
# 
# # median
# 
# mcar_mar_b_med.men.cl.15.25
# mcar_mar_b_med.women.cl.15.25
# mcar_mar_b_med.men.cl.25.40
# mcar_mar_b_med.women.cl.15.25
# mcar_mar_b_med.men.cl.40.50
# mcar_mar_b_med.women.cl.40.50
# 
# 
# # standard deviation
# 
# mcar_mar_b_sd.men.cl.15.25
# mcar_mar_b_sd.women.cl.15.25
# mcar_mar_b_sd.men.cl.25.40
# mcar_mar_b_sd.women.cl.15.25
# mcar_mar_b_sd.men.cl.40.50
# mcar_mar_b_sd.women.cl.40.50
# 
# 
# # c
# 
# 
# # mean
# 
# mcar_mar_c_mean.men.cl.15.25
# mcar_mar_c_mean.women.cl.15.25
# mcar_mar_c_mean.men.cl.25.40
# mcar_mar_c_mean.women.cl.15.25
# mcar_mar_c_mean.men.cl.40.50
# mcar_mar_c_mean.women.cl.40.50
# 
# # median
# 
# mcar_mar_c_med.men.cl.15.25
# mcar_mar_c_med.women.cl.15.25
# mcar_mar_c_med.men.cl.25.40
# mcar_mar_c_med.women.cl.15.25
# mcar_mar_c_med.men.cl.40.50
# mcar_mar_c_med.women.cl.40.50
# 
# 
# # standard deviation
# 
# mcar_mar_c_sd.men.cl.15.25
# mcar_mar_c_sd.women.cl.15.25
# mcar_mar_c_sd.men.cl.25.40
# mcar_mar_c_sd.women.cl.15.25
# mcar_mar_c_sd.men.cl.40.50
# mcar_mar_c_sd.women.cl.40.50
# 



final_table <- matrix(c(
  
  # a
  
  mcar_mar_a_prop.men15.25.F.15.25,
  mcar_mar_a_prop.women15.25.M.15.25,
  mcar_mar_a_prop.men25.40.F.15.25,
  mcar_mar_a_prop.women15.25.M.25.40,
  mcar_mar_a_prop.men25.40.F.25.40,
  mcar_mar_a_prop.women25.40.M.25.40,
  mcar_mar_a_prop.men40.50.F.15.25,
  mcar_mar_a_prop.women15.25.M.40.50,
  mcar_mar_a_prop.men40.50.F.25.40,
  mcar_mar_a_prop.women25.40.M.40.50,
  
  # b
  
  mcar_mar_b_prop.men15.25.F.15.25,
  mcar_mar_b_prop.women15.25.M.15.25,
  mcar_mar_b_prop.men25.40.F.15.25,
  mcar_mar_b_prop.women15.25.M.25.40,
  mcar_mar_b_prop.men25.40.F.25.40,
  mcar_mar_b_prop.women25.40.M.25.40,
  mcar_mar_b_prop.men40.50.F.15.25,
  mcar_mar_b_prop.women15.25.M.40.50,
  mcar_mar_b_prop.men40.50.F.25.40,
  mcar_mar_b_prop.women25.40.M.40.50,
  
  
  # c
  
  mcar_mar_c_prop.men15.25.F.15.25,
  mcar_mar_c_prop.women15.25.M.15.25,
  mcar_mar_c_prop.men25.40.F.15.25,
  mcar_mar_c_prop.women15.25.M.25.40,
  mcar_mar_c_prop.men25.40.F.25.40,
  mcar_mar_c_prop.women25.40.M.25.40,
  mcar_mar_c_prop.men40.50.F.15.25,
  mcar_mar_c_prop.women15.25.M.40.50,
  mcar_mar_c_prop.men40.50.F.25.40,
  mcar_mar_c_prop.women25.40.M.40.50,
  
  
  # a
  
  # mean
  
  mcar_mar_a_mean.men.cl.15.25,
  mcar_mar_a_mean.women.cl.15.25,
  mcar_mar_a_mean.men.cl.25.40,
  mcar_mar_a_mean.women.cl.25.40,
  mcar_mar_a_mean.men.cl.40.50,
  mcar_mar_a_mean.women.cl.40.50,
  
  # median
  
  mcar_mar_a_med.men.cl.15.25,
  mcar_mar_a_med.women.cl.15.25,
  mcar_mar_a_med.men.cl.25.40,
  mcar_mar_a_med.women.cl.25.40,
  mcar_mar_a_med.men.cl.40.50,
  mcar_mar_a_med.women.cl.40.50,
  
  
  # standard deviation
  
  mcar_mar_a_sd.men.cl.15.25,
  mcar_mar_a_sd.women.cl.15.25,
  mcar_mar_a_sd.men.cl.25.40,
  mcar_mar_a_sd.women.cl.25.40,
  mcar_mar_a_sd.men.cl.40.50,
  mcar_mar_a_sd.women.cl.40.50,
  
  
  # b
  
  
  # mean
  
  mcar_mar_b_mean.men.cl.15.25,
  mcar_mar_b_mean.women.cl.15.25,
  mcar_mar_b_mean.men.cl.25.40,
  mcar_mar_b_mean.women.cl.25.40,
  mcar_mar_b_mean.men.cl.40.50,
  mcar_mar_b_mean.women.cl.40.50,
  
  # median
  
  mcar_mar_b_med.men.cl.15.25,
  mcar_mar_b_med.women.cl.15.25,
  mcar_mar_b_med.men.cl.25.40,
  mcar_mar_b_med.women.cl.25.40,
  mcar_mar_b_med.men.cl.40.50,
  mcar_mar_b_med.women.cl.40.50,
  
  
  # standard deviation
  
  mcar_mar_b_sd.men.cl.15.25,
  mcar_mar_b_sd.women.cl.15.25,
  mcar_mar_b_sd.men.cl.25.40,
  mcar_mar_b_sd.women.cl.25.40,
  mcar_mar_b_sd.men.cl.40.50,
  mcar_mar_b_sd.women.cl.40.50,
  
  
  # c
  
  
  # mean
  
  mcar_mar_c_mean.men.cl.15.25,
  mcar_mar_c_mean.women.cl.15.25,
  mcar_mar_c_mean.men.cl.25.40,
  mcar_mar_c_mean.women.cl.25.40,
  mcar_mar_c_mean.men.cl.40.50,
  mcar_mar_c_mean.women.cl.40.50,
  
  # median
  
  mcar_mar_c_med.men.cl.15.25,
  mcar_mar_c_med.women.cl.15.25,
  mcar_mar_c_med.men.cl.25.40,
  mcar_mar_c_med.women.cl.25.40,
  mcar_mar_c_med.men.cl.40.50,
  mcar_mar_c_med.women.cl.40.50,
  
  
  # standard deviation
  
  mcar_mar_c_sd.men.cl.15.25,
  mcar_mar_c_sd.women.cl.15.25,
  mcar_mar_c_sd.men.cl.25.40,
  mcar_mar_c_sd.women.cl.25.40,
  mcar_mar_c_sd.men.cl.40.50,
  mcar_mar_c_sd.women.cl.40.50
),
ncol = 13,
byrow = TRUE)



rownames(final_table) <- c(
  
  "mcar_mar_a_prop.M.15.24.F.15.24",
  "mcar_mar_a_prop.F.15.24.M.15.24",
  "mcar_mar_a_prop.M.25.39.F.15.24",
  "mcar_mar_a_prop.F.15.24.M.25.39",
  "mcar_mar_a_prop.M.25.39.F.25.39",
  "mcar_mar_a_prop.F.25.39.M.25.39",
  "mcar_mar_a_prop.M.40.49.F.15.24",
  "mcar_mar_a_prop.F.15.24.M.40.49",
  "mcar_mar_a_prop.M.40.49.F.25.39",
  "mcar_mar_a_prop.F.25.39.M.40.49",
  
  "mcar_mar_b_prop.M.15.24.F.15.24",
  "mcar_mar_b_prop.F.15.24.M.15.24",
  "mcar_mar_b_prop.M.25.39.F.15.24",
  "mcar_mar_b_prop.F.15.24.M.25.39",
  "mcar_mar_b_prop.M.25.39.F.25.39",
  "mcar_mar_b_prop.F.25.39.M.25.39",
  "mcar_mar_b_prop.M.40.49.F.15.24",
  "mcar_mar_b_prop.F.15.24.M.40.49",
  "mcar_mar_b_prop.M.40.49.F.25.39",
  "mcar_mar_b_prop.F.25.39.M.40.49",
  
  "mcar_mar_c_prop.M.15.24.F.15.24",
  "mcar_mar_c_prop.F.15.24.M.15.24",
  "mcar_mar_c_prop.M.25.39.F.15.24",
  "mcar_mar_c_prop.F.15.24.M.25.39",
  "mcar_mar_c_prop.M.25.39.F.25.39",
  "mcar_mar_c_prop.F.25.39.M.25.39",
  "mcar_mar_c_prop.M.40.49.F.15.24",
  "mcar_mar_c_prop.F.15.24.M.40.49",
  "mcar_mar_c_prop.M.40.49.F.25.39",
  "mcar_mar_c_prop.F.25.39.M.40.49",
  
  "mcar_mar_a_mean.M.15.24",
  "mcar_mar_a_mean.F.15.24",
  "mcar_mar_a_mean.M.25.39",
  "mcar_mar_a_mean.F.25.39",
  "mcar_mar_a_mean.M.40.49",
  "mcar_mar_a_mean.F.40.49",
  
  
  "mcar_mar_a_med.M.15.24",
  "mcar_mar_a_med.F.15.24",
  "mcar_mar_a_med.M.25.39",
  "mcar_mar_a_med.F.25.39",
  "mcar_mar_a_med.M.40.49",
  "mcar_mar_a_med.F.40.49",
  
  
  "mcar_mar_a_sd.M.15.24",
  "mcar_mar_a_sd.F.15.24",
  "mcar_mar_a_sd.M.25.39",
  "mcar_mar_a_sd.F.25.39",
  "mcar_mar_a_sd.M.40.49",
  "mcar_mar_a_sd.F.40.49",
  
  
  "mcar_mar_b_mean.M.15.24",
  "mcar_mar_b_mean.F.15.24",
  "mcar_mar_b_mean.M.25.39",
  "mcar_mar_b_mean.F.25.39",
  "mcar_mar_b_mean.M.40.49",
  "mcar_mar_b_mean.F.40.49",
  
  
  "mcar_mar_b_med.M.15.24",
  "mcar_mar_b_med.F.15.24",
  "mcar_mar_b_med.M.25.39",
  "mcar_mar_b_med.F.25.39",
  "mcar_mar_b_med.M.40.49",
  "mcar_mar_b_med.F.40.49",
  
  
  "mcar_mar_b_sd.M.15.24",
  "mcar_mar_b_sd.F.15.24",
  "mcar_mar_b_sd.M.25.39",
  "mcar_mar_b_sd.F.25.39",
  "mcar_mar_b_sd.M.40.49",
  "mcar_mar_b_sd.F.40.49",
  
  
  "mcar_mar_c_mean.M.15.24",
  "mcar_mar_c_mean.F.15.24",
  "mcar_mar_c_mean.M.25.39",
  "mcar_mar_c_mean.F.25.39",
  "mcar_mar_c_mean.M.40.49",
  "mcar_mar_c_mean.F.40.49",
  
  
  "mcar_mar_c_med.M.15.24",
  "mcar_mar_c_med.F.15.24",
  "mcar_mar_c_med.M.25.39",
  "mcar_mar_c_med.F.25.39",
  "mcar_mar_c_med.M.40.49",
  "mcar_mar_c_med.F.40.49",
  
  
  "mcar_mar_c_sd.M.15.24",
  "mcar_mar_c_sd.F.15.24",
  "mcar_mar_c_sd.M.25.39",
  "mcar_mar_c_sd.F.25.39",
  "mcar_mar_c_sd.M.40.49",
  "mcar_mar_c_sd.F.40.49")


colnames(final_table) <- c("cov.35", "cov.40", "cov.45", "cov.50", "cov.55",
                           "cov.60", "cov.65", "cov.70", "cov.75", "cov.80",
                           "cov.85", "cov.90", "cov.95")

saveRDS(final_table, file = "/home/david/age_mixing_patterns_phylogenetic/results/MCAR_MAR_comparison/compar_mcar_mar_table.RDS")



final_table <- readRDS("/home/david/age_mixing_patterns_phylogenetic/results/MCAR_MAR_comparison/compar_mcar_mar_table.RDS")


# final_table %>%
#   kable() %>%
#   kable_styling("striped") # Commented OCTOBER


final_table_proportions <- final_table[1:30,]

# final_table_proportions <- round(final_table_proportions, digits = 3)


# write.csv(final_table_proportions, file = "/home/david/age_mixing_patterns_phylogenetic/results/MCAR_MAR_comparison/Table_comp_proportions_mcar_mar.csv")
# 
# 
# final_table_proportions <- read.csv("/home/david/age_mixing_patterns_phylogenetic/results/MCAR_MAR_comparison/Table_comp_proportions_mcar_mar.csv")
# 


final_table_age_difference <- final_table[31:nrow(final_table),]

# final_table_age_difference <- round(final_table_age_difference, digits = 3)



# write.csv(final_table_age_difference, file = "/home/david/age_mixing_patterns_phylogenetic/results/MCAR_MAR_comparison/Table_comp_age_difference_mcar_mar.csv")
# 

# final_table_age_difference <- read.csv("/home/david/age_mixing_patterns_phylogenetic/results/MCAR_MAR_comparison/Table_comp_age_difference_mcar_mar.csv")
# 



final_table_proportions <- as.data.frame(final_table_proportions)
final_table_age_difference <- as.data.frame(final_table_age_difference)


# props_mar_a <-  final_table_proportions %>%
#   select(contains("mcar_mar_a_")) 


# For proportions of pairings
#############################

props_mar_a <-  as.matrix(final_table_proportions[1:10,])
y_1 <- rownames(props_mar_a)
x_1 <- sub("mcar_mar_a_prop.", "", y_1)
rownames(props_mar_a) <- x_1


props_mar_b <-  as.matrix(final_table_proportions[11:20,])

rownames(props_mar_b) <- x_1


props_mar_c <-  as.matrix(final_table_proportions[21:30,])

rownames(props_mar_c) <- x_1



# For age difference
####################

ad_a <- as.matrix(final_table_age_difference[1:18,])

ad_b <- as.matrix(final_table_age_difference[19:36,])

ad_c <- as.matrix(final_table_age_difference[37:54,])


mean_ad_a <- as.matrix(ad_a[1:6,])
sd_ad_a <- as.matrix(ad_a[13:18,])

y_2 <- rownames(mean_ad_a)
x_2 <- sub("mcar_mar_a_mean.", "", y_2)
rownames(mean_ad_a) <- x_2
rownames(sd_ad_a) <- x_2


mean_ad_b <- as.matrix(ad_b[1:6,])
sd_ad_b <- as.matrix(ad_b[13:18,])

rownames(mean_ad_b) <- x_2
rownames(sd_ad_b) <- x_2


mean_ad_c <- as.matrix(ad_c[1:6,])
sd_ad_c <- as.matrix(ad_c[13:18,])

rownames(mean_ad_c) <- x_2
rownames(sd_ad_c) <- x_2



# Plot heatmaps


library(ComplexHeatmap) # from Bioconductor  
library(circlize)


heatmap_props_mar_a <- Heatmap(props_mar_a, 
                               name = "P-values",
                               show_heatmap_legend = TRUE,
                               show_row_names = TRUE,
                               cluster_rows = FALSE,
                               cluster_columns = FALSE)

heatmap_props_mar_b <- Heatmap(props_mar_b, 
                               # name = "P-values",
                               show_heatmap_legend = FALSE,
                               show_row_names = FALSE,
                               cluster_rows  =FALSE,
                               cluster_columns = FALSE)

heatmap_props_mar_c <- Heatmap(props_mar_c, 
                               # name = "P-values",
                               show_heatmap_legend = FALSE,
                               show_row_names = FALSE,
                               cluster_rows  =FALSE,
                               cluster_columns = FALSE)


heatmap_mean_mar_a <- Heatmap(mean_ad_a, 
                              name = "P-values",
                              show_heatmap_legend = TRUE,
                              show_row_names = TRUE,
                              cluster_rows  =FALSE,
                              cluster_columns = FALSE)

heatmap_mean_mar_b <- Heatmap(mean_ad_b, 
                              # name = "P-values",
                              show_heatmap_legend = FALSE,
                              show_row_names = FALSE,
                              cluster_rows  =FALSE,
                              cluster_columns = FALSE)


heatmap_mean_mar_c <- Heatmap(mean_ad_c, 
                              # name = "P-values",
                              show_heatmap_legend = FALSE,
                              show_row_names = FALSE,
                              cluster_rows  =FALSE,
                              cluster_columns = FALSE)



heatmap_sd_mar_a <- Heatmap(sd_ad_a, 
                            name = "P-values",
                            show_heatmap_legend = TRUE,
                            show_row_names = TRUE,
                            cluster_rows  =FALSE,
                            cluster_columns = FALSE)

heatmap_sd_mar_b <- Heatmap(sd_ad_b, 
                            # name = "P-values",
                            show_heatmap_legend = FALSE,
                            show_row_names = FALSE,
                            cluster_rows  =FALSE,
                            cluster_columns = FALSE)

heatmap_sd_mar_c <- Heatmap(sd_ad_c, 
                            # name = "P-values",
                            show_heatmap_legend = FALSE,
                            show_row_names = FALSE,
                            cluster_rows  =FALSE,
                            cluster_columns = FALSE)


# b - 30, c - 50, a - 70

plot.heatmapsproportions <- heatmap_props_mar_b + heatmap_props_mar_c + heatmap_props_mar_a

# 12 x 6 inches

ggsave(filename = "Plot_00_heatmapsproportions.pdf",
       plot = plot.heatmapsproportions,
       path = "/home/david/age_mixing_patterns_phylogenetic/results/MCAR_MAR_comparison",
       width = 16, height = 10, units = "cm")


plot.heatmapsmeanagegap <- heatmap_mean_mar_b + heatmap_mean_mar_c + heatmap_mean_mar_a

ggsave(filename = "Plot_00_heatmapsmeanagegap.pdf",
       plot = plot.heatmapsmeanagegap,
       path = "/home/david/age_mixing_patterns_phylogenetic/results/MCAR_MAR_comparison",
       width = 16, height = 10, units = "cm")


plot.heatmapsSDagegap <- heatmap_sd_mar_b + heatmap_sd_mar_c + heatmap_sd_mar_a

ggsave(filename = "Plot_00_heatmapsSDagegap.pdf",
       plot = plot.heatmapsSDagegap,
       path = "/home/david/age_mixing_patterns_phylogenetic/results/MCAR_MAR_comparison",
       width = 16, height = 10, units = "cm")

