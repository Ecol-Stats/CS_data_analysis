library(ggplot2)
library(lubridate)
library(tidyverse)
library(INLA)
library(fmesher)
library(sf)

rm(list = ls())

inla.setOption(num.threads="2:4")

#these I have downloaded from the github repo of the RSS paper
load(here::here("Data/ringletdata.rda"))
load(here::here("Data/maxListLAreaRinglet.rda"))


inla.Occupancy_detCov <- function(X_det){
  
  if(class(X_det)=="list"){
    if(length(X_det)>10){
      warning("exceeded number of detection covariates, numerical issues may occur")
    }
    
    if(lapply(X_det, ncol)%>%unlist()%>%unique()%>%length()>2){
      stop("inconsistent number of visits in provided detection covariates")
    }
    if(length(lapply(X_det, nrow) %>% unlist() %>% unique())>1){
      stop("inconsistent number of sites in provided detection covariates")
    }
    K<- lapply(X_det, ncol) %>% unlist() %>% max() # Max num of visits
    M<- lapply(X_det, nrow) %>% unlist() %>% unique() # Number of sites
    P <- length(X_det)
    
    if(lapply(X_det, ncol)%>%unlist()%>%unique()%>%length()==2 & 
       1 %in% lapply(X_det, ncol)%>%unlist()%>%unique()){
      warning(paste("At least one covariate of dimension [",M,",1] has been provided, values for this covariate will be repeated over the max numver of visits",sep=""))
      for(l in which(lapply(X_det, ncol) %>% unlist() < K)){
        X_det[[l]] <- do.call("cbind",replicate(K,X_det[[l]]))
        
      }
    }
    covariates <- do.call("cbind", lapply(1:K, function(i) {
      do.call("cbind", lapply(X_det, function(mat) mat[, i]))
    }))
    
  }
  
  if(is.data.frame(X_det)|is.matrix(X_det)){
    K<- ncol(X_det)
    M<- nrow(X_det)
    P <- 1
    covariates <- as.matrix(X_det)
  }
  
  X_mat <- matrix(NA,nrow=M,ncol=K*(P+1))
  X_mat[,seq(1,(K*(P+1)),by=(P+1))]<-1 # add Intercept at the begining of each visit-specific covariate matrix
  X_mat[, which(!(1:(K*(P+1)) %in% seq(1,(K*(P+1)),by=(P+1))))] <- covariates
  return(X_mat)
  
}

# data cleaning -----------------------------------------------------------

{
  data$maxListlArea <- maxListlArea
  
  data <- data %>%
    select(-X) %>%
    rename(Site = Gridref) %>% 
    mutate(SamplUnit = paste(Site, Year, sep = " - "),
           Date = as.Date(Date)) %>% 
    mutate(month  = month(Date),
           jDate  = yday(Date)) %>%
    arrange(Year, Site) %>%
    mutate(x = EAST/1000,
           y = NORTH/1000)
  
  # Add detection covariates
  data$relativeListLength <- data$listL / data$maxListlArea
  
  data$relativeListLength <- (data$relativeListLength-
                                mean(data$relativeListLength)) / sd(data$relativeListLength)
  
  
  
}



## Occupancy data ----------------------------------------------------------


data_filter = data %>% 
  filter(month %in% c(5,6,7,8)) %>%
  group_by(SamplUnit) %>% 
  mutate(yy = seq_along(x)) %>% 
  ungroup() 

data_filter  = data_filter %>%
  mutate(JulianDate1 = scale(jDate),
         JulianDate2 = scale(jDate)^2,
         JulianDate3 = scale(jDate)^3,
         relativeListLength = scale(relativeListLength))

Occ = data_filter%>%
  select(SamplUnit, Year, x,y, Occ, yy) %>%
  pivot_wider(names_from = yy, values_from = Occ)


Y =  Occ[,-c(1:4)]
X_Occ = Occ[,1:4]




## detection data ----------------------------------------------------------


julian1 = data_filter%>%
  select(SamplUnit,JulianDate1, yy) %>%
  pivot_wider(names_from = yy, values_from = JulianDate1)
julian2 = data_filter%>%
  select(SamplUnit,JulianDate2, yy) %>%
  pivot_wider(names_from = yy, values_from = JulianDate2)
julian3 = data_filter%>%
  select(SamplUnit,JulianDate3, yy) %>%
  pivot_wider(names_from = yy, values_from = JulianDate3)
relativeListLength = data_filter%>%
  select(SamplUnit,relativeListLength, yy) %>%
  pivot_wider(names_from = yy, values_from = relativeListLength)

X_det <- list(
  as.matrix(julian1[,-1]) ,
  as.matrix(julian2[,-1]),
  # as.matrix(julian3[,-1]),
  as.matrix(relativeListLength[,-1])
) %>% 
  inla.Occupancy_detCov()



# temporal model without spatial effect  ------------------------------------------------------------

if(0)
{
  # This is hte model in the RSS-A paper 
  
  X_det <- list(
    as.matrix(julian1[,-1]) ,
    as.matrix(julian2[,-1]),
    #as.matrix(julian3[,-1]),
    as.matrix(relativeListLength[,-1])
  ) %>% 
    inla.Occupancy_detCov()
  
  
  
  X_occ2 <- Occ[,c(1:4)] %>% 
    mutate(Int_occ =1) 
  
  
  data_list = as.list(X_occ2)
  data_list$Y = Y
  data_list$X = X_det
  
  formula_time <- inla.mdata(Y,X) ~  -1 + Int_occ +  f(Year, model = "iid", 
                                                       fixed = F, initial  = 0, constr = T)
  model_time <- inla(formula_time,    # model formula
                     data= data_list,        # data 
                     family= 'occupancy', # model likelihood
                     #control.inla = list(int.strategy = "eb"),
                     
                     # priors
                     control.fixed =  list(prec = 1, prec.intercept = 1),
                     control.mode = list(theta = c(0.13521223 , 0.09538584, -5.31917496,
                                                   1.75581385,  0)),
                     # compute WAIC and DIC
                     control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                     verbose = TRUE,
                     # choose link functions for:
                     # (i) the state process (control.link)
                     # (ii) the observation process (link.simple)
                     control.family = list(control.link = list(model = "logit"),
                                           link.simple = "logit",
                                           # priors for hyperparameters
                                           hyper = list(
                                             # Detection intercept prior 
                                             beta1 = list(param = c(0,100)),
                                             # Covariate 1 effect prior
                                             beta2 = list(param = c(0,100)),
                                             # Covariate 2 effect prior 
                                             beta3 = list(param = c(0,100)),
                                             beta4 = list(param = c(0,100)),
                                             beta5 = list(param = c(0,100))
                                           )
                     ))
  
  samples_time = inla.posterior.sample(300, model_time)
  
  ## year spec intercept -----------------------------------------------------
  
  func_int = function(...)
  {
    out = inla.link.invlogit(Year)
    out
  }
  
  out = inla.posterior.sample.eval(func_int, samples_time) 
  time_effect = data.frame(time = sort(unique(data$Year)),
                           mean = apply(out,1,quantile, 0.5),
                           q1 = apply(out,1,quantile, 0.025),
                           q2 = apply(out,1,quantile, 0.975))
  time_effect %>%
    ggplot() + 
    geom_point(aes(time,mean)) +
    geom_errorbar(aes(time, ymin = q1, 
                      ymax =  q2)) +
    ylim(0,1) + ggtitle("Year specific intercept") +
    ylab("Probability of Occurrence") + xlab("")
  
  ## Detection function ------------------------------------------------------
  
  
  xx = yday(as.Date("01/05/2012", "%d/%m/%y")):yday(as.Date("31/08/2012", "%d/%m/%y"))
  dd = cbind(1, 
             scale(xx),
             scale(xx)^2,
             scale(xx)^3)    
  
  
  data.frame(x = parse_date_time(x = xx, orders = "j") ,
             m3 = inla.link.invlogit(dd[,1:3] %*% model_time$summary.hyperpar$mean[1:3]),
             q12 = inla.link.invlogit(dd[,1:3] %*% model_time$summary.hyperpar$`0.025quant`[1:3]),
             q22 = inla.link.invlogit(dd[,1:3] %*% model_time$summary.hyperpar$`0.025quant`[1:3])) %>%
    ggplot() + 
    geom_line(aes(x,m3)) +
    geom_ribbon(aes(x, ymin = q12, ymax = q22), alpha = 0.5) +
    ylab("Probability") + xlab("") + ggtitle("Probability of Detection")
  
}




# create the Mesh -------------------------------------------------------

data_sf = st_as_sf(data_filter %>% distinct(Site, .keep_all = TRUE), coords = c("x","y"))
boundary = fm_nonconvex_hull(data_sf, convex = -0.03)

mesh = fm_mesh_2d(boundary = boundary, max.edge = c(25,120),
                  offset = c(10,80))
plot(mesh)




# FULL MODEL  --------------------------------------------------------------
# This model adds a spatial effect for each year

# NB: This takes quite some memory, so it is recommended that one runs it on a "new" R session
# without running the time model
matern <- inla.spde2.pcmatern(mesh,
                              prior.range = c(150, 0.5),
                              prior.sigma = c(1, 0.5))


X_det <- list(
  as.matrix(julian1[,-1]) ,
  as.matrix(julian2[,-1]),
  # as.matrix(julian3[,-1]),
  as.matrix(relativeListLength[,-1])
) %>% 
  inla.Occupancy_detCov()



X_occ2 <- X_Occ %>% 
  mutate(Int_occ =1,
         spatial_field = rep(NA,nrow(X_Occ)),
         time = Year-min(Year)+1) %>%
  mutate(x_scale = scale(  x)[,1],
         y_scale = scale(  y)[,1]) %>%
  mutate(Year1 = Year,
         Year2  = Year)

A_sp2 <- inla.spde.make.A(mesh = mesh, 
                          group = X_occ2$time,
                          loc = cbind(X_occ2$x, X_occ2$y))


data_list = as.list(X_occ2)
data_list$Y = Y
data_list$X = X_det




# INLA implementation -----------------------------------------------------

formula2 <- inla.mdata(Y,X) ~   -1 + #Int_occ +  
  f(spatial_field, model=matern,  A.local = A_sp2, replicate = time) +
  f(Year, x_scale, model= "iid", initial = 0, fixed = T,  constr = F) + 
  f(Year1, y_scale, model= "iid", initial = 0, fixed = T,  constr = F) + 
  f(Year2, model= "iid", initial = 0, fixed = T,  constr = F)

model2 <- inla(formula2,    # model formula
               data= data_list,        # data 
               family= 'occupancy', # model likelihood
               # priors
               #control.mode = list(theta = c(-0.02448643 , 0.73152176, -3.92154786,
               #                                1.18437337, 4.96780413, 1.07614563)),
               control.fixed =  list(prec = 1, prec.intercept = 1),
               #control.inla=list(int.strategy='eb'),
               # compute WAIC and DIC
               control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
               verbose = TRUE,
               # choose link functions for:
               # (i) the state process (control.link)
               # (ii) the observation process (link.simple)
               control.family = list(control.link = list(model = "logit"),
                                     link.simple = "logit",
                                     # priors for hyperparameters
                                     hyper = list(
                                       # Detection intercept prior
                                       beta1 = list(param = c(0,500)),
                                       # Covariate 1 effect prior
                                       beta2 = list(param = c(0,500)),
                                       # Covariate 2 effect prior
                                       beta3 = list(param = c(0,500)),
                                       beta4 = list(param = c(0,500)),
                                       # Covariate 1 effect prior
                                       beta5 = list(param = c(0,500)),
                                       # Covariate 2 effect prior
                                       beta6 = list(param = c(0,500))
                                     ))
)

