library(tidyverse)
library(mirt)
library(psych)


###############################################################################           
##############################    sub functions  ############################## 
###############################################################################           
P_LCDM <- function(par, Theta, ncat){
  # 1-7 intercept and main effects
  # 8-15 interaction effects
  pstar <- c(par[1:7], 
             par[8] - min(par[c(2, 6)]),
             par[9] - min(par[c(2, 7)]),
             par[10] - min(par[c(3, 6)]),
             par[11] - min(par[c(3, 7)]),
             par[12] - min(par[c(4, 6)]),
             par[13] - min(par[c(4, 7)]),
             par[14] - min(par[c(5, 6)]),
             par[15] - min(par[c(5, 7)])
  )
  
  z <- as.vector(Theta %*% pstar)
  p <- plogis(z)
  cbind(1-p, p)
}

get_index <- function(attrs, all_attrs){
  if(missing(all_attrs)){
    all_attrs <- colnames(Q_mat_6_lcdm)
  }
  index <- NULL
  for(attr in attrs){
    index <- c(index, 
               which(all_attrs == attr))
  }
  return(index)
}

set_P_LCDM <- function(attr){
  n_col_full = 1 + n_main + n_interaction
  index <- get_index(attr)
  par <- rep(0, n_col_full)
  names(par) <- paste('a', 1:n_col_full, sep='')
  est <- rep(FALSE, n_col_full)
  lbound <- rep(-Inf, n_col_full)
  est[1] <- TRUE
  
  par[index[index <= n_main + 1]] <- 1
  par[index[index > n_main + 1]] <- .5
  est[index] <- TRUE
  lbound[index] <- 0
  
  return(list(par, est, lbound))
}

createItem_customize <- function(attrs){
  if(length(attrs) > 1){
    attrs <- c(attrs, paste(attrs[1], attrs[2], sep ="_"))
  }
  item_names <- paste('LCDM', attrs[length(attrs)], sep = '_')
  
  item <- createItem(item_names,
                     par = set_P_LCDM(attrs)[[1]],
                     est = set_P_LCDM(attrs)[[2]],
                     P=P_LCDM,
                     lbound = set_P_LCDM(attrs)[[3]])
}



###############################################################################           
##############################       main        ############################## 
############################################################################### 

load('data.RData')

# the total number of main effects and interaction effects may be contained
# within the lcdm model
n_main <<- 6 # 6 main effects/number of attributes
n_interaction <<- 8 # number of potential interaction effects


Q_mat_6_lcdm <<- Q_mat_6_new

theta6 <- rep(1, 2^5)
for(i in 0:5){
  theta6 <- cbind(theta6, 
                  rep(c(0,1), each = 2^5/(2^i), times = 2^i))
}

theta_lcdm <- theta6
names_lcdm <- colnames(Q_mat_6_new)

# set the full Q_matrix contained interaction effects
for(i in 2:5){
  for(j in 6:7){
    Q_mat_6_lcdm <- cbind(Q_mat_6_lcdm,
                          Q_mat_6_lcdm[,i]*Q_mat_6_lcdm[,j])
    names_lcdm <- c(names_lcdm, paste(names_lcdm[i], names_lcdm[j], sep="_"))
    theta_lcdm <- cbind(theta_lcdm, theta6[,i]*theta6[,j])
  }
}
colnames(Q_mat_6_lcdm) <- names_lcdm

# specify the mirt.model
str_attributes <- NULL

for(name in colnames(Q_mat_6_lcdm)[-1]){
  str_attribute <- paste(name, ' = ', 
                         paste(which(Q_mat_6_lcdm[, name] == 1), 
                               collapse = ','),
                         sep = "")
  if(is.null(str_attributes)){
    str_attributes <- str_attribute
  } else {
    str_attributes <- paste(str_attributes, 
                            str_attribute, sep = "\n")
  }
}

s6_lcdm <- paste('intercept = 1-80', str_attributes,
                 sep = "\n")

# set the customized item

itemtype <- NULL
for(i in 1:ncol(response)){
  attrs <- colnames(Q_mat_6_lcdm)[Q_mat_6_lcdm[i, ] == 1][-1]
  itemtype <- c(itemtype,
                paste('LCDM', attrs[length(attrs)], sep = '_'))
}

customItems <- list(
  LCDM_fraction = createItem_customize('fraction'),
  LCDM_decimal = createItem_customize('decimal'),
  LCDM_t2_fraction = createItem_customize(c('t2', 'fraction')),
  LCDM_t2_decimal = createItem_customize(c('t2', 'decimal')),
  LCDM_t3_fraction = createItem_customize(c('t3', 'fraction')),
  LCDM_t3_decimal = createItem_customize(c('t3', 'decimal')),
  LCDM_t4_fraction = createItem_customize(c('t4', 'fraction')),
  LCDM_t4_decimal = createItem_customize(c('t4', 'decimal')),
  LCDM_t5_fraction = createItem_customize(c('t5', 'fraction')),
  LCDM_t5_decimal = createItem_customize(c('t5', 'decimal'))
)

mod_lcdm <- mdirt(response*1, model = s6_lcdm,
                  customTheta = theta_lcdm, 
                  itemtype = itemtype,
                  customItems=customItems, 
                  TOL = 1e-3,
                  technical = list(NCYCLES = 1000))

#saveRDS(mod_lcdm, file = file.path(save_path, '6_model_lcdm.rds'))

## ***************** classification ***************** ##
mod_lcdm_pv <- fscores(mod_lcdm)[,-1] > .5
mod_lcdm_pv <- as.data.frame(mod_lcdm_pv)
names(mod_lcdm_pv) <- colnames(Q_mat_6_new)[-1]
#write.csv(mod_lcdm_pv[,1:6], row.names = FALSE,
#          file.path(save_path, "m6_lcdm_classification.csv"))

## ************** tetrachoric correlation ************** ##
col_lcdm_pv <- tetrachoric(mod_lcdm_pv, delete = FALSE)$rho %>%
  round(4)
print(col_lcdm_pv[1:6, 1:6])

## ************** item parameters ************** ##
cfs_mod_lcdm <- coef(mod_lcdm, simplify=TRUE)$items
colnames(cfs_mod_lcdm) <- names_lcdm

