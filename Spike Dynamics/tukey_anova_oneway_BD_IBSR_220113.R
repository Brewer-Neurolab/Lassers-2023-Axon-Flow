
# Well BD -----------------------------------------------------------------



#Burst Duration
#ec

ec<-read.csv("D:\\Brewer lab data\\HFS\\well_5SD_500maxSD_BD_moment_tbl_EC.csv")

Mean <- c(ec$mean)
SD <- c(ec$sd)
SampleSize <- c(ec$n)

gen_data <- function(means, sds, samplesizes){
  n.grp <- length(means)
  grps <- factor(rep(1:n.grp, samplesizes))
  dat <- lapply(1:n.grp, function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
  y <- do.call(rbind, dat)
  out <- data.frame(group = grps, y = y)
  out
}

simulated_data <- gen_data(Mean, SD,SampleSize)
av <- aov(y ~ group, data = simulated_data)
summary(av)
ancova_tbl<-TukeyHSD(av)
ancova_tbl$group

#dg
dg<-read.csv("D:\\Brewer lab data\\HFS\\well_5SD_500maxSD_BD_moment_tbl_DG.csv")

Mean <- c(dg$mean)
SD <- c(dg$sd)
SampleSize <- c(dg$n)

gen_data <- function(means, sds, samplesizes){
  n.grp <- length(means)
  grps <- factor(rep(1:n.grp, samplesizes))
  dat <- lapply(1:n.grp, function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
  y <- do.call(rbind, dat)
  out <- data.frame(group = grps, y = y)
  out
}

simulated_data <- gen_data(Mean, SD,SampleSize)
av <- aov(y ~ group, data = simulated_data)
summary(av)
ancova_tbl<-TukeyHSD(av)
ancova_tbl$group

#ca3
ca3<-read.csv("D:\\Brewer lab data\\HFS\\well_5SD_500maxSD_BD_moment_tbl_CA3.csv")

Mean <- c(ca3$mean)
SD <- c(ca3$sd)
SampleSize <- c(ca3$n)

gen_data <- function(means, sds, samplesizes){
  n.grp <- length(means)
  grps <- factor(rep(1:n.grp, samplesizes))
  dat <- lapply(1:n.grp, function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
  y <- do.call(rbind, dat)
  out <- data.frame(group = grps, y = y)
  out
}

simulated_data <- gen_data(Mean, SD,SampleSize)
av <- aov(y ~ group, data = simulated_data)
summary(av)
ancova_tbl<-TukeyHSD(av)
ancova_tbl$group

#ca1
ca1<-read.csv("D:\\Brewer lab data\\HFS\\well_5SD_500maxSD_BD_moment_tbl_CA1.csv")

Mean <- c(ca1$mean)
SD <- c(ca1$sd)
SampleSize <- c(ca1$n)

gen_data <- function(means, sds, samplesizes){
  n.grp <- length(means)
  grps <- factor(rep(1:n.grp, samplesizes))
  dat <- lapply(1:n.grp, function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
  y <- do.call(rbind, dat)
  out <- data.frame(group = grps, y = y)
  out
}

simulated_data <- gen_data(Mean, SD,SampleSize)
av <- aov(y ~ group, data = simulated_data)
summary(av)
ancova_tbl<-TukeyHSD(av)
ancova_tbl$group

# FF BD -------------------------------------------------------------------
#Burst Duration
#ec-dg

ec<-read.csv("D:\\Brewer lab data\\HFS\\ff_BD_moment_tbl_EC-DG.csv")

Mean <- c(ec$mean)
SD <- c(ec$sd)
SampleSize <- c(ec$n)

gen_data <- function(means, sds, samplesizes){
  n.grp <- length(means)
  grps <- factor(rep(1:n.grp, samplesizes))
  dat <- lapply(1:n.grp, function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
  y <- do.call(rbind, dat)
  out <- data.frame(group = grps, y = y)
  out
}

simulated_data <- gen_data(Mean, SD,SampleSize)
av <- aov(y ~ group, data = simulated_data)
summary(av)
TukeyHSD(av)

#dg-ca3
dg<-read.csv("D:\\Brewer lab data\\HFS\\ff_BD_moment_tbl_DG-CA3.csv")

Mean <- c(dg$mean)
SD <- c(dg$sd)
SampleSize <- c(dg$n)

gen_data <- function(means, sds, samplesizes){
  n.grp <- length(means)
  grps <- factor(rep(1:n.grp, samplesizes))
  dat <- lapply(1:n.grp, function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
  y <- do.call(rbind, dat)
  out <- data.frame(group = grps, y = y)
  out
}

simulated_data <- gen_data(Mean, SD,SampleSize)
av <- aov(y ~ group, data = simulated_data)
summary(av)
TukeyHSD(av)

#ca3-ca1
ca3<-read.csv("D:\\Brewer lab data\\HFS\\ff_BD_moment_tbl_CA3-CA1.csv")

Mean <- c(ca3$mean)
SD <- c(ca3$sd)
SampleSize <- c(ca3$n)

gen_data <- function(means, sds, samplesizes){
  n.grp <- length(means)
  grps <- factor(rep(1:n.grp, samplesizes))
  dat <- lapply(1:n.grp, function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
  y <- do.call(rbind, dat)
  out <- data.frame(group = grps, y = y)
  out
}

simulated_data <- gen_data(Mean, SD,SampleSize)
av <- aov(y ~ group, data = simulated_data)
summary(av)
TukeyHSD(av)

#ca1-ec
ca1<-read.csv("D:\\Brewer lab data\\HFS\\ff_BD_moment_tbl_CA1-EC.csv")

Mean <- c(ca1$mean)
SD <- c(ca1$sd)
SampleSize <- c(ca1$n)

gen_data <- function(means, sds, samplesizes){
  n.grp <- length(means)
  grps <- factor(rep(1:n.grp, samplesizes))
  dat <- lapply(1:n.grp, function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
  y <- do.call(rbind, dat)
  out <- data.frame(group = grps, y = y)
  out
}

simulated_data <- gen_data(Mean, SD,SampleSize)
av <- aov(y ~ group, data = simulated_data)
summary(av)
TukeyHSD(av)

# FB BD -------------------------------------------------------------------

#ec-dg

ec<-read.csv("D:\\Brewer lab data\\HFS\\fb_BD_moment_tbl_DG-EC.csv")

Mean <- c(ec$mean)
SD <- c(ec$sd)
SampleSize <- c(ec$n)

gen_data <- function(means, sds, samplesizes){
  n.grp <- length(means)
  grps <- factor(rep(1:n.grp, samplesizes))
  dat <- lapply(1:n.grp, function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
  y <- do.call(rbind, dat)
  out <- data.frame(group = grps, y = y)
  out
}

simulated_data <- gen_data(Mean, SD,SampleSize)
av <- aov(y ~ group, data = simulated_data)
summary(av)
TukeyHSD(av)

#dg-ca3
dg<-read.csv("D:\\Brewer lab data\\HFS\\fb_BD_moment_tbl_CA3-DG.csv")

Mean <- c(dg$mean)
SD <- c(dg$sd)
SampleSize <- c(dg$n)

gen_data <- function(means, sds, samplesizes){
  n.grp <- length(means)
  grps <- factor(rep(1:n.grp, samplesizes))
  dat <- lapply(1:n.grp, function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
  y <- do.call(rbind, dat)
  out <- data.frame(group = grps, y = y)
  out
}

simulated_data <- gen_data(Mean, SD,SampleSize)
av <- aov(y ~ group, data = simulated_data)
summary(av)
TukeyHSD(av)

#ca3-ca1
ca3<-read.csv("D:\\Brewer lab data\\HFS\\fb_BD_moment_tbl_CA1-CA3.csv")

Mean <- c(ca3$mean)
SD <- c(ca3$sd)
SampleSize <- c(ca3$n)

gen_data <- function(means, sds, samplesizes){
  n.grp <- length(means)
  grps <- factor(rep(1:n.grp, samplesizes))
  dat <- lapply(1:n.grp, function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
  y <- do.call(rbind, dat)
  out <- data.frame(group = grps, y = y)
  out
}

simulated_data <- gen_data(Mean, SD,SampleSize)
av <- aov(y ~ group, data = simulated_data)
summary(av)
TukeyHSD(av)

#ca1-ec
ca1<-read.csv("D:\\Brewer lab data\\HFS\\fb_BD_moment_tbl_EC-CA1.csv")

Mean <- c(ca1$mean)
SD <- c(ca1$sd)
SampleSize <- c(ca1$n)

gen_data <- function(means, sds, samplesizes){
  n.grp <- length(means)
  grps <- factor(rep(1:n.grp, samplesizes))
  dat <- lapply(1:n.grp, function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
  y <- do.call(rbind, dat)
  out <- data.frame(group = grps, y = y)
  out
}

simulated_data <- gen_data(Mean, SD,SampleSize)
av <- aov(y ~ group, data = simulated_data)
summary(av)
TukeyHSD(av)


# Well IBSR ---------------------------------------------------------------



#IBSR
#ec

ec<-read.csv("D:\\Brewer lab data\\HFS\\well_5SD_500maxSD_IBSR_moment_tbl_EC.csv")

Mean <- c(ec$mean)
SD <- c(ec$sd)
SampleSize <- c(ec$n)

gen_data <- function(means, sds, samplesizes){
  n.grp <- length(means)
  grps <- factor(rep(1:n.grp, samplesizes))
  dat <- lapply(1:n.grp, function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
  y <- do.call(rbind, dat)
  out <- data.frame(group = grps, y = y)
  out
}

simulated_data <- gen_data(Mean, SD,SampleSize)
av <- aov(y ~ group, data = simulated_data)
summary(av)
ancova_tbl<-TukeyHSD(av)
ancova_tbl$group

#dg
dg<-read.csv("D:\\Brewer lab data\\HFS\\well_5SD_500maxSD_IBSR_moment_tbl_DG.csv")

Mean <- c(dg$mean)
SD <- c(dg$sd)
SampleSize <- c(dg$n)

gen_data <- function(means, sds, samplesizes){
  n.grp <- length(means)
  grps <- factor(rep(1:n.grp, samplesizes))
  dat <- lapply(1:n.grp, function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
  y <- do.call(rbind, dat)
  out <- data.frame(group = grps, y = y)
  out
}

simulated_data <- gen_data(Mean, SD,SampleSize)
av <- aov(y ~ group, data = simulated_data)
summary(av)
ancova_tbl<-TukeyHSD(av)
ancova_tbl$group

#ca3
ca3<-read.csv("D:\\Brewer lab data\\HFS\\well_5SD_500maxSD_IBSR_moment_tbl_CA3.csv")

Mean <- c(ca3$mean)
SD <- c(ca3$sd)
SampleSize <- c(ca3$n)

gen_data <- function(means, sds, samplesizes){
  n.grp <- length(means)
  grps <- factor(rep(1:n.grp, samplesizes))
  dat <- lapply(1:n.grp, function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
  y <- do.call(rbind, dat)
  out <- data.frame(group = grps, y = y)
  out
}

simulated_data <- gen_data(Mean, SD,SampleSize)
av <- aov(y ~ group, data = simulated_data)
summary(av)
ancova_tbl<-TukeyHSD(av)
ancova_tbl$group

#ca1
ca1<-read.csv("D:\\Brewer lab data\\HFS\\well_5SD_500maxSD_IBSR_moment_tbl_CA1.csv")

Mean <- c(ca1$mean)
SD <- c(ca1$sd)
SampleSize <- c(ca1$n)

gen_data <- function(means, sds, samplesizes){
  n.grp <- length(means)
  grps <- factor(rep(1:n.grp, samplesizes))
  dat <- lapply(1:n.grp, function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
  y <- do.call(rbind, dat)
  out <- data.frame(group = grps, y = y)
  out
}

simulated_data <- gen_data(Mean, SD,SampleSize)
av <- aov(y ~ group, data = simulated_data)
summary(av)
ancova_tbl<-TukeyHSD(av)
ancova_tbl$group

# FF IBSR -----------------------------------------------------------------
#ec-dg

ec<-read.csv("D:\\Brewer lab data\\HFS\\ff_IBSR_moment_tbl_EC-DG.csv")

Mean <- c(ec$mean)
SD <- c(ec$sd)
SampleSize <- c(ec$n)

gen_data <- function(means, sds, samplesizes){
  n.grp <- length(means)
  grps <- factor(rep(1:n.grp, samplesizes))
  dat <- lapply(1:n.grp, function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
  y <- do.call(rbind, dat)
  out <- data.frame(group = grps, y = y)
  out
}

simulated_data <- gen_data(Mean, SD,SampleSize)
av <- aov(y ~ group, data = simulated_data)
summary(av)
TukeyHSD(av)

#dg-ca3
dg<-read.csv("D:\\Brewer lab data\\HFS\\ff_IBSR_moment_tbl_DG-CA3.csv")

Mean <- c(dg$mean)
SD <- c(dg$sd)
SampleSize <- c(dg$n)

gen_data <- function(means, sds, samplesizes){
  n.grp <- length(means)
  grps <- factor(rep(1:n.grp, samplesizes))
  dat <- lapply(1:n.grp, function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
  y <- do.call(rbind, dat)
  out <- data.frame(group = grps, y = y)
  out
}

simulated_data <- gen_data(Mean, SD,SampleSize)
av <- aov(y ~ group, data = simulated_data)
summary(av)
TukeyHSD(av)

#ca3-ca1
ca3<-read.csv("D:\\Brewer lab data\\HFS\\ff_IBSR_moment_tbl_CA3-CA1.csv")

Mean <- c(ca3$mean)
SD <- c(ca3$sd)
SampleSize <- c(ca3$n)

gen_data <- function(means, sds, samplesizes){
  n.grp <- length(means)
  grps <- factor(rep(1:n.grp, samplesizes))
  dat <- lapply(1:n.grp, function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
  y <- do.call(rbind, dat)
  out <- data.frame(group = grps, y = y)
  out
}

simulated_data <- gen_data(Mean, SD,SampleSize)
av <- aov(y ~ group, data = simulated_data)
summary(av)
TukeyHSD(av)

#ca1-ec
ca1<-read.csv("D:\\Brewer lab data\\HFS\\ff_IBSR_moment_tbl_CA1-EC.csv")

Mean <- c(ca1$mean)
SD <- c(ca1$sd)
SampleSize <- c(ca1$n)

gen_data <- function(means, sds, samplesizes){
  n.grp <- length(means)
  grps <- factor(rep(1:n.grp, samplesizes))
  dat <- lapply(1:n.grp, function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
  y <- do.call(rbind, dat)
  out <- data.frame(group = grps, y = y)
  out
}

simulated_data <- gen_data(Mean, SD,SampleSize)
av <- aov(y ~ group, data = simulated_data)
summary(av)
TukeyHSD(av)

# FB IBSR -----------------------------------------------------------------
#ec-dg

ec<-read.csv("D:\\Brewer lab data\\HFS\\fb_IBSR_moment_tbl_EC-DG.csv")

Mean <- c(ec$mean)
SD <- c(ec$sd)
SampleSize <- c(ec$n)

gen_data <- function(means, sds, samplesizes){
  n.grp <- length(means)
  grps <- factor(rep(1:n.grp, samplesizes))
  dat <- lapply(1:n.grp, function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
  y <- do.call(rbind, dat)
  out <- data.frame(group = grps, y = y)
  out
}

simulated_data <- gen_data(Mean, SD,SampleSize)
av <- aov(y ~ group, data = simulated_data)
summary(av)
TukeyHSD(av)

#dg-ca3
dg<-read.csv("D:\\Brewer lab data\\HFS\\fb_IBSR_moment_tbl_DG-CA3.csv")

Mean <- c(dg$mean)
SD <- c(dg$sd)
SampleSize <- c(dg$n)

gen_data <- function(means, sds, samplesizes){
  n.grp <- length(means)
  grps <- factor(rep(1:n.grp, samplesizes))
  dat <- lapply(1:n.grp, function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
  y <- do.call(rbind, dat)
  out <- data.frame(group = grps, y = y)
  out
}

simulated_data <- gen_data(Mean, SD,SampleSize)
av <- aov(y ~ group, data = simulated_data)
summary(av)
TukeyHSD(av)

#ca3-ca1
ca3<-read.csv("D:\\Brewer lab data\\HFS\\fb_IBSR_moment_tbl_CA3-CA1.csv")

Mean <- c(ca3$mean)
SD <- c(ca3$sd)
SampleSize <- c(ca3$n)

gen_data <- function(means, sds, samplesizes){
  n.grp <- length(means)
  grps <- factor(rep(1:n.grp, samplesizes))
  dat <- lapply(1:n.grp, function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
  y <- do.call(rbind, dat)
  out <- data.frame(group = grps, y = y)
  out
}

simulated_data <- gen_data(Mean, SD,SampleSize)
av <- aov(y ~ group, data = simulated_data)
summary(av)
TukeyHSD(av)

#ca1-ec
ca1<-read.csv("D:\\Brewer lab data\\HFS\\fb_IBSR_moment_tbl_CA1-EC.csv")

Mean <- c(ca1$mean)
SD <- c(ca1$sd)
SampleSize <- c(ca1$n)

gen_data <- function(means, sds, samplesizes){
  n.grp <- length(means)
  grps <- factor(rep(1:n.grp, samplesizes))
  dat <- lapply(1:n.grp, function(i) {scale(rnorm(samplesizes[i]))*sds[i] + means[i]})
  y <- do.call(rbind, dat)
  out <- data.frame(group = grps, y = y)
  out
}

simulated_data <- gen_data(Mean, SD,SampleSize)
av <- aov(y ~ group, data = simulated_data)
summary(av)
TukeyHSD(av)
