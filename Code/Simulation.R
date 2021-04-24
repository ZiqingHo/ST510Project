# Simulation 1
# 10 informative with all contaminated, 100 noise with 10 contaminated noise
# 10% contamination 
s <- seq(1.1,sqrt(110),0.5)
set.seed(8)
##################################
# K = 3
k<-2:7
result.K.3<-c()
n_iter <- 0
while(n_iter<50){
  size <- c(10,10,20)
  test<-SimData(K = 3, cluster_size = size,
                p_inf = 10, p_noise = 100,
                p_inf_out = 10, out_type = F, inf_out = 0.10,
                p_noise_out = 10, noise_out = 0.10,
                unif_interval = c(-12,-6,6,12),
                mu_interval = c(-6,-3,3,6),
                scatter_interval = c(3,10),
                rho_interval = c(0.1,0.9))
  x.scale<-scale(test$X)
  test.W <- matrix(0,nrow=length(k),ncol=length(s))
  for (i in 1:length(k)){
    test.B <- wrskB(x.scale, K = k[i], s = s)$W_sk/k[i]
    lb <- length(test.B)
    test.W[i,1:lb] <- test.B
  }
  colnames(test.W) <- paste0("s", 1:length(s))
  rownames(test.W) <- paste0("K = ", k)
  test.Ball<-ball_list(test.W)
  result.K.3<-c(result.K.3,ball_index(test.Ball,k=k))
  n_iter<-n_iter+1
}
Range.K3<-as.data.frame(result.K.3[1:50])
colnames(Range.K3)<-"x"
Bar_ball1<-ggplot(Range.K3, aes(x=x))+
  geom_bar(color="black", fill="white")+
  scale_x_continuous(lim=c(2,7))+
  ylim(0,50)+
  labs(y="Count")+
  ggtitle("K = 3") +
  theme(axis.title.y = element_blank(),axis.title.x=element_blank())


#################################


##################################
# K = 4
k<-2:7
result.K.4<-c()
n_iter <- 0
while(n_iter<50){
  size <- c(10,10,10,20)
  test<-SimData(K = 4, cluster_size = size,
                p_inf = 10, p_noise = 100,
                p_inf_out = 10, out_type = F, inf_out = 0.10,
                p_noise_out = 10, noise_out = 0.10,
                unif_interval = c(-12,-6,6,12),
                mu_interval = c(-6,-3,3,6),
                scatter_interval = c(3,10),
                rho_interval = c(0.1,0.9))
  x.scale<-scale(test$X)
  test.W <- matrix(0,nrow=length(k),ncol=length(s))
  for (i in 1:length(k)){
    test.B <- wrskB(x.scale, K = k[i], s = s)$W_sk/k[i]
    lb <- length(test.B)
    test.W[i,1:lb] <- test.B
  }
  test.Ball<-ball_list(test.W)
  result.K.4<-c(result.K.4,ball_index(test.Ball,k=k))
  n_iter<-n_iter+1
}

Range.K4<-as.data.frame(result.K.4[1:50])
colnames(Range.K4)<-"x"
Bar_ball2<-ggplot(Range.K4, aes(x=x))+
  geom_bar(color="black", fill="white")+
  scale_x_continuous(lim=c(2,7))+
  ylim(0,50)+
  ggtitle("K = 4")+
  labs(y="Count")+
  theme(axis.title.y = element_blank(),axis.title.x=element_blank())


#################################

##################################

# K = 10
k<-7:12
set.seed(10)
result.K.10<-c()
n_iter <- 0
while(n_iter<50){
  size <- c(rep(10,9),20)
  test<-SimData(K = 10, cluster_size = size,
                p_inf = 10, p_noise = 100,
                p_inf_out = 10, out_type = F, inf_out = 0.10,
                p_noise_out = 10, noise_out = 0.10,
                unif_interval = c(-12,-6,6,12),
                mu_interval = c(-6,-3,3,6),
                scatter_interval = c(3,10),
                rho_interval = c(0.1,0.9))
  x.scale<-scale(test$X)
  test.W <- matrix(0,nrow=length(k),ncol=length(s))
  for (i in 1:length(k)){
    test.B <- wrskB(x.scale, K = k[i], s = s)$W_sk/k[i]
    lb <- length(test.B)
    test.W[i,1:lb] <- test.B
  }
  test.Ball<-ball_list(test.W)
  result.K.10<-c(result.K.10,ball_index(test.Ball,k=k))
  n_iter<-n_iter+1
}
Range.K10<-as.data.frame(result.K.10)
colnames(Range.K10)<-"x"
Bar_ball10<-ggplot(Range.K10, aes(x=x))+
  geom_bar(color="black", fill="white")+
  scale_x_continuous(lim=c(7,12))+
  ylim(0,50)+
  labs(y="Count")+
  ggtitle("K = 10")+
  theme(axis.title.y = element_blank(),axis.title.x=element_blank())

#################################
# Figure 1
Bar_Ball<-ggarrange(Bar_ball1,Bar_ball2,Bar_ball10,nrow=3)
annotate_figure(Bar_Ball,bottom = "Number of clusters selected",left = "Count")

##################################
# Simulation 2
# Estimation of K and s for K = 3

set.seed(12)
s <- seq(1.1,sqrt(110),0.5)
k<-3:5
k.h<-2:5
result.G1.par<-c()
n_varW_G1<-c()
result.CH1<-c()
n_varW_CH1<-c()
result.H1<-c()
n_varW_H1<-c()
index.s.G1<-c()
index.s.CH1<-c()
index.s.H1<-c()
n_iter <- 0
while(n_iter<50){
  size <-  c(rep(10,2),20)
  test<-SimData(K = 3, cluster_size = size,
                p_inf = 10, p_noise = 100,
                p_inf_out = 10, out_type = F, inf_out = 0.10,
                p_noise_out = 10, noise_out = 0.10,
                unif_interval = c(-12,-6,6,12),
                mu_interval = c(-6,-3,3,6),
                scatter_interval = c(3,10),
                rho_interval = c(0.1,0.9))
  x.scale<-scale(test$X)
  # Gap
  
  gap.G1<-c()
  n_varW.G<-c()
  index.s.G<-c()
  
  # CH
  
  index.s.CH<-c()
  CH.G<-c()
  n_varW.CH<-c()
  
  # Hartigan
  
  n_varW.H <- matrix(NA,nrow = length(k.h),ncol = length(s))
  test.H <- matrix(0,nrow = length(k.h),ncol = length(s))
  
  for (i in 1:length(k)){
    #Gap
    test.G <- wrskGap_par(x.scale, K = k[i], s = s)
    index.max<-which.max(test.G$Gap)
    index.opt<-which(test.G$Gap[1:index.max]>(test.G$Gap[index.max]-test.G$se[index.max]))[1]
    n_varW.G<-c(n_varW,sum(test.G$result[[index.opt]]$varweights!=0))
    index.s.G<-c(index.s.G,index.opt)
    gap.G1<-c(gap.G1,test.G$Gap[index.opt])
    #CH
    test.CH <- wrskCH(x.scale, K = k[i], s = s)
    index.max.CH<-which.max(test.CH$CH)
    index.s.CH<-c(index.s.CH,index.max.CH)
    n_varW.CH<-c(n_varW.CH,sum(test.CH$result[[index.max.CH]]$varweights!=0))
    CH.G<-c(CH.G,test.CH$CH[index.max.CH])
  }
  
  #Gap
  id.max.gap<-which.max(gap.G1)
  n_varW_G1<-c(n_varW_G1,n_varW.G[id.max.gap])
  result.G1.par<-c(result.G1.par,k[id.max.gap])
  index.s.G1<-c(index.s.G1,s[id.max.gap])
  #CH
  id.max.ch<-which.max(CH.G)
  index.s.CH1<-c(index.s.CH1,s[id.max.ch])
  n_varW_CH1<-c(n_varW_CH1,n_varW.CH[id.max.ch])
  result.CH1<-c(result.CH1,k[id.max.ch])
  
  for (i in 1:length(k.h)){
    test <- wrskH(x.scale, K = k.h[i], s = s)
    lh <- length(test$W_sk)
    test.H[i,1:lh] <- test$W_sk
    for(j in 1:lh){
      n_varW.H[i,j]<-sum(test$result[[j]]$varweights!=0)
    }
  }
  # H
  test.hart<-Hart.sk(test.H,k.h)
  result.Hart<-Hart.sk.p2(test.hart,k.h)
  index.s.H1<-c(index.s.H1,s[result.Hart$s])
  result.H1<-c(result.H1,k.h[result.Hart$K])
  n_varW_H1<-c(n_varW_H1,n_varW.H[result.Hart$K,result.Hart$s])
  
  n_iter<-n_iter+1
}

optimal.K1<-majority.K(result.G1.par,result.CH1,result.H1)

##################################
# Estimation of K and s for K = 4
set.seed(12)
s <- seq(1.1,sqrt(110),0.5)
k<-3:5
k.h<-2:5
result.G10.par<-c()
n_varW_G10<-c()
result.CH10<-c()
n_varW_CH10<-c()
result.H10<-c()
n_varW_H10<-c()
index.s.G10<-c()
index.s.CH10<-c()
index.s.H10<-c()
n_iter <- 0
while(n_iter<50){
  size <-  c(rep(10,3),20)
  test<-SimData(K = 4, cluster_size = size,
                p_inf = 10, p_noise = 100,
                p_inf_out = 10, out_type = F, inf_out = 0.10,
                p_noise_out = 10, noise_out = 0.10,
                unif_interval = c(-12,-6,6,12),
                mu_interval = c(-6,-3,3,6),
                scatter_interval = c(3,10),
                rho_interval = c(0.1,0.9))
  x.scale<-scale(test$X)
  # Gap
  
  gap.G1<-c()
  n_varW.G<-c()
  index.s.G<-c()
  # CH
  
  index.s.CH<-c()
  CH.G<-c()
  n_varW.CH<-c()
  
  # Hartigan
  
  n_varW.H <- matrix(NA,nrow = length(k.h),ncol = length(s))
  test.H <- matrix(0,nrow = length(k.h),ncol = length(s))
  
  for (i in 1:length(k)){
    #Gap
    test.G <- wrskGap_par(x.scale, K = k[i], s = s)
    index.max<-which.max(test.G$Gap)
    index.opt<-which(test.G$Gap[1:index.max]>(test.G$Gap[index.max]-test.G$se[index.max]))[1]
    n_varW.G<-c(n_varW,sum(test.G$result[[index.opt]]$varweights!=0))
    index.s.G<-c(index.s.G,index.opt)
    gap.G1<-c(gap.G1,test.G$Gap[index.opt])
    #CH
    test.CH <- wrskCH(x.scale, K = k[i], s = s)
    index.max.CH<-which.max(test.CH$CH)
    index.s.CH<-c(index.s.CH,index.max.CH)
    n_varW.CH<-c(n_varW.CH,sum(test.CH$result[[index.max.CH]]$varweights!=0))
    CH.G<-c(CH.G,test.CH$CH[index.max.CH])
  }
  
  #Gap
  id.max.gap<-which.max(gap.G1)
  n_varW_G10<-c(n_varW_G1,n_varW.G[id.max.gap])
  result.G10.par<-c(result.G10.par,k[id.max.gap])
  index.s.G10<-c(index.s.G10,s[id.max.gap])
  #CH
  id.max.ch<-which.max(CH.G)
  index.s.CH10<-c(index.s.CH10,s[id.max.ch])
  n_varW_CH10<-c(n_varW_CH10,n_varW.CH[id.max.ch])
  result.CH10<-c(result.CH10,k[id.max.ch])
  
  #H
  for (i in 1:length(k.h)){
    test <- wrskH(x.scale, K = k.h[i], s = s)
    lh <- length(test$W_sk)
    test.H[i,1:lh] <- test$W_sk
    for(j in 1:lh){
      n_varW.H[i,j]<-sum(test$result[[j]]$varweights!=0)
    }
  }
  test.hart<-Hart.sk(test.H,k.h)
  result.Hart<-Hart.sk.p2(test.hart,k.h)
  index.s.H10<-c(index.s.H10,s[result.Hart$s])
  result.H10<-c(result.H10,k.h[result.Hart$K])
  n_varW_H10<-c(n_varW_H10,n_varW.H[result.Hart$K,result.Hart$s])
  
  n_iter<-n_iter+1
}

optimal.K2<-majority.K(result.G10.par,result.CH10,result.H10)

###########################
# Figure 2
n_K<-c(rep("3",50),rep("4",50))
combined.ch<-c(result.CH1,
               result.CH10)
combined.gap<-c(result.G1.par,result.G10.par)
combined.h<-c(result.H1,result.H10)
combined.m<-c(optimal.K1,optimal.K2)

data.sim2.ch<-data.frame(combined.ch,n_K)
data.sim2.gap<-data.frame(combined.gap,n_K)
data.sim2.h<-data.frame(combined.h,n_K)
data.sim2.m<-data.frame(combined.m,n_K)

bar.sim2.ch<-ggplot(data.sim2.ch,aes(x=combined.ch))+
  geom_bar(aes(fill=as.factor(n_K)),position="dodge",color="black")+
  scale_x_continuous(lim=c(1.5,5.5))+
  ylim(0,40)+
  ggtitle("CH")+
  labs(fill="K")+
  theme(axis.title.y = element_blank(),axis.title.x=element_blank())+
  scale_fill_manual(labels = c("3", "4"), values = c("royalblue1", "brown3")) 

bar.sim2.gap<-ggplot(data.sim2.gap,aes(x=combined.gap))+
  geom_bar(aes(fill=as.factor(n_K)),position="dodge",color="black")+
  scale_x_continuous(lim=c(1.5,5.5))+
  ylim(0,40)+
  ggtitle("Gap")+
  labs(fill="K")+
  theme(axis.title.y = element_blank(),axis.title.x=element_blank())+
  scale_fill_manual(labels = c("3", "4"), values = c("royalblue1", "brown3")) 

bar.sim2.h<-ggplot(data.sim2.h,aes(x=combined.h))+
  geom_bar(aes(fill=as.factor(n_K)),position="dodge",color="black")+
  scale_x_continuous(lim=c(1.5,5.5))+
  ylim(0,40)+
  ggtitle("Hartigan")+
  labs(fill="K")+
  theme(axis.title.y = element_blank(),axis.title.x=element_blank())+
  scale_fill_manual(labels = c("3", "4"), values = c("royalblue1", "brown3")) 

bar.sim2.m<-ggplot(data.sim2.m,aes(x=combined.m))+
  geom_bar(aes(fill=as.factor(n_K)),position="dodge",color="black")+
  scale_x_continuous(lim=c(1.5,5.5))+
  ylim(0,40)+
  ggtitle("Majority")+
  labs(fill="K")+
  theme(axis.title.y = element_blank(),axis.title.x=element_blank())+
  scale_fill_manual(labels = c("3", "4"), values = c("royalblue1", "brown3")) 

bar.Sim2.all<-ggarrange(bar.sim2.ch,bar.sim2.gap,bar.sim2.h,bar.sim2.m,nrow=2,ncol=2,
                        common.legend = T,legend="top")
annotate_figure(bar.Sim2.all,bottom = "Number of clusters selected",left = "Count")

###########################
# Simulation 3

# 2%
set.seed(12)
s.opt=9.6
cer.wrsk<-c()
cer.sk<-c()
cer.rskc<-c()
cer.kmeans<-c()
n_iter = 0
while(n_iter<50){
  # size <- round(runif(4,15,150))
  size<-c(50,50,100,150)
  
  dat<-SimData(K = 4, cluster_size = size,
               p_inf = 80, p_noise = 3920,
               p_inf_out = 20, out_type = T, inf_out = 0.20,
               p_noise_out = 790, noise_out = 0.10,
               unif_interval = c(-12,-6,6,12),
               mu_interval = c(-6,-3,3,6),
               scatter_interval = c(3,10),
               rho_interval = c(0.1,0.9))
  x.scale <- scale(dat$X)
  y.true <- dat$y
  
  ### optimal s = 9.6
  # test.G <- wrskGap_par(x.scale, K = 4, s = s3[1:30])
  # index.max<-which.max(test.G$Gap)
  # index.opt<-which(test.G$Gap[1:index.max]>(test.G$Gap[index.max]-test.G$se[index.max]))[1]
  
  # WRSK
  result.wrsk <- WRSK(x.scale, K = 4, s = s.opt)
  y.wrsk<- result.wrsk$clusters
  #CER(y.wrsk,y.true)
  #sum(result.wrsk$varweights!=0)
  cer.wrsk<-c(cer.wrsk,CER(y.wrsk,y.true))
  
  # standard Kmeans
  result.kc <- kmeans(x.scale, centers=4)
  y.kc <- result.kc$cluster  
  cer.kmeans <- c(cer.kmeans, CER(y.kc,y.true))
  
  # SC
  result.sk <- KMeansSparseCluster(x.scale, K = 4, wbounds = s.opt)
  y.sk <- result.sk[[1]]$Cs
  cer.sk <- c(cer.sk, CER(y.sk,y.true))
  
  # RSKC
  result.rskc <- RSKC(x.scale, ncl = 4, alpha = 0.3, L1 = s.opt)
  y.rskc <- result.rskc$labels 
  cer.rskc <- c(cer.rskc, CER(y.rskc,y.true))
  
  n_iter <- n_iter + 1 
}

# 3%
set.seed(11)
s.opt=9.6
cer.wrsk4<-c()
cer.sk4<-c()
cer.rskc4<-c()
cer.kmeans4<-c()
n_iter = 0
while(n_iter<50){
  # size <- round(runif(4,15,150))
  size<-c(50,50,100,150)
  
  dat<-SimData(K = 4, cluster_size = size,
               p_inf = 120, p_noise = 3880,
               p_inf_out = 30, out_type = T, inf_out = 0.20,
               p_noise_out = 785, noise_out = 0.10,
               unif_interval = c(-12,-6,6,12),
               mu_interval = c(-6,-3,3,6),
               scatter_interval = c(3,10),
               rho_interval = c(0.1,0.9))
  x.scale <- scale(dat$X)
  y.true <- dat$y
  
  ### optimal s = 9.6
  # test.G <- wrskGap_par(x.scale, K = 4, s = s3[1:30])
  #  index.max<-which.max(test.G$Gap)
  # index.opt<-which(test.G$Gap[1:index.max]>(test.G$Gap[index.max]-test.G$se[index.max]))[1]
  
  # WRSK
  result.wrsk <- WRSK(x.scale, K = 4, s = s.opt)
  y.wrsk<- result.wrsk$clusters
  #CER(y.wrsk,y.true)
  #sum(result.wrsk$varweights!=0)
  cer.wrsk4<-c(cer.wrsk4,CER(y.wrsk,y.true))
  
  # standard Kmeans
  result.kc <- kmeans(x.scale, centers=4)
  y.kc <- result.kc$cluster  
  cer.kmeans4 <- c(cer.kmeans4, CER(y.kc,y.true))
  
  # SC
  result.sk <- KMeansSparseCluster(x.scale, K = 4, wbounds = s.opt)
  y.sk <- result.sk[[1]]$Cs
  cer.sk4 <- c(cer.sk4, CER(y.sk,y.true))
  
  # RSKC
  result.rskc <- RSKC(x.scale, ncl = 4, alpha = 0.3, L1 = s.opt)
  y.rskc <- result.rskc$labels 
  cer.rskc4 <- c(cer.rskc4, CER(y.rskc,y.true))
  
  n_iter <- n_iter + 1 
}

# 5%
set.seed(10)
s.opt=12.1
cer.wrsk2<-c()
cer.sk2<-c()
cer.rskc2<-c()
cer.kmeans2<-c()
n_iter = 0

while(n_iter<50){
  size<-c(50,50,100,150)
  
  dat<-SimData(K = 4, cluster_size = size,
               p_inf = 200, p_noise = 3800,
               p_inf_out = 50, out_type = T, inf_out = 0.20,
               p_noise_out = 770, noise_out = 0.10,
               unif_interval = c(-12,-6,6,12),
               mu_interval = c(-6,-3,3,6),
               scatter_interval = c(3,10),
               rho_interval = c(0.1,0.9))
  x.scale <- scale(dat$X)
  y.true <- dat$y
  
  ### optimal s = 12.1
  # test.G <- wrskGap_par(x.scale, K = 4, s = s3[1:30])
  #  index.max<-which.max(test.G$Gap)
  #  index.opt<-which(test.G$Gap[1:index.max]>(test.G$Gap[index.max]-test.G$se[index.max]))[1]
  
  # WRSK
  result.wrsk <- WRSK(x.scale, K = 4, s = s.opt)
  y.wrsk<- result.wrsk$clusters
  #CER(y.wrsk,y.true)
  #sum(result.wrsk$varweights!=0)
  cer.wrsk2<-c(cer.wrsk2,CER(y.wrsk,y.true))
  
  # standard Kmeans
  result.kc <- kmeans(x.scale, centers = 4)
  y.kc <- result.kc$cluster  
  cer.kmeans2 <- c(cer.kmeans2, CER(y.kc,y.true))
  
  # SC
  result.sk <- KMeansSparseCluster(x.scale, K = 4, wbounds = s.opt)
  y.sk <- result.sk[[1]]$Cs
  cer.sk2 <- c(cer.sk2, CER(y.sk,y.true))
  
  # RSKC
  result.rskc <- RSKC(x.scale, ncl = 4, alpha = 0.3, L1 = s.opt)
  y.rskc <- result.rskc$labels 
  cer.rskc2 <- c(cer.rskc2, CER(y.rskc,y.true))
  
  n_iter <- n_iter + 1 
}

# Figure 4
labels.f4<-c(rep("2%",200),rep("3%",200),rep("5%",200))
measure.f4<-rep(c(rep("WRSK",50),rep("KC",50),rep("SK",50),rep("RSKC",50)),3)
cer.f4<-c(cer.wrsk,cer.kmeans,cer.sk,cer.rskc,
          cer.wrsk4,cer.kmeans4,cer.sk4,cer.rskc4,cer.wrsk2,cer.kmeans2,
          cer.sk2,cer.rskc2)
#cer.wrsk3,cer.kmeans3,cer.sk3,cer.rskc3)
data.f4<-data.frame(labels.f4,measure.f4,cer.f4)

colnames(data.f4)<-c("Label","Method","CER")
data.f4$Label<-as.factor(data.f4$Label)
ggplot(data = data.f4, aes(x=Method, y=CER)) +
  geom_boxplot(aes(fill=Label))+
  scale_fill_manual(labels = c("2 %", "3 %","5 %"), 
                    values = c("royalblue1", "brown3","green4")) +
  labs(fill="Proportion of informative variables")+
  ylim(0,0.75)+
  ylab("Classification error rate")+
  theme(axis.title.x = element_blank(),legend.position = "top")

