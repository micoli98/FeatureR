###1. Breakpoint count per 5MB
chrlen <- read.table(system.file("extdata", "hg38.chrom.sizes.txt", package = "FeatureR"),
                     sep="\t",stringsAsFactors = F)[1:24,]
centromeres <- read.table(system.file("extdata", "hg38.centromeres.txt", package = "FeatureR"),
                     sep="\t",stringsAsFactors = F)

#Extraction of the 5MB segments and counts per each one
bp5MB_extract <- function(samples, segmentation, chromosome_lengths = chrlen)
{
  out<-c()
  samps<- samples$sample
  for(i in samps)
  {
    segTab<-segmentation %>%
      filter(sample == i)

    chrs<-unique(segTab$chromosome)
    allBPnum<-c()
    for(c in chrs)
    {
      currseg<-segTab[segTab$chromosome==c,]
      intervals<-seq(1,chrlen[chrlen[,1]==c,2]+25000000,25000000)
      res <- hist(as.numeric(currseg$end[-nrow(currseg)]),breaks=intervals,plot=FALSE)$counts
      allBPnum<-c(allBPnum,res)
    }
    out<-rbind(out,cbind(sample=rep(i,length(allBPnum)),value=allBPnum))
  }
  rownames(out)<-NULL
  out <- data.frame(out,stringsAsFactors = F)
  out$value <- as.numeric(out$value)

  return(out)
}

#Discretization of the counts in classes
bp5MB_discretize <- function(df, seed, plot)
{
  #discretization of the dataset through Poisson mixture models
  if(!is.na(seed))
  {
    set.seed(seed)
  }

  dat <- df$value
  control<-new("FLXcontrol")
  control@minprior<-0.001
  control@iter.max<-1000

  fit<-stepFlexmix(dat ~ 1,model = flexmix::FLXMCmvpois(),k=1:10,nrep=1,control=control,verbose=F)
  fit<-getModel(fit,which="BIC")

  posteriors <- flexmix::posterior(fit) %>% as.data.frame()

  #Calculate the model of the fit (Poisson)
  metrics_df <- data.frame(value = df$value, group = fit@cluster)
  metrics_df <- metrics_df %>%
    group_by(group) %>%
    summarize(mean=mean(value), min=min(value), max=max(value), .groups = "drop") %>%
    arrange(mean) %>%
    mutate(line = (min - lag(max, default = 0))/2 + max) %>%
    as.data.frame()
  metrics_df$group <- paste0(rep("fiveMB", nrow(metrics_df)), seq(1, nrow(metrics_df)))

  #create plotting elements
  lines <- log1p(metrics_df$line)
  lines <- c(0.4, lines[2:(length(lines)-1)])
  metrics_df <- metrics_df %>% dplyr::select(-min, -max, -line)

  #Understand which group every observation belongs to and boundaries of the groups:
  colnames(posteriors) <- paste0(rep("fiveMB", ncol(posteriors)), seq(1, ncol(posteriors)))

  posteriors$group <- colnames(posteriors)[apply(posteriors,1,which.max)]

  df$group <- posteriors$group
  discr <- df %>%
    group_by(sample, group) %>%
    summarise(count = n(), .groups = "drop") %>%
    as.data.frame() %>%
    spread(group, count, fill = 0)

  results <- c(list(metrics_df), list(discr))

  if(plot)
  {
    df$log <- log1p(df$value)

    p <- ggplot(df, aes(x=log)) +
      labs(title = paste0("bp5Mb ", nrow(metrics_df), " components"),
           x = "log") +
      #ggtitle(paste0(names_features[f])) +
      geom_density(color="#7f7f7f", fill="#a5a5a5") +
      theme_classic() +
      theme(plot.title = element_text(hjust=0.5, size = 10), axis.title.y = element_blank(),
            axis.title.x = element_text(color = "#666666", size = 8), axis.text.x = element_text(size=8))

    for (x_val in lines) {
      p <- p + geom_vline(xintercept = x_val, linetype = "longdash", color = "red")
    }

    print(p)
  }

  return(results)
}

###2. Oscillation CN chains
oscil_extract <- function(samples, segmentation, chromosome_lengths = chrlen)
{
  samps <- samples$sample
  out<-c()
  for(i in samps)
  {
    segTab<-segmentation %>%
      filter(sample == i)

    chrs<-unique(segTab$chromosome)
    oscCounts<-c()
    for(c in chrs)
    {
      currseg<-segTab[segTab$chromosome==c,"copyNumber"]
      currseg<-round(as.numeric(currseg))
      if(length(currseg)>3)
      {
        prevval<-currseg[1]
        count=0
        for(j in 3:length(currseg))
        {
          if(currseg[j]==prevval&currseg[j]!=currseg[j-1])
          {
            count<-count+1
          }else{
            oscCounts<-c(oscCounts,count)
            count=0
          }
          prevval<-currseg[j-1]
        }
      }
    }
    out<-rbind(out,cbind(sample=rep(i,length(oscCounts)),value=oscCounts))
    if(length(oscCounts)==0)
    {
      out<-rbind(out,cbind(ID=i,value=0))
    }
  }
  rownames(out)<-NULL
  out <- data.frame(out,stringsAsFactors = F)
  out$value <- as.numeric(out$value)

  return(out)
}

###3. Breakpoint count per arm
bpArm_extract <- function(samples, segmentation, chromosome_lengths = chrlen, centromere_positions = centromeres)
{
  out<-c()
  samps<-samples$sample
  for(i in samps)
  {

    segTab<-segmentation %>%
      filter(sample == i) %>%
      dplyr::select(chromosome, start, end, copyNumber)

    chrs<-unique(segTab$chromosome)
    all_dists<-c()
    for(c in chrs)
    {
      if(nrow(segTab)>1)
      {
        starts<-as.numeric(segTab[segTab$chromosome==c,2])[-1]
        segstart<-as.numeric(segTab[segTab$chromosome==c,2])[1]
        ends<-as.numeric(segTab[segTab$chromosome==c,3])
        segend<-ends[length(ends)]
        ends<-ends[-length(ends)]
        centstart<-as.numeric(centromeres[centromeres$chr == c, 'start']) #centromeres[substr(centromeres[,2],4,5)==c,3])
        centend<-as.numeric(centromeres[centromeres$chr == c, 'end'])
        #centend<-as.numeric(centromeres[substr(centromeres[,2],4,5)==c,4])
        chrend<-chrlen[chrlen[,1]==c,2]
        ndist<-cbind(rep(NA,length(starts)),rep(NA,length(starts)))
        ndist[starts<=centstart,1]<-(centstart-starts[starts<=centstart])/(centstart-segstart)*-1
        ndist[starts>=centend,1]<-(starts[starts>=centend]-centend)/(segend-centend)
        ndist[ends<=centstart,2]<-(centstart-ends[ends<=centstart])/(centstart-segstart)*-1
        ndist[ends>=centend,2]<-(ends[ends>=centend]-centend)/(segend-centend)
        ndist<-apply(ndist,1,min)
        ndist<-na.omit(ndist)
        all_dists<-rbind(all_dists,sum(ndist>0))
        all_dists<-rbind(all_dists,sum(ndist<=0))
      }
    }
    if(nrow(all_dists)>0)
    {
      out<-rbind(out,cbind(sample=i,value=all_dists[,1]))
    }
  }
  rownames(out)<-NULL
  out <- data.frame(out,stringsAsFactors = F)
  out$value <- as.numeric(out$value)

  return(out)
}

###4. Changepoint
changp_extract <- function(samples, segmentation)
{
  out<-c()
  samps<-samples$sample
  for(i in samps)
  {
    segTab <- segmentation %>%
      filter(sample == i) %>%
      dplyr::select(chromosome, start, end, copyNumber)

    segTab$copyNumber[as.numeric(segTab$copyNumber)<0]<-0
    chrs<-unique(segTab$chromosome)
    allcp<-c()
    for(c in chrs)
    {
      currseg<-as.numeric(segTab[segTab$chromosome==c,"copyNumber"])
      allcp<-c(allcp,abs(currseg[-1]-currseg[-length(currseg)]))
    }
    if(length(allcp)==0)
    {
      allcp<-0 #if there are no changepoints
    }
    out<-rbind(out,cbind(sample=rep(i,length(allcp)),value=allcp))
  }
  rownames(out)<-NULL
  out <- data.frame(out,stringsAsFactors = F)
  out$value <- as.numeric(out$value)

  return(out)
}

