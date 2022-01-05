##
## use in conjunction with kallisto output, which generates 1 directory per sample, each with identical contents
## rename directories A thru whatever then adjust names here accordingly once imported to R
##

library(sleuth)

base_dir<- "~/kallisto_output/"

samples <- paste0("sample", c("A","B", "C","D","E","F","G","H","I","J","K"))

kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))

## these names need to be adjusted accordingly
s2c <- data.frame(path=kal_dirs, sample=samples,
                  sample = c( "wt.1" ,   "wt.2",
                              "h1.1","h1.2","h1.3",
                              "ddm1.1","ddm1.2","ddm1.3",
                              "h1ddm1.1","h1ddm1.2","h1ddm1.3"
                             ),
                  stringsAsFactors=FALSE)

# get replicate information from directory names
s2c_thing<-str_split_fixed(s2c$sample.1, '\\.', 2)
s2c<-data.frame(s2c, s2c_thing)
names(s2c)[4]="genotype"
names(s2c)[5]="replicate"


## isolate the pair of genotypes to compare
s2c<-s2c[s2c$genotype=="ddm1" | s2c$genotype=='wt',]

# run the program
so<-sleuth_prep(s2c, ~genotype)
so <- sleuth_fit(so, ~genotype, "full")
so <- sleuth_fit(so, ~1, "reduced")
so <- sleuth_lrt(so, 'reduced', 'full')
lrt_results <- sleuth_results(so, 'reduced:full', test_type = 'lrt')

tmp_1<- ifelse(grepl ("[A|B]", so$obs_norm_filt$sample), 'wt',  'ddm1')

a_q<-lrt_results[,c(1:4)]

a<-data.frame(so$obs_norm_filt, tmp_1)
names(a)[7]="genotype"

a_q_test<-merge(a, a_q, by="target_id")

head(a_q_test)
kOut1<-a_q_test[,c(1,3:10)]
ko.agg<-aggregate(kOut1$tpm, list(gt=kOut1$genotype, id=kOut1$target_id, qval=kOut1$qval), mean)

ko.new<-dcast(ko.agg, qval+id ~ gt, value.var="x")

ko.new=transform(ko.new, fold.ch=((ddm1+1)/(wt+1)))

ko.new=transform(ko.new, mean=(ddm1+wt)/2)

names(ko.new)[7]='feature'

###subset if desired to remove rRNA snoRNA etc
ko.new.clean<-ko.new[ko.new$feature!='rRNA' & ko.new$feature!='snoRNA' & ko.new$feature!='snRNA',]



#plotting kallisto output for h1ddm1 and ddm1
# 'clean' indicates removed rRNA snoRNA snRNA as above
ggplot (subset(ko.new.clean, feature="Note=protein_coding_gene"), aes(log2(wt+1),log2(h1+1)))+
  #facet_wrap(~feature)+
  geom_point(aes((log2(wt+1)),log2(h1+1),
                 col=qval<0.05),
             size=6, alpha=I(1/3) )+
  
  scale_colour_manual(name = '',
                      values = setNames(c('red','black'),c(T, F)))+
  #ylim(0,15)+xlim(0,15)+
  #geom_abline(slope=1, intercept = 0, col='dodgerblue', lwd=1.5, alpha=I(1/3))+
  
  geom_abline(slope=1, intercept = 0, col='dodgerblue', lwd=1.5, alpha=I(1/3))+
  geom_abline(slope=1, intercept = -1, col='orange', lwd=1, lty=2, alpha=I(1/2))+
  geom_abline(slope=1, intercept = 1, col='orange', lwd=1, lty=2, alpha=I(1/2))+
  theme(axis.text=element_text(size=22))+
  theme_bw(base_size = 22 )+
  #ylim(0,15)+
  #xlim(0,15)+
  labs(x='WT TPM',
       y='h1 TPM')
