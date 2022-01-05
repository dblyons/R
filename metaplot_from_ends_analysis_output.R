## plot 5' and 3' average methylation from a set of features


## define this function for later
simple.colmeans<-function(y){
          sapply (y, function(x){ 
          colMeans(x[,-1], na.rm=TRUE)})}

# name your output matrix files like this if you want this script to work
## col 1: genotype
## col 2: generation
## col 5: context
## col 7: end (fiveEnds or threeEnds, f comes before t)

# import your files
filenames<-list.files("~/my_ends_output_matrices/",
                      pattern='.*Ends$', full.names = TRUE)


h1 <- lapply(filenames, read.table)

#make an empty list 
resultlis<-list()


for (i in 1:length(filenames)) {  
  
  # populate the list with previously matched of loci
  # 
  
  resultlis[[i]]<-simple.colmeans(h1[[i]])
  
  
}

together<-do.call(cbind, resultlis)
together<-as.data.frame(together)
names(together)<-basename(filenames)
together<-data.frame(pos=seq(-2500,2450,50), together)
melt.together<-melt(together, id='pos')

melt_thing<-str_split_fixed(melt.together$variable, '\\.', 10)
meltbig<-data.frame(melt.together, end=melt_thing[,7], 
                    generation=melt_thing[,2], genotype=melt_thing[,1], context=melt_thing[,5])

meltbig<-data.frame(meltbig, genotypeBygeneration=paste(meltbig$genotype, 
                                                        meltbig$generation))

ggplot(subset(meltbig, generation=="f2_1"|generation=="f4_1"|
                generation=='f7_1'|generation=='f11_1' ), aes(pos, value, color=genotypeBygeneration))+
  geom_line(lwd=1)+
  theme_bw()+
  facet_grid(context~end, scales='free_y')+
  theme_bw(base_size = 22 )+#remove background
  
  theme(legend.title=element_blank())+
  theme( panel.margin = unit(0, "in"))+#squish panels together
  
  theme(axis.text=element_text(size=22))+ #
  theme(strip.text.x = element_blank())+#remove panel headers
  theme(strip.text.y = element_blank())+
  scale_x_continuous(expand = c(0, 0), limits=c(-2000,2000))+#set limits of x
  ylim(0,0.75)+
  geom_vline(xintercept = 0,  col='grey10', lty=2)+
  labs(y='', x='')+#add axis lab
  theme(legend.title = element_text(size=0, face="bold"))+
  scale_color_manual(values =
                       c( "orange3",
                          "steelblue",
                          "firebrick",
                          "forestgreen"
                       ))

