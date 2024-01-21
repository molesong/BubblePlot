
library(readr)
library(stringr)
library(ggplot2)
library(tibble)
library(magrittr)
library(dplyr)
library(ggpubr)
library(reshape2)
library(cowplot)
library(gridExtra) #grid.arrange function
library(data.table)
library(tidyverse)

conflicts_prefer(cowplot::get_legend)
conflicts_prefer(reshape2::melt)


percentage_plot_fn = function(df,Altitude_){
  the_plot_DF <<- df %>% filter(Altitude == Altitude_)
  chainLengthSumDF = the_plot_DF %>% 
    group_by(FA.length) %>% 
    summarise(sum_same_length_intensity  = sum(median_intensity))
  percentage_array = the_plot_DF %>% 
    apply(.,1, function(x){
      x <<- x
      FA_length = x[which(names(x) == 'FA.length')]
      the_median_intensity = x[which(names(x) == 'median_intensity')] %>% as.numeric()
      sum_same_length_intensity = chainLengthSumDF %>% filter(FA.length == FA_length)  %>% mutate(sum_same_length_intensity = as.numeric(sum_same_length_intensity))
      result = the_median_intensity/sum_same_length_intensity[1,2]
      return(result)
    }) %>% unlist
  the_plot_DF$percentage = percentage_array
  the_plot_DF$FA.length = factor(the_plot_DF$FA.length)
  the_plot_DF <<- the_plot_DF
  bar_color <<- 'YlGnBu'
  p <- ggplot(the_plot_DF) +
    geom_col(aes(x = FA.length, y = percentage,fill = FA.bond))+
    scale_fill_distiller(name="bond number",palette=bar_color,direction = 1)+
    labs(title=Altitude_, x='FA Chain Length', y='Median Intensity') +
    theme_bw()+
    xlim(c('14','15','16','17','18','20','22'))+
    theme( axis.title=element_text(size=10,face="plain",color="black"),
           axis.text = element_text(size=10,face="plain",color="black"),
           panel.grid.minor = element_blank(),
           panel.grid.major = element_blank(),
           legend.background = element_blank(),
           plot.title = element_text(hjust = 0.5))+
    theme(legend.position = 'none')
  return(p)
}

theme_georgia <- function(...) {
  theme_gray(base_family = "Georgia", ...) + 
    theme(plot.title = element_text(face = "bold"))
}

bubble_plot_FA_fn = function(df,Altitude_){
  the_plot_DF <<-df  %>% filter(Altitude == Altitude_) %>% 
    mutate(FA.length = as.character(FA.length))
  FA_length_array = theQuanDF_combineFA$FA.length %>% table %>% names
  bubble_plot_color =  brewer.pal(length(FA_length_array), "Spectral")
  names(bubble_plot_color) = FA_length_array
  bubble_p= ggplot(the_plot_DF) +
    aes(x = FA.length, y = FA.bond, colour = FA.length,    size = median_intensity  ) +  #atttention: is colour, not color!
    geom_point(shape = "circle") +
    scale_color_manual(values=bubble_plot_color)+
    labs(title=Altitude_, x='FA Chain Length', y='FA Chain Bond') +
    theme_bw()+
    theme( axis.title=element_text(size=10,face="plain",color="black"),
           axis.text = element_text(size=10,face="plain",color="black"),
           legend.background = element_blank(),
           plot.title = element_text(hjust = 0.5))+
    theme(legend.position = 'none')+
    scale_x_discrete(limits=theQuanDF_combineFA$FA.length %>% table %>% names) #get the FA length array,and set x axis
  return(bubble_p)
}

bubble_plot_CL_fn = function(df,Altitude_){
  the_plot_DF <<- df  %>% filter(Altitude == Altitude_) %>% 
    mutate(CL_total_length = as.character(CL_total_length))
  CL_length_array = theQuanDF_combineCL$CL_total_length %>% table %>% names
  bubble_plot_color =  brewer.pal(length(CL_length_array), "Spectral")
  names(bubble_plot_color) = CL_length_array
  bubble_plot_cl= ggplot(the_plot_DF) +
    aes(x = CL_total_length , y = CL_total_bond , colour = CL_total_length,    size = percentage  ) +  #atttention: is colour, not color!
    geom_point(shape = "circle") +
    scale_color_manual(values=bubble_plot_color)+
    labs(title=Altitude_, x='FA Chain Length', y='FA Chain Bond') +
    theme_bw()+
    ylim(c(4,15))+
    theme( axis.title=element_text(size=10,face="plain",color="black"),
           axis.text = element_text(size=10,face="plain",color="black"),
           legend.background = element_blank(),
           plot.title = element_text(hjust = 0.5))+
    theme(legend.position = 'none')+
    scale_x_discrete(limits=CL_length_array) #get the FA length array,and set x axis
  return(bubble_plot_cl)
}
#read file


#Mingjun's project.
path = './input' 
files_array = list.files(path) %>% paste0(path,'/',.)
input_df = openxlsx::read.xlsx(files_array[1],sheet = 1)

input_df =input_df %>% 
  filter(grepl('\\d{2}',Species)) %>% 
  filter(!grepl('carnitine',Species)) %>% 
  mutate(rowname = Species) %>% 
  column_to_rownames(var = 'rowname')
  
target_lipid_class = c('PC','PE','PI','GluCer','SM','LacCer','GM3','CL')
target_lipid_class_char =  target_lipid_class%>% paste(., collapse = '\\b|\\b') %>% paste('\\b',.,'\\b',sep = '')

filtered_class_df = input_df %>% 
  mutate(lipid_class=str_extract(Species,'GM3|S1P|^[[:alpha:]]+\\W[[:alpha:]]{2,9}|^[[:alpha:]]+')) %>%  #a regular expression for selecting lipid class,
  filter(grepl(target_lipid_class_char,lipid_class))

sample_group_df = input_df %>% names %>% as.data.frame %>% 
  rename('sample_name' = ".") %>% 
  filter(sample_name != 'Species') %>% 
  separate(col = 'sample_name', into = c('group',NA), sep = '-',remove = F) %>% 
  filter(group != 'QC')
  
sample_group_array = sample_group_df$group %>% unique

median_group_level_DF = sample_group_array %>% lapply(., function(sample_group_name){  #calculate the median value of each lipid in each sample group.
  sample_group_name <<- sample_group_name
  df_ <<- filtered_class_df %>% select(contains(sample_group_name)) %>% 
    mutate(median_ = pmap_dbl(., ~ median(c(...)))) %>%   #calculate median value.
    select(median_)
}) %>% do.call(cbind,.) 
names(median_group_level_DF) = sample_group_array
median_group_level_DF$Species = rownames(median_group_level_DF)
median_group_level_DF$lipid_class = filtered_class_df$lipid_class
median_group_level_DF = median_group_level_DF %>% 
  mutate(lipid_class = filtered_class_df$lipid_class) %>% 
  relocate(Species,lipid_class)

median_group_level_DF_ = median_group_level_DF %>%  
  rowwise %>% 
  mutate(FA1 = case_when(!grepl('[_/]',Species) ~ Species,
                          grepl('[_/]',Species) ~ str_split(Species,pattern = '[_/]')[[1]][1])) %>% 
  mutate(FA2 = case_when(!grepl('[_/]',Species) ~ NA,
                         grepl('[_/]',Species) ~ str_split(Species,pattern = '[_/]')[[1]][2]))
  
# aa = median_group_level_DF_ %>% filter(lipid_class == 'CL') %>% select(Species) %>%unlist %>%  lapply(., function(x){
#   str_split(x, '\\(')[[1]][1] 
# }) %>% do.call(rbind,.) %>% as.array()   #for test if exist non-unique CL




# 创建一个包含文本的向量
text <- c("This is a PC", "There are LPC here", "The PC is not here", "The PS is not here")

# 使用grep函数进行全词匹配
matches <- grep("\\bPC\\b|\\bPS\\b", text, value=TRUE)

# 输出匹配的结果
print(matches)



#
theLipidClass = 'TAG'
SpeciesArrayOfTheLipidClass = lipid_class_df %>% filter(lipid_class == theLipidClass) %>% select(Species) %>% unlist
the_lipid_df = input_df %>% filter(Species %in% SpeciesArrayOfTheLipidClass) 

# dirs_array = dirs_array[-c(1,2)]
allDF =dirs_array %>% lapply(., function(x){        
  theOrganDirsFileList = list.files(x)
  result <<- grep('_cut95_conc_all.csv',theOrganDirsFileList,value = T)   %>% paste0(x,'/',.) #get the csv path
  result_ <<- read_csv(result)       #read the csv.
  names(result_)[1] = 'CL_FullName'
  return(result_)
}) 
group_fileName = './metaData/sample_group_information.csv'
group_df = read.csv(group_fileName,check.names = F)

#deal IDs using split function.
#get the four FAs on cardiolipin.
#get total chain length and total number of double bonds.
organ_name_array = dirs_array %>% str_match(.,'\\.\\/resultNew\\/(\\w+)') %>% .[,2]
1:length(organ_name_array) %>% lapply(.,function(index_organ){
  index_organ <<- index_organ
  organ_name = organ_name_array[index_organ]
  inputDF = allDF[[index_organ]]
 
  theQuanDF <<- allDF[[index_organ]] %>% melt(.,id.vars = 'CL_FullName') %>%    #second, deal quan table.
    filter(!is.na(value)) %>% 
    rename(sample_name = 'variable', intensity = 'value') %>% 
    rowwise() %>% 
    mutate(sample_number = str_extract(sample_name,'\\d{1,2}')) %>% 
    merge(group_df,.,by = 'sample_number') %>% 
    separate(col = 'CL_FullName',into = c('CL_ID','fourFA'),sep = '\\(',remove = F,fill = 'left') %>%
    separate(col = 'fourFA',into = c('FA1','FA2','FA3','FA4'),sep = '_',remove = T,fill = 'left') %>%
    separate(col = 'FA4',into = c('FA4',NA),sep = '\\)',remove = T,fill = 'left') %>% 
    select(-c(sample_number,sample_long_name,species,
              Gender,Body_Weight,Heart_Weight,Altitude_index,sample_name)) 
  
  theQuanDF_forFA= theQuanDF %>% 
    select(-c(CL_ID)) %>% 
    melt(.,id.vars = c('CL_FullName','Altitude','intensity')) %>%  #A wide matrix becomes a long matrix
    rename(FA_length_bond = 'value',FA_index = 'variable') %>% 
    group_split(.,Altitude) 
  
  theQuanDF_combineFA <<- lapply(theQuanDF_forFA,function(x){
    result = x %>% group_by(FA_length_bond) %>% 
          summarise(median_intensity=median(intensity),.groups = 'drop') %>%  #get the median of all CL's same FAs
          as.data.frame()  
    result$Altitude = x$Altitude[1]
    return(result)
    }) %>% do.call(rbind,.) %>% 
    rowwise() %>% 
    mutate(FA.length = as.numeric(str_extract(FA_length_bond,'\\d{2}')) )%>% 
    mutate(FA.bond = as.numeric(str_split(FA_length_bond,':')[[1]][2])) 
  
  
  #title diagram
  title_gg <<- ggdraw() +draw_label(organ_name,  x = 0.5, hjust = 0.5) #set title of the page
  
  #bar diagram 
  percentage_plot_list <<- theQuanDF_combineFA$Altitude %>% table %>% names %>% map(~ percentage_plot_fn(theQuanDF_combineFA,.x))
  p_forlegend <<- ggplot(the_plot_DF) +
    geom_col(aes(x = FA.length, y = percentage,fill = FA.bond))+
    scale_fill_distiller(name="bond number",palette=bar_color,direction = 1) 

  combined_plot <<- grid.arrange(grobs= percentage_plot_list,nrow = 2,right = get_legend(p_forlegend))  # combine all plot to one page 
  result_plot <<- plot_grid(title_gg,combined_plot,  ncol = 1, rel_heights = c(0.1, 1))
  ggsave(paste('./mls_result/barPlot/',organ_name,'.pdf',sep = ''),result_plot,width = 8,height = 5)
  
  
  #bubble diagram at FA level.
  bubble_plot_list <- theQuanDF_combineFA$Altitude %>% table %>% names %>% map(~ bubble_plot_FA_fn(theQuanDF_combineFA,.x))
  bubble_plot_forlegend = ggplot(the_plot_DF) +
    aes(x = FA.length, y = FA.bond, colour = FA.length,    size = median_intensity  ) +  #atttention: is colour, not color!
    geom_point(shape = "circle") +
    scale_color_brewer(palette = "Spectral", direction = 1)+
    theme_bw()+
    scale_x_discrete(limits=theQuanDF_combineFA$FA.length %>% table %>% names)+
    guides(size = 'none')  #only delete the legend of size(median_intensity)
  
  combined_bubble_plot <- grid.arrange(grobs= bubble_plot_list,nrow = 2,right = get_legend(bubble_plot_forlegend)) 
  result_bubble_plot <<- plot_grid(title_gg,combined_bubble_plot,  ncol = 1, rel_heights = c(0.1, 1))
  ggsave(paste('./mls_result/bubblePlot/',organ_name,'_FA_level.pdf',sep = ''),result_bubble_plot,width = 8,height = 6)
  
  
  #bubble diagran at CL level
  theQuanDF_combineCL <<-  theQuanDF %>% 
    select(-c(CL_FullName,FA1,FA2,FA3,FA4)) %>% 
    group_split(.,Altitude) %>% 
    lapply(.,function(x){
    result <<- x %>% group_by(CL_ID) %>% 
      summarise(median_intensity=median(intensity),.groups = 'drop') %>%  #get the median of all CL's same FAs
      as.data.frame()  
    result$Altitude = x$Altitude[1]
    sum_median_intensity = result$median_intensity %>% sum
    result = result %>% rowwise %>% mutate(percentage = median_intensity/sum_median_intensity) %>% select(-median_intensity)
    return(result)
  }) %>% do.call(rbind,.) %>% 
    rowwise() %>% 
    separate(col='CL_ID',
             sep = ':',
             into = c('CL_total_length','CL_total_bond'), 
             fill   = "right",
             remove = F) %>% 
    rowwise() %>% 
    mutate(CL_total_length = sub(x=CL_total_length, pattern = 'CL',replacement = '')) %>% 
    mutate(CL_total_length = as.numeric(CL_total_length)) %>% 
    mutate(CL_total_bond = as.numeric(CL_total_bond))
  
  ################

  ###################
  bubble_CL_plot_list <- theQuanDF_combineFA$Altitude %>% table %>% names %>% map(~ bubble_plot_CL_fn(theQuanDF_combineCL,.x))
  combined_bubble_CL_plot <- grid.arrange(grobs= bubble_CL_plot_list,nrow = 2) 
  result_bubble_CL_plot <<- plot_grid(title_gg,combined_bubble_CL_plot,  ncol = 1, rel_heights = c(0.1, 1))
  ggsave(paste('./mls_result/bubblePlot/',organ_name,'_CL_level.pdf',sep = ''),result_bubble_CL_plot,width = 8,height = 6)
  
  })


# 1:length(organ_name_array) %>% lapply(.,function(index_organ){
#   index_organ <<- index_organ
#   organ_name = organ_name_array[index_organ]
#   input_df = allDF[[index_organ]]
#   CL_ID_SplitDF <-  input_df %>% separate(col=CL_FullName,
#                                           sep = '\\(',
#                                           into = c('CL_total_chain_length','meiyong'), 
#                                           fill   = "right",
#                                           remove = TRUE) %>% 
#     select(-meiyong) %>% melt %>% 
#     mutate(sample_number = str_extract(sample_name,'\\d{1,2}')) %>% 
#     merge(group_df,.,by = 'sample_number') 
# })






FA_length_array = theQuanDF_combineFA$FA.length %>% table %>% names
df1 <- df <- data.frame(x=LETTERS[1:5],y=1+rnorm(5))
p3 <- ggplot(df,aes(x,y)) + geom_point() + scale_x_discrete(limits=c("E","D","C","B","A","F"))

percentage_array = the_plot_DF %>% 
  apply(.,1, function(x){
    x <<- x
    FA_length = x[which(names(x) == 'FA.length')]
    the_median_intensity = x[which(names(x) == 'median_intensity')] %>% as.numeric()
    sum_same_length_intensity = chainLengthSumDF %>% filter(FA.length == FA_length)  %>% mutate(sum_same_length_intensity = as.numeric(sum_same_length_intensity))
    result = the_median_intensity/sum_same_length_intensity[1,2]
    return(result)
  }) %>% unlist
the_plot_DF$percentage = percentage_array
the_plot_DF$FA.length = factor(the_plot_DF$FA.length)
the_plot_DF <<- the_plot_DF
bar_color <<- 'Spectral'

ggplot(the_plot_DF) +
  aes(x = FA.length, y = FA.bond, colour = FA.length,    size = median_intensity  ) +  #atttention: is colour, not color!
  geom_point(shape = "circle") +
  scale_color_brewer(palette = "Spectral", direction = 1)+
  labs(title=Altitude_, x='FA Chain Length', y='FA Chain Bond') +
  theme_bw()+
  theme( axis.title=element_text(size=10,face="plain",color="black"),
         axis.text = element_text(size=10,face="plain",color="black"),
         legend.background = element_blank(),
         plot.title = element_text(hjust = 0.5))+
  theme(legend.position = 'none')

return(p)



ggplot(the_plot_DF) +
  aes(
    x = FA.length,
    y = FA.bond,
    colour = FA.length,
    size = median_intensity
  ) +
  geom_point(shape = "circle") +
  scale_color_brewer(palette = "Spectral", direction = 1) +
  theme_bw()







# 
# 
# multi.page = egg::ggarrange(plot_list = list(p,p),
#           # labels = c("A", "B", "C"),
#           ncol = 3, nrow = 2, common.legend = TRUE)
# 



          # nrow = 2)

# p1 = ggplot(theQuanDF_combineFA) +
#   geom_point(aes(x = FA.length, y=FA.bond,size = median_intensity,col = FA.length ))+
#   theme_bw()+
#   theme(legend.position="none")+
#   # scale_colour_gradientn(colours = terrain.colors(10))
#   # scale_colour_gradient2()
#   scale_color_distiller(palette = "Spectral")
# 
# theQuanDF_combineFA$FA.length = factor(theQuanDF_combineFA$FA.length,levels = theQuanDF_combineFA$FA.length %>% table %>% names)
# ggscatterhist(
#   theQuanDF_combineFA, x = "FA.length", y = "FA.bond",  
#   shape=21,color ="black",fill= "median_intensity", size ='median_intensity', alpha = 1,
#   # palette = c("#00AFBB", "#E7B800", "#FC4E07","black",'green','grey'),
#   # margin.plot =  "boxplot",
#   margin.params = list(fill = "Altitude", color = "black", size = 0.2),
#   legend = c(0.82,0.15),
#   ggtheme = theme)
# aa =theQuanDF_forFA %>% 
#   select(CL_FullName, FA1, FA2,FA3,FA4) %>% 
#   melt(.,id.vars = 'CL_FullName') %>% 
#   merge()





# CL_length_and_bond_statistic = CL_ID_SplitDF$CL_ID %>% table %>% as.data.frame() %>% mutate(chain_length = str_replace(`.`, 'CL(\\d{2})\\:\\d{1,2}','\\1')%>% as.numeric()) %>% 
#   mutate(bond_number = str_replace(`.`, 'CL\\d{2}\\:(\\d{1,2})','\\1') %>% as.numeric()) 
# 
# p1 = ggplot(CL_length_and_bond_statistic) +
#   geom_point(aes(x = chain_length, y=bond_number,size = Freq,col = chain_length ))+
#   scale_x_continuous(breaks = seq(60,80,by = 4)) +
#   scale_y_continuous(breaks = seq(3,15,by = 2)) +
#   theme_bw()+
#   theme(legend.position="none")

# p2 = ggplot()


# ggscatterhist(aa , x = 'chain_length', y='bond_number')

# str_sub('CL68:5','CL(\\d{2})','\\1')
# str_match('CL68:5', 'CL(\\d{2})')[2]
# grep()



#make bar plot for exhibiting bond and chain length
fileName = './resultNew/all_organ_CL_result.xlsx'
df_ = xlsx::read.xlsx(fileName,sheetIndex = 1)  #read quantification information
names(df_)[1] = 'CL_name'
df_t = df_ %>% t %>% as.data.frame
names(df_t) = df_t[1,]; df_t %<>%  .[-1,]
df_t = df_t %>% tibble::rownames_to_column(var="sample_name") #please remember following code,I've looked it up a million times!
df_t = df_t %>% rowwise %>% mutate(organ = stringr::str_split(sample_name,'\\d')[[1]][1]) %>%   #get organ class in each row.
  mutate(sample_number = stringr::str_split(sample_name,'[:alpha:]')[[1]][2]) %>%
  as.data.frame() %>% 
  relocate(sample_name,organ,sample_number) %>% 
  merge(group_df, df_t, by = 'sample_number')

#loop species and organ to calculate intensity median of each altitude.
calculate_median = function(x){
  if((sum(is.na(x))/length(x))<=0.5){
    return(median(x %>% as.numeric(), na.rm=TRUE))
  }else{return(NA)}}

make_FA_length_or_bond_fn <- function(df_t,species_,organ_,length_or_bond){
  meiyong = df_t %>% filter(organ == organ_) %>%   #select organ and species
    select(c(-sample_number,-sample_long_name,-species,
             -Gender, -Body_Weight,-Heart_Weight,-sample_name,-organ,-Altitude))
  meiyong_ =  meiyong %>% t %>% as.data.frame() 
  colnames(meiyong_) = meiyong_[1,]
  meiyong_= meiyong_[-1,] 
  meiyong_$'CL_name' = rownames(meiyong_)
  meiyong_ = meiyong_ %>% relocate(CL_name)
  
  meiyong_long = melt(meiyong_,id.vars='CL_name')   #turn wide data into long data
  
  chain_and_bond_DF = meiyong_long[,1] %>% lapply(.,function(x){
    xx = str_replace(x,pattern = 'CL\\d{2}\\:\\d{1,2}\\(','')  #delete the first one
    xx = str_replace(xx,pattern = '\\)','')
    result_ = c(x,str_split(xx,'_')[[1]], str_split(xx,'_|:')[[1]]) #split  : and _
  }) %>% do.call(rbind,.) %>% as.data.frame()
  names(chain_and_bond_DF) = c('CL_name','FA1','FA2','FA3','FA4', 'FA1_length', 'FA1_bond','FA2_length', 'FA2_bond', 'FA3_length', 'FA3_bond','FA4_length', 'FA4_bond') 
  
  meiyong_long_ = cbind(meiyong_long, chain_and_bond_DF[,-1]) %>% rename(.,altitude_index =variable) %>%
    rowwise() %>% mutate(altitude = str_split(altitude_index,'_')[[1]][1]) %>% 
    select(-altitude_index) %>% relocate(CL_name, altitude)
  
  FA_length_or_bond_intensity_Df = meiyong_long_ %>% 
    filter(!is.na(value)) %>% as_tibble() %>%  group_split(altitude) %>%  #loop altitude
    lapply(.,function(x){
      x <<- x
      result  <<- x %>% select(c(CL_name,value,length,bond)) %>% 
        apply(.,1,function(z){
          return( data.frame(c(z[3],z[2]),c(z[4],z[2]),c(z[5],z[2]),c(z[6],z[2])) )
        }) %>% do.call(cbind,.) %>% t %>% as.data.frame() %>% rename('FA' = paste('FA1_',length_or_bond,sep = '')) %>%
        group_split(.,FA) %>% 
        lapply(.,function(y){
          FA_name <-  y[1,1] %>% as.character
          sum_ <- y %>% .$value %>% as.numeric() %>% sum
          return(c(FA_name, sum_))
        })  %>% 
        do.call(rbind, .) %>% as.data.frame() %>% rename('FA_name' = 'V1', 'intensity' = 'V2') %>% 
        mutate(intensity = as.numeric(intensity)) %>% 
        mutate(intensity = intensity / sum(intensity))     # 把intensity变成百分比
      
      result$altitude =  as.character(x[1,2])
      
      return(result)
    }) %>% do.call(rbind,.)
  
  match_list = list('H'='Heart', 'K' = 'Kidney', 'L' = 'Liver', 'S'='spleen')
  p = ggplot(FA_length_or_bond_intensity_Df) +
    geom_point(aes(x = FA_name,y=intensity,color= altitude)) +
    xlab(paste('FA',length_or_bond))+
    ylab(paste('Percentage'))+
    ggtitle(paste(species_, '   ',match_list[organ_] %>% as.character()))+
    theme_classic()
  
  ggsave(filename = paste(species_,'__',organ_,'__',length_or_bond,'.jpg',sep = ''),plot = p, width = 4,height = 3)
  
}


for (species_ in unique(df_t $species)){
  for(organ_ in unique(df_t$organ)){
    for(length_or_bond in c('length','bond')){
      make_FA_length_or_bond_fn(df_t,species_,organ_,length_or_bond)
    }
  }
}
species_ = 'Cansus zokor'
organ_ = 'L'
length_or_bond = 'length'

calculateLengthOrBondPercentage <- function(x,length_or_bond){
    length_result  <<- x %>% select(c(CL_name,value,contains(length_or_bond))) %>% 
      apply(.,1,function(z){
        return( data.frame(c(z[3],z[2]),c(z[4],z[2]),c(z[5],z[2]),c(z[6],z[2])) )
      }) %>% 
      do.call(cbind,.) %>% 
      t %>% as.data.frame() %>% 
      rename('FA' = paste('FA1_',length_or_bond,sep = '')) %>%  #Attention: the FA1 is NOT only FA1 column, it is called FA1. But it have four chain information.
      group_split(.,FA) %>% 
      lapply(.,function(y){
        FA_name <-  y[1,1] %>% as.character
        sum_ <- y %>% .$value %>% as.numeric() %>% sum
        return(c(FA_name, sum_))
      })  %>% 
      do.call(rbind, .) %>% as.data.frame() %>% rename('FA_name' = 'V1', 'intensity' = 'V2') %>% 
      mutate(intensity = as.numeric(intensity)) %>% 
      mutate(intensity = intensity / sum(intensity))     # 把intensity变成百分比
    result$altitude =  as.character(x[1,2])
    
    return(result)
  }
make_FA_length_or_bond_fn <- function(df_t,species_,organ_,length_or_bond){
  meiyong = df_t %>% filter(organ == organ_,species == species_) %>%   #select organ and species
    select(c(-sample_number,-sample_long_name,-species,-Gender, 
             -Body_Weight,-Heart_Weight,-sample_name,-organ,-Altitude))
  meiyong_ =  meiyong %>% t %>% as.data.frame() 
  colnames(meiyong_) = meiyong_[1,]
  meiyong_= meiyong_[-1,] 
  meiyong_$'CL_name' = rownames(meiyong_)
  meiyong_ = meiyong_ %>% relocate(CL_name)
  
  meiyong_long = melt(meiyong_,id.vars='CL_name')   #turn wide data into long data
  
  chain_and_bond_DF = meiyong_long[,1] %>% lapply(.,function(x){
    xx = str_replace(x,pattern = 'CL\\d{2}\\:\\d{1,2}\\(','')  #delete the first one
    xx = str_replace(xx,pattern = '\\)','')
    result_ = c(x,str_split(xx,'_')[[1]], str_split(xx,'_|:')[[1]]) #split  : and _
  }) %>% do.call(rbind,.) %>% as.data.frame()
  names(chain_and_bond_DF) = c('CL_name','FA1','FA2','FA3','FA4', 'FA1_length', 'FA1_bond','FA2_length', 'FA2_bond', 'FA3_length', 'FA3_bond','FA4_length', 'FA4_bond') 
  
  meiyong_long_ = cbind(meiyong_long, chain_and_bond_DF[,-1]) %>% rename(.,altitude_index =variable) %>%
    rowwise() %>% mutate(altitude = str_split(altitude_index,'_')[[1]][1]) %>% 
    select(-altitude_index) %>% relocate(CL_name, altitude)
  
  FA_length_or_bond_intensity_Df = meiyong_long_ %>% 
    filter(!is.na(value)) %>% as_tibble() %>%  group_split(altitude) %>%  #loop altitude
    lapply(.,calculateLengthOrBondPercentage,length_or_bond = 'length') %>% do.call(rbind,.)
  
  p = ggplot(FA_length_or_bond_intensity_Df) +
    geom_point(aes(x = FA_name,y=intensity,color= altitude)) +
    xlab(paste('FA',length_or_bond))+
    ylab(paste('Percentage'))+
    theme_classic()
  
  ggsave(filename = paste(species_,'__',organ_,'__',length_or_bond,'.jpg',sep = ''),plot = p, width = 4,height = 3)
  
}


meiyong_median = meiyong %>%          #calculate median  
  group_split(Altitude) %>% 
  sapply(.,function(x){
    apply(x[,-1], 2, calculate_median)
  }) %>% t %>% as.data.frame
row_name_ = meiyong %>%  group_split(Altitude) %>% sapply(.,function(x){
  x[,1] %>% unique %>% as.character()})# %>% paste(.,'_median',sep = '')
row.names(meiyong_median) = row_name_

meiyong_median_pencentage  = meiyong_median %>% apply(.,1, function(x){
  return(x/sum(x,na.rm = T))
}) %>%  as.data.frame() %>% filter(!if_all(.fns = is.na)) %>%  .[complete.cases(.),] %>% 
  tibble::rownames_to_column('CL_name')

# meiyong_median_pencentage = meiyong_median_pencentage %>% tibble::rownames_to_column('altitude') %>% reshape2::melt()
# meiyong_median_pencentage = meiyong_median_pencentage %>% filter(!is.na(value)) %>% mutate(value = as.numeric(value))

meiyong_median_pencentage$ymin<-apply(meiyong_median_pencentage[,c(2,3)], 1, min)
meiyong_median_pencentage$ymax<-apply(meiyong_median_pencentage[,c(2,3)], 1, max)

meiyong_median_pencentage$ymin1<-meiyong_median_pencentage$ymin
meiyong_median_pencentage$ymin1[as.integer((meiyong_median_pencentage$`1335m`-meiyong_median_pencentage$`2795m`)>0)]=NA

meiyong_median_pencentage$ymax1<-meiyong_median_pencentage$ymax
meiyong_median_pencentage$ymax1[as.integer((meiyong_median_pencentage$`1335m`-meiyong_median_pencentage$`2795m`)>0)==0]=NA

meiyong_median_pencentage$ymin2<-meiyong_median_pencentage$ymin
meiyong_median_pencentage$ymin2[as.integer((meiyong_median_pencentage$`1335m`-meiyong_median_pencentage$`2795m`)<=0)==0]=NA

meiyong_median_pencentage$ymax2<-meiyong_median_pencentage$ymax
meiyong_median_pencentage$ymax2[as.integer((meiyong_median_pencentage$`1335m`-meiyong_median_pencentage$`2795m`)<=0)==0]=NA

meiyong_median_pencentage = meiyong_median_pencentage %>% mutate(ID = 1:n())
ggplot(meiyong_median_pencentage, aes(x =ID))+
  # geom_ribbon( aes(ymin=ymin, ymax=ymax),alpha=0.5,fill="white",color=NA)+
  geom_ribbon( aes(ymin=ymin1, ymax=ymax1),alpha=0.5,fill="#FF6B5E",color=NA)+#,fill = AMZN > AAPL
  geom_ribbon( aes(ymin=ymin2, ymax=ymax2),alpha=0.5,fill="#00B2F6",color=NA)+#,fill = AMZN > AAPL
  geom_line(aes(y=`1335m`,color="#FF6B5E"),size=0.75)+#color="black",
  geom_line(aes(y=`2795m`,color="#00B2F6"),size=0.75)+#color="black",
  # scale_x_date(date_labels = "%Y",date_breaks = "2 year")+
  xlab("CL name")+ 
  ylab("Pencentage")+
  scale_colour_manual(name = "Variable",
                      labels = c("1335m", "2795m"),
                      values = c("#FF6B5E", "#00B2F6"))+
  theme( axis.title=element_text(size=10,face="plain",color="black"),
         axis.text = element_text(size=10,face="plain",color="black"),
         legend.position = c(0.15,0.8),
         legend.background = element_blank())

#二维散点图与箱型图的组合图表
library(ggpubr)

theme<-theme_minimal()+theme(
  axis.title=element_text(size=14,face="plain",color="black"),
  axis.text = element_text(size=12,face="plain",color="black"),
  legend.text= element_text(size=12,face="plain",color="black"),
  legend.title=element_text(size=12,face="plain",color="black"),
  legend.background=element_rect(fill=NA,colour=NA)
)

ggscatterhist(
  iris, x = "Sepal.Length", y = "Sepal.Width",  #iris
  shape=21,color ="black",fill= "Species", size ='Petal.Width', alpha = 1,
  palette = c("#00AFBB", "#E7B800", "#FC4E07"),
  margin.plot =  "boxplot",
  margin.params = list(fill = "Species", color = "black", size = 0.2),
  legend = c(0.82,0.15),
  ggtheme = theme)


heart_df = allDF[[1]]
