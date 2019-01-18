require("readxl")
require("dplyr")
require("ggplot2")


#####CHAIN LENGTH STANDARD MIXTURES USED#####

#####Set working directory#####
working.directory <- "C:/folder1/folder2/folder3/"

#####PREPARE DATASET FOR PATTERN RECONSTRUCTION#####
{
#Open input file
input <- read_xlsx(paste(working.directory,"input_ref.xlsx", sep = ""))
input$STD_code <- as.factor(input$STD_code)


#Create all possible binary combinations between given calibration sets
Combinations <- combn(x = levels(input$STD_code),m = 2,FUN = NULL,simplify = TRUE)

#Store sum RFs for each chain length standard
input <- input %>% group_by(Reference_standard) %>% mutate(Sum_response_factor = sum(Response_factor, na.rm = TRUE))
input[c(1:5)] <- lapply(input[c(1:5)], as.factor) 
input$Response_factor[is.na(input$Response_factor)] <- 0
}


#####SECTION FOR SAMPLE PROCESSING#####


####Open sample data####
{
sample <- read_xlsx(paste(working.directory,"sample.xlsx", sep = ""))
####Set sample name####
sample.name <- "Generic example"
####RUN PATTERN RECONSTRACTION FOR SELECTED (LOADED) SAMPLE####
{
sample$Chain_length <- as.factor(sample$Chain_length)
sample$Area <- as.numeric(sample$Area)
#Calculate relative area distribution within each homologue group
sample$Relative_distribution <- NA
sample$Area[is.na(sample$Area)] <- 0
sample <- sample %>% group_by(Chain_length) %>% mutate(Relative_distribution = Area/sum(Area, na.rm = TRUE))
results <- sample
results[c("Comp_1","Comp_2","Fraction_Comp_1","Simulated_pattern")] <- NA
#Deconvolation of homologue patterns

for (i in 1:length(levels(sample$Chain_length))) {
REF <- sample$Relative_distribution[sample$Chain_length == levels(input$Chain_length)[i]]
Distance <- 100
for (z in 1:length(Combinations[1,])){                                     
C_1 <- subset(input,subset = (STD_code == Combinations[1,z] & Chain_length == levels(sample$Chain_length)[i]))
C_2 <- subset(input,subset = (STD_code == Combinations[2,z] & Chain_length == levels(sample$Chain_length)[i]))
for (j in 1:100) {
  Combo <- (C_1$Response_factor*j+C_2$Response_factor*(100-j))/sum((C_1$Response_factor*j+C_2$Response_factor*(100-j)), na.rm=TRUE)
  if (Distance > (sum(sqrt((REF-Combo)^2)))) {
    results$Comp_1[results$Chain_length == levels(sample$Chain_length)[i]] <- as.character(C_1$STD_code)
    results$Comp_2[results$Chain_length == levels(sample$Chain_length)[i]] <- as.character(C_2$STD_code)
    results$Fraction_Comp_1 <- j
    results$Simulated_pattern[results$Chain_length == levels(input$Chain_length)[i]] <- Combo
    Distance <- sum(sqrt((REF-Combo)^2))
  }
}
}
}

#Calculate concentrations (ng per microliter)

results[c("RF_1st", "RF_2nd", "Concentration")] <- NA
results[c("RF_1st", "RF_2nd", "Concentration")] <- sapply(results[c("RF_1st", "RF_2nd", "Concentration")], as.numeric)
results$Fraction_Comp_1<- as.numeric(results$Fraction_Comp_1)
results$RF_1st[results$Chain_length == levels(results$Chain_length)[1]] <- input$Sum_response_factor[input$STD_code == results$Comp_1 & input$Chain_length == levels(results$Chain_length)[1]]

for (i in 1:length(levels(results$Chain_length))) {
  results$RF_1st[results$Chain_length == levels(results$Chain_length)[i]] <- input$Response_factor[input$STD_code == results$Comp_1 & input$Chain_length == levels(results$Chain_length)[i]]
  results$RF_2nd[results$Chain_length == levels(results$Chain_length)[i]] <- input$Response_factor[input$STD_code == results$Comp_2 & input$Chain_length == levels(results$Chain_length)[i]]
  }

results <- results %>% group_by(Chain_length) %>% mutate(Concentration= sum(Area)/sum(RF_1st*(Fraction_Comp_1/100)+RF_2nd*((100-Fraction_Comp_1)/100)))

#Visualization of results

plot_table<-data.frame(Distribution = c(results$Relative_distribution,results$Simulated_pattern),Homologue = results$Homologue, Chain_length = results$Chain_length, Origin = rep(as.factor(c("Measured","Simulated")), each = nrow(results)))
plot_table$Homologue <- factor(plot_table$Homologue, levels=unique(plot_table$Homologue))

plot <- ggplot(plot_table, aes(x = Homologue,y = Distribution*100, fill = Origin, colour = Origin))+
  geom_bar(stat="identity",position = position_dodge(width = 0.9), width = 0.8, size = .4)+
  theme(panel.background = element_blank())+
  scale_fill_manual(values=c("darkolivegreen3", "darkslategray4"))+
  scale_color_manual(values=c("darkolivegreen4", "darkslategray"))+
  ggtitle(label = paste(sample.name," - Distribution of CP homologues")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0))+
  xlab("") + ylab("Relative area distribution, %")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.key.size =  unit(0.15, "in"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10),
        panel.background = element_rect(fill = "white"),panel.border = element_rect(fill = NA, colour = "grey20"),
        panel.grid.major.y = element_line(colour = "grey50"),
        panel.grid.minor.y = element_line(colour = "grey80"),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",strip.background = element_rect(fill = "grey20"), strip.text = element_text(colour = "white", face = "bold"))+
  facet_wrap(.~ Chain_length, scales = "free",nrow = 4, ncol = 4)


results_output <-results %>% group_by(Chain_length)  %>% summarise(median(Concentration))
results_output$Type <- results$Type[1]
results_output$Sample <- sample.name
results_output$Chain_length <- as.character(results_output$Chain_length)
results_output <- cbind(results_output,results %>% group_by(Chain_length) %>% summarise(Comp_1 = first(Comp_1), Comp_2 = first(Comp_2), Fraction_Comp_1 = first(Fraction_Comp_1)))
results_output$Comment <- paste("The best match:",results_output$Fraction_Comp_1[1],"% of ",results_output$Comp_1[1]," and ",100-results_output$Fraction_Comp_1[1],"% of ",results_output$Comp_2[1], sep="")
results_output <- rbind(results_output,c(paste("Sum ",results_output$Type[1],"s",sep=""),sum(results_output$`median(Concentration)`),results_output$Type[1],sample.name))
results_output <- results_output[-c(5:8)]
colnames(results_output) <- c("Chain length","Concentration, ng/\U00B5L","Type","Sample name","Comment")
}
}
####VIEW RESULTS####
#View Overview
print(results_output)
#View full results
print(results)
#View graph
plot

####SAVE RESULTS####
#Set save directory
save.directory <- "C:/folder1/folder2/folder3/"

#Save Overview
write.table(results_output, paste(save.directory,sample.name,"-overview", ".txt", sep = ""), sep="\t") 
#Save full results
write.table(results, paste(save.directory,sample.name,"-full_results", ".txt", sep = ""), sep="\t") 
#Save distribution plot
ggsave(filename = paste(sample.name,"-CP_distribution_plot",".tiff", sep =""),device = "tiff",plot = plot, path = save.directory,width = unit(10,"cm"),height = unit(5,"cm"))



#References
#Hadley Wickham and Jennifer Bryan (2018). readxl: Read Excel Files. R package version 1.1.0. URL:https://CRAN.R-project.org/package=readxl
#Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2018). dplyr: A Grammar of Data Manipulation. R package version 0.7.6. URL:https://CRAN.R-project.org/package=dplyr
#H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
