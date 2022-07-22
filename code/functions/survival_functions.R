KM_plot <- function(metadata, fit_model, col, title, ylab, legendlabs){
  kmplot <- ggsurvplot(fit_model, data = metadata,
                       palette = as.vector(col), conf.int=FALSE, 
                       xlim = c(0,5000), break.x.by = 1000,
                       #Risk table
                       risk.table = TRUE,
                       # pvalue
                       pval = TRUE, pval.size = 10, pval.coord = c(200, 0.1),
                       # legend
                       legend.title="", font.legend = 25, legend.labs = legendlabs, legend = c(0.75, 0.9),
                       # fonts
                       font.main = 30, font.x = 30, font.y = 30, font.tickslab = 30, 
                       # titles
                       title = title, xlab = "Time (days)", ylab = ylab)
  
  
  
  kmplot$table <- ggrisktable(fit_model, data = metadata, 
                              color = "strata", palette = as.vector(col),
                              fontsize = 10, risk.table.title = "",
                              xlim = c(0,5000),break.time.by = 1000,
                              y.text = TRUE, ylab = "",  xlab = "",legend.labs = legendlabs,
                              tables.theme = theme_cleantable(), font.tickslab = 20)
  
  return(kmplot)
}


#Stacked barplot related to Fig2A-B
stacked_plots <- function(freqtab, var1, var2){
  
  myplot <- ggplot(freqtab, aes(fill = eval(as.name(var2)), y = freq, x = eval(as.name(var1)))) + 
    geom_bar(position="stack", stat="identity") +
    #   scale_fill_brewer(palette = color) + 
    theme(axis.title.y = element_text(size = 35),
          axis.title.x = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text.x = element_text(size = 35,angle = 45, hjust = 1, color = "black"),
          axis.text.y = element_text(size = 35, color = "black")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "transparent", colour = NA),
          plot.title = element_text(size = 35)) +
    theme(legend.position = "right", legend.direction="vertical",
          legend.text = element_text(size = 35),
          legend.key.height= unit(2, 'cm'),
          legend.key.width= unit(2, 'cm'),
          legend.title = element_blank())
  
  return(myplot)
}