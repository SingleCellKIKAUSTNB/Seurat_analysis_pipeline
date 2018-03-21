#Blank chart
blankPlot <- ggplot()+geom_blank(aes(1,1)) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

#Visualize data for filtering propouse
seurat_data_plot<- function(output_pdf_name) {
  pdf(file=output_pdf_name, width=12, height=6)
  #
  for (n in 1:length(file_list)){
    vln=VlnPlot(
      object = list_all_Seurat[[n]], 
      features.plot = c("nGene", "nUMI"), 
      nCol = 2
    )
    scatterPlot = ggplot(list_all_Seurat[[n]]@meta.data, aes(x=nUMI, y=nGene)) +
      geom_point() + geom_rug()
    
    xdensity = ggplot(list_all_Seurat[[n]]@meta.data, aes(nUMI)) + 
      geom_density(alpha=.5) + 
      theme(legend.position = "none") +
      geom_vline(xintercept = mean(list_all_Seurat[[n]]@meta.data$nUMI), color='red') +
      geom_vline(xintercept = median(list_all_Seurat[[n]]@meta.data$nUMI), color='blue')
    
    ydensity = ggplot(list_all_Seurat[[n]]@meta.data, aes(nGene)) + 
      geom_density(alpha=.5) + 
      theme(legend.position = "none")+
      geom_vline(xintercept = mean(list_all_Seurat[[n]]@meta.data$nGene), color='red') +
      geom_vline(xintercept = median(list_all_Seurat[[n]]@meta.data$nGene), color='blue')
    
    scater_gen_umi = gridExtra::arrangeGrob(xdensity, blankPlot, scatterPlot, ydensity, 
                                  ncol=2, nrow=2, widths=c(4, 2), heights=c(1.4, 4))
    
    plot(ggarrange(vln, scater_gen_umi,
              labels = c("A", "B"),
              ncol = 2, nrow = 1,
              widths = c(1, 1.5))
         )
  }
  dev.off()
}

















