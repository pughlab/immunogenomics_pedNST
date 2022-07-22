col = c("Missense" = "#d53e4f", 
        "In frame Ins/Del" = "#66bd63",
        "Nonsense" = "#f46d43",
        "Frameshift Ins/Del" = "#2166ac",
        "Nonstop" = "#bf812d", 
        "Fusion" = "black")

alter_fun = list(
  background = function(x, y, w, h){
    grid.rect(x, y, w, h, 
              gp = gpar(fill = "#c8cacf", col = "white"))},
  "Missense" = function(x, y, w, h){
    grid.rect(x, y, w*0.9, h*0.9, 
              gp = gpar(fill = col["Missense"], 
                        col = NA))},
  "In frame Ins/Del" = function(x, y, w, h){ 
    grid.rect(x, y, w*0.9, h*0.9, 
              gp = gpar(fill = col["In frame Ins/Del"], 
                        col = NA))},
  "Nonsense" = function(x,y,w,h){
    grid.rect(x,y,w*0.9,h*0.9, 
              gp = gpar(fill = col["Nonsense"], 
                        col = NA))}, 
  "Frameshift Ins/Del" = function(x,y,w,h){
    grid.rect(x,y,w*0.9,h*0.9, 
              gp = gpar(fill = col["Frameshift Ins/Del"], 
                        col = NA))},
  "Fusion" = function(x,y,w,h){
    grid.rect(x,y,w*0.9,h*0.4, 
              gp = gpar(fill = col["Fusion"], 
                        col = NA))}
)

IC_oncoprint <- function(IC, dta, mat,showrownames){
  myIC <- mat[,dta$sample_id[dta$immune_cluster == IC]]
  myIC[is.na(myIC)] <- ""
  # number of samples in each cluster
  n <- ncol(myIC)
  myIC_oncoprint = oncoPrint(myIC,
                             alter_fun = alter_fun,
                             col = col, 
                             show_heatmap_legend = FALSE,
                             show_column_names = FALSE,
                             show_row_names = showrownames,
                             row_names_gp = gpar(fontsize = 50),
                             remove_empty_rows = FALSE,
                             show_pct = FALSE,
                             right_annotation = NULL, 
                             row_names_side = "left",
                             top_annotation = NULL,
                             row_order = rownames(mat),
                             column_title = paste0(IC,"\n","(n = ", n, ")"), column_title_gp = gpar(fontsize = 50)) 
  
  return(myIC_oncoprint)
}


