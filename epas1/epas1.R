library(R.utils)
library(matrixStats)
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))


# Custom Epas1 transcript extended by 5-RACE
Epas1_sequence <- unlist(strsplit(readLines(file("Epas1_long.fasta", open = "r"), n=2)[[2]], ""))

liver  <- vector("list", length=length(list.files(path="./liver_long", pattern='coverage')))
kidney <- vector("list", length=length(list.files(path="./kidney_long", pattern='coverage')))

names(liver) <- list.files(path="./liver_long", pattern='coverage')
names(kidney) <- list.files(path="./kidney_long", pattern='coverage')

for(i in 1:length(liver)) {
  handle <- file(paste0(getwd(), "/liver_long/", names(liver)[i]), open = "r");
  line <- readLines(handle, n=1)
  line <- unlist(strsplit(line, "\t"))
  close(handle)
  liver[[i]] <- tail(line, n = -1)
}
for(i in 1:length(kidney)) {
  handle <- file(paste0(getwd(), "/kidney_long/", names(kidney)[i]), open = "r");
  line <- readLines(handle, n=1)
  line <- unlist(strsplit(line, "\t"))
  close(handle)
  kidney[[i]] <- tail(line, n = -1)
}

liver_table  <- read.table("./liver_long/liver_table.txt", header = TRUE, stringsAsFactors = FALSE)
kidney_table <- read.table("./kidney_long/kidney_table.txt", header = TRUE, stringsAsFactors = FALSE)

liver_structure <- vector("list", length = length(unique(liver_table[,2]))) 
names(liver_structure) <- as.character(sort(unique(liver_table[,2])))
for(i in names(liver_structure)) {    
  liver_structure[[i]] <-  vector("list", length = nrow(liver_table[liver_table[2] == i,]))
  names(liver_structure[[i]]) <- liver_table[liver_table[2] == i, 1]
  #fill structure with data
  for(j in names(liver_structure[[i]])) {   
    liver_structure[[i]][[j]] <- as.numeric(liver[[j]])
  }
}

kidney_structure <- vector("list", length = length(unique(kidney_table[,2]))) 
names(kidney_structure) <- as.character(sort(unique(kidney_table[,2])))
for(i in names(kidney_structure)) {    
  kidney_structure[[i]] <-  vector("list", length = nrow(kidney_table[kidney_table[2] == i,]))
  names(kidney_structure[[i]]) <- kidney_table[kidney_table[2] == i, 1]
  #fill structure with data
  for(j in names(kidney_structure[[i]])) {   
    kidney_structure[[i]][[j]] <- as.numeric(kidney[[j]])
  }
}

#normalize replicates by total number of reads
liver_data1 <- lapply(liver_structure, function(a) { lapply(a, function(c) { output <- c / sum(c, na.rm = TRUE); return(output)  })})
kidney_data1 <- lapply(kidney_structure, function(a) { lapply(a, function(c) { output <- c / sum(c, na.rm = TRUE); return(output)  })})
liver_data <- lapply(liver_data1, function(b) { if(length(b) > 1) {rowMeans(simplify2array(b))} else {return(as.numeric(unlist(b)))}  })  
kidney_data <- lapply(kidney_data1, function(b) { if(length(b) > 1) {rowMeans(simplify2array(b))} else {return(as.numeric(unlist(b)))}  })  

plot_data <-    list( "liver" = data.frame(cbind("position" = c(1:4000),  as.data.frame(liver_data,  col.names=names(liver_data))[1:4000,]   )),
                      "kidney" = data.frame(cbind("position" = c(1:4000), as.data.frame(kidney_data, col.names=names(kidney_data))[1:4000,]  )),
                      "sequence" = Epas1_sequence
)


library(shiny)
library(plotly)
library(DT)
library(webshot)

    server <- function(input, output) {

      plotInput <- reactive({
        p1 <- plot_ly(plot_data[["liver"]], hoverinfo='none') %>%
             add_trace(x = ~position, y = ~X0, name = 'coverage', type = "bar", marker = list(color="#9ecae1"), showlegend = FALSE) %>%
             layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
                     xaxis = list(range=c(-nrow(plot_data[["liver"]])*0.05, nrow(plot_data[["kidney"]])),
                                  title = 'nucleotide position',
                                  autorange=FALSE,
                                  ticks='outside',
                                  tickwidth=2),
                 bargap = 0,
                 barmode = "overlay",
                 plot_bgcolor = 'rgba(235,235,235,1)',
                 legend = list(orientation = "h", bgcolor = "#08519c"),
                 margin = list(l = 0, r = 0, t = 0, b = 0)
          )
        p2 <- plot_ly(plot_data[["liver"]], hoverinfo='none') %>%
          add_trace(x = ~position, y = ~X15, name = 'coverage', type = "bar", marker = list(color="#6baed6"), showlegend = FALSE) %>%
          layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
                  xaxis = list(range=c(-nrow(plot_data[["liver"]])*0.05, nrow(plot_data[["kidney"]])),
                               title = 'nucleotide position',
                               autorange=FALSE,
                               ticks='outside',
                               tickwidth=2),
                  bargap = 0,
                  barmode = "overlay",
                  plot_bgcolor = 'rgba(235,235,235,1)',
                  legend = list(orientation = "h", bgcolor = "#08519c"),
                  margin = list(l = 0, r = 0, t = 0, b = 0)
          )
        p3 <- plot_ly(plot_data[["liver"]], hoverinfo='none') %>%
          add_trace(x = ~position, y = ~X30, name = 'coverage', type = "bar", marker = list(color="#4292c6"), showlegend = FALSE) %>%
          layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
                  xaxis = list(range=c(-nrow(plot_data[["liver"]])*0.05, nrow(plot_data[["kidney"]])),
                               title = 'nucleotide position',
                               autorange=FALSE,
                               ticks='outside',
                               tickwidth=2),
                  bargap = 0,
                  barmode = "overlay",
                  plot_bgcolor = 'rgba(235,235,235,1)',
                  legend = list(orientation = "h", bgcolor = "#08519c"),
                  margin = list(l = 0, r = 0, t = 0, b = 0)
          )
        p4 <- plot_ly(plot_data[["liver"]], hoverinfo='none') %>%
          add_trace(x = ~position, y = ~X45, name = 'coverage', type = "bar", marker = list(color="#2171b5"), showlegend = FALSE) %>%
          layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
                  xaxis = list(range=c(-nrow(plot_data[["liver"]])*0.05, nrow(plot_data[["kidney"]])),
                               title = 'nucleotide position',
                               autorange=FALSE,
                               ticks='outside',
                               tickwidth=2),
                  bargap = 0,
                  barmode = "overlay",
                  plot_bgcolor = 'rgba(235,235,235,1)',
                  legend = list(orientation = "h", bgcolor = "#08519c"),
                  margin = list(l = 0, r = 0, t = 0, b = 0)
          )
        p5 <- plot_ly(plot_data[["liver"]], hoverinfo='none') %>%
          add_trace(x = ~position, y = ~X300, name = 'coverage', type = "bar", marker = list(color="#084594"), showlegend = FALSE) %>%
          layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
                  xaxis = list(range=c(-nrow(plot_data[["liver"]])*0.05, nrow(plot_data[["kidney"]])),
                               title = 'nucleotide position',
                               autorange=FALSE,
                               ticks='outside',
                               tickwidth=2),
                  bargap = 0,
                  barmode = "overlay",
                  plot_bgcolor = 'rgba(235,235,235,1)',
                  legend = list(orientation = "h", bgcolor = "#08519c"),
                  margin = list(l = 0, r = 0, t = 0, b = 40)
          )
        p6 <- plot_ly(plot_data[["kidney"]], hoverinfo='none') %>%
          add_trace(x = ~position, y = ~X0, name = 'coverage', type = "bar", marker = list(color="#fc9272"), showlegend = FALSE) %>%
          layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
                  xaxis = list(range=c(-nrow(plot_data[["liver"]])*0.05, nrow(plot_data[["kidney"]])),
                               title = 'nucleotide position',
                               autorange=FALSE,
                               ticks='outside',
                               tickwidth=2),
                  bargap = 0,
                  barmode = "overlay",
                  plot_bgcolor = 'rgba(235,235,235,1)',
                  legend = list(orientation = "h", bgcolor = "#08519c"),
                  margin = list(l = 0, r = 0, t = 0, b = 0)
          )
        p7 <- plot_ly(plot_data[["kidney"]], hoverinfo='none') %>%
          add_trace(x = ~position, y = ~X15, name = 'coverage', type = "bar", marker = list(color="#fb6a4a"), showlegend = FALSE) %>%
          layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
                  xaxis = list(range=c(-nrow(plot_data[["liver"]])*0.05, nrow(plot_data[["kidney"]])),
                               title = 'nucleotide position',
                               autorange=FALSE,
                               ticks='outside',
                               tickwidth=2),
                  bargap = 0,
                  barmode = "overlay",
                  plot_bgcolor = 'rgba(235,235,235,1)',
                  legend = list(orientation = "h", bgcolor = "#08519c"),
                  margin = list(l = 0, r = 0, t = 0, b = 0)
          )
        p8 <- plot_ly(plot_data[["kidney"]], hoverinfo='none') %>%
          add_trace(x = ~position, y = ~X30, name = 'coverage', type = "bar", marker = list(color="#ef3b2c"), showlegend = FALSE) %>%
          layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
                  xaxis = list(range=c(-nrow(plot_data[["liver"]])*0.05, nrow(plot_data[["kidney"]])),
                               title = 'nucleotide position',
                               autorange=FALSE,
                               ticks='outside',
                               tickwidth=2),
                  bargap = 0,
                  barmode = "overlay",
                  plot_bgcolor = 'rgba(235,235,235,1)',
                  legend = list(orientation = "h", bgcolor = "#08519c"),
                  margin = list(l = 0, r = 0, t = 0, b = 0)
          )
        p9 <- plot_ly(plot_data[["kidney"]], hoverinfo='none') %>%
          add_trace(x = ~position, y = ~X45, name = 'coverage', type = "bar", marker = list(color="#cb181d"), showlegend = FALSE) %>%
          layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
                  xaxis = list(range=c(-nrow(plot_data[["liver"]])*0.05, nrow(plot_data[["kidney"]])),
                               title = 'nucleotide position',
                               autorange=FALSE,
                               ticks='outside',
                               tickwidth=2),
                  bargap = 0,
                  barmode = "overlay",
                  plot_bgcolor = 'rgba(235,235,235,1)',
                  legend = list(orientation = "h", bgcolor = "#08519c"),
                  margin = list(l = 0, r = 0, t = 0, b = 0)
          )
        p10 <- plot_ly(plot_data[["kidney"]], hoverinfo='none') %>%
          add_trace(x = ~position, y = ~X300, name = 'coverage', type = "bar", marker = list(color="#99000d"), showlegend = FALSE) %>%
          layout( yaxis = list(title = 'read count', autorange = TRUE, showgrid = FALSE, ticks="", showticklabels = FALSE),
                  xaxis = list(range=c(-nrow(plot_data[["liver"]])*0.05, nrow(plot_data[["kidney"]])),
                               title = 'nucleotide position',
                               autorange=FALSE,
                               ticks='outside',
                               tickwidth=2),
                  bargap = 0,
                  barmode = "overlay",
                  plot_bgcolor = 'rgba(235,235,235,1)',
                  legend = list(orientation = "h", bgcolor = "#08519c"),
                  margin = list(l = 0, r = 0, t = 0, b = 40)
          )
        p <- subplot(p1, p6, p2, p7, p3, p8, p4, p9, p5, p10, nrows = 5, shareY = FALSE, shareX = TRUE, margin = c(0, 0.01, 0, 0.02))
      })
      
    output$coverage_plot <- renderPlotly({ plotInput()  })
    
    #export image locally
    output$imagesave1 <- downloadHandler(
                            filename = function() { paste0("Epas1.pdf")},
                            content  = function(file) { export(plotInput(), file = file)  })
    
    output$imagesave2 <- downloadHandler(
                            filename = function() { paste0("Epas1.png")},
                            content  = function(file) { export(plotInput(), file = file, vheight = 400, vwidth = 800, zoom = 4)  })


    }
    ui <- fluidPage(   fluidRow(div(column(width=12, plotlyOutput("coverage_plot", height = 400))), style = "padding:25px;"),
                       fluidRow(column(width=4, tags$head(tags$style("#gene_table  {white-space: nowrap;  }")), div(dataTableOutput("gene_table"), style = "font-size: 80%")),
                                column(width=4, selectInput(inputId = "data_set", label = "Select the data set:", choices = c("overlap", "full_gene_set"), selected = "overlap"),
                                                downloadButton("imagesave1", "Save the image as pdf", style = "background-color:#084594; color:white"),
                                                downloadButton("imagesave2", "Save the image as png", style = "background-color:#005a32; color:white")))
                       
          )

shinyApp(ui = ui, server = server)



write.table(x = plot_data[["liver"]], file = "Epas1_liver_plot_data.txt", sep = "\t")
write.table(x = plot_data[["kidney"]], file = "Epas1_kidney_plot_data.txt", sep = "\t")
write.table(x = plot_data[["sequence"]], file = "Epas1_sequence_plot_data.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)



