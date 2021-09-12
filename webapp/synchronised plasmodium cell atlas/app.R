library(Seurat)
library(shiny)
library(DT)
library(pheatmap)
library(rintrojs)

options(stringsAsFactors=FALSE)

# setwd("C:/Users/HP/Desktop/sem 8/project/IG/plasmodium/synchronised plasmodium cell atlas")
seurat_obj <- readRDS("combined.rds")
gene_data <- readRDS("gene_counts_new.rds")
data_monocle_controls = read.csv("monocle_controls_new.csv",row.names = 1)
data_monocle_treatment = read.csv("monocle_treated_new.csv", row.names = 1)

sample_file = read.csv("test_file.txt", header = FALSE)
tt=substr(rownames(seurat_obj@meta.data),17,18)

df_frequency_controls = as.data.frame(table(seurat_obj@meta.data[which(substr(rownames(seurat_obj@meta.data),17,18)=="_1"),6]))
colnames(df_frequency_controls) = c("Cluster Number", "Number of Cells (Controls)")

df_frequency_treated = as.data.frame(table(seurat_obj@meta.data[which(substr(rownames(seurat_obj@meta.data),17,18)=="_2"),6]))
colnames(df_frequency_treated) = c("Cluster Number", "Number of Cells (Treated)")

df_frequency = cbind(df_frequency_controls, df_frequency_treated)
df_frequency = df_frequency[-3]

gene1 = read.csv("cluster_specifc_dge.csv")
gene2 = read.csv("essential_markers.csv")
gene3 = read.csv("Stage_specific_markers.csv")
gene4 = read.csv("Target_AP2 transcription factor.csv")
gene5 = read.csv("Target_Gametogenesis.csv")
gene6 = read.csv("Target_Histone acetyltransferase.csv")
gene7 = read.csv("Target_Histone deacetylase.csv")
gene8 = read.csv("Target_Histone genes.csv")
gene9 = read.csv("Target_PEXEL MOTIF GENES.csv")
gene10 = read.csv("Target_Resistance related genes.csv")
gene11 = read.csv("Target_Ribosomal proteins.csv")
gene12 = read.csv("Target_RIFIN and STEVOR.csv")
gene13 = read.csv("Target_stress responsive.csv")
gene14 = read.csv("Target_Var gene list.csv")

gene_set1 = gene1[,1]
gene_set2 = gene2[,1]
gene_set3 = unique(gene3[,1])
gene_set4 = gene4[,1]
gene_set5 = gene5[,1]
gene_set6 = gene6[,1]
gene_set7 = gene7[,1]
gene_set8 = gene8[,1]
gene_set9 = gene9[,1]
gene_set10 = gene10[,1]
gene_set11 = gene11[,1]
gene_set12 = gene12[,1]
gene_set13 = gene13[,1]
gene_set14 = gene14[,1]

##### READING MCS DATA #################
mseurat_obj <- readRDS("MalariaCellAtlas.rds")
mgene_data = readRDS("malaria_gene_counts.rds")

mdf_frequency = as.data.frame(table(mseurat_obj@meta.data[,6]))
colnames(mdf_frequency) = c("Cluster Number", "Number of Cells")



######### UI Code ########

# Define UI for application
ui <- shinyUI(fluidPage(
  introjsUI(),
  introBox(
  # Header or Title Panel 
  titlePanel("Plasmodium Synchronized Cell Atlas"),
  data.step = 1,
  data.intro = "This is the Plasmodium Synchronized Cell Atlas. 
  To investigate cellular heterogeneity, single cell RNA-sequencing (scRNA-seq) was perfomed for 4949 and 6873 
  synchronized Plasmodium cells in control and under temperature stress condition 
  (phenocopying the cyclic bouts of fever experienced during malarial infection). This web-application allows for easy visualisation of the dataset."
  ),
  # Main Panel
  mainPanel(
    tabsetPanel(type="pills",
                tabPanel("Control vs Treated ",
                         fluidRow(
                           br(),
                           column(4,
                                  introBox(
                                    selectInput("var1_1", "1. Select/Type the gene name", choices =c(rownames(seurat_obj))),
                                    data.step = 2,
                                    data.intro ="Select a gene from the dropdown menu. You can also type the name of the gene directly")
                                  ),
                           column(4,
                                  selectInput("var1_2", "2. Select cluster number", choices = c(0:8))
                                  ),
                           
                         ), 
                         mainPanel( 
                           br(),
                           introBox(
                           DT::dataTableOutput("table1", width = "1000px"),
                           data.step = 3,
                           data.intro = "This table displays the RNA count for the selected gene"
                           ),
                           br(),
                           introBox(
                           plotOutput("plot1_1", click = "plot1_1", width = "1000px", height = "350px"),
                           data.step = 4,
                           data.intro = "This plot shows the expression level of RNA for each gene cluster"
                           ),
                           tags$hr(),
                           introBox(
                           plotOutput("plot1_2", click = "plot1_2", width = "1000px", height = "350px"),
                           data.step = 5,
                           data.intro = "This plot is colored according the the expression level of the gene across different cells in the dataset"
                           ),
                           br(),
                           br(),
                         )
                ),
                tabPanel("QC Panel",
                         fluidRow(
                           br(),
                           column(6,
                                  introBox(
                                    selectInput("var2_1", "1. Select the data type to be projected", choices =c("tsne", "umap")),
                                    data.step = 6,
                                    data.intro = "You can select between t-distributed Stochastic Neighbour Embedding (t-SNE) and Umiform Manifold Approximation and Projection(UMAP)"
                                  ),
                           ),
                           column(6,
                                  introBox(
                                    selectInput("var2_2", "2. Select one", choices =c("nCount_RNA", "nFeature_RNA","percent.mt"), selected = "nCount_RNA"),
                                    data.step = 7,
                                    data.intro="The data quality can be ascertained by choosing the different options from this dropdown"
                                  ),
                           ),
                         ),
                         mainPanel(
                           br(),
                           plotOutput("plot2_1", click="plot2_1", width = "1050px", height = "400px"),
                           tags$hr(),
                           plotOutput("plot2_2", click = "plot2_2", width = "1050px", height = "400px"),
                           br(),
                           br()
                         )
                ),
                tabPanel("Gene Set Information",
                         br(),
                         sidebarPanel(
                           introBox(
                           selectInput("var3_1", "1. Select the gene set", 
                                       choices =c("gene_set1","gene_set2", "gene_set3", "gene_set4", "gene_set5", 
                                                  "gene_set6", "gene_set7", "gene_set8", "gene_set9", "gene_set10", 
                                                  "gene_set11", "gene_set12", "gene_set13", "gene_set14")),
                           data.step = 8,
                           data.intro = "A gene set can be selected from the dropdown list"
                           ),
                           
                         ),
                         mainPanel(
                           br(),
                           introBox(
                            plotOutput("heat3_1"),
                            data.step = 9, 
                            data.intro = "This heatmap depicts the average gene expression for each gene of the gene-set cluster wise"
                           ),
                           br(),
                           introBox(
                           DT::dataTableOutput("table3_1"),
                           data.step = 10,
                           data.intro = "This table shows the values corresponding to the heatmap above"
                           ),
                           br(),
                           br(),
                         )
                ),
                tabPanel("Pseudotime analysis",
                         fluidRow(
                           br(),
                           column(4,
                                  introBox(
                                    selectInput("var4", "1. Select the cluster number",
                                                choices = c(0:8)),
                                    data.step = 11,
                                    data.intro = "Select a gene cluster to see its pseudotime trajectory")
                           )
                         ),
                         mainPanel(
                           fluidRow(
                             column(7,plotOutput("plot4_1", width = "700px", height = "450px")),
                             column(7,plotOutput("plot4_2", width = "700px", height = "450px")),
                             br(),
                             br(),
                           )
                         )
                ),
                tabPanel("Upload Gene Set",
                         br(),
                         # Sidebar layout with input and output definitions ----
                         sidebarLayout(
                           
                           # Sidebar panel for inputs ----
                           sidebarPanel(
                             # Input: Select a file ----
                             introBox(
                             fileInput("file1", "Choose .txt file with one gene name in each line",
                                       multiple = FALSE,
                                       accept = c("text/csv",
                                                  "text/comma-separated-values,text/plain",
                                                  ".csv")),
                             data.step = 12,
                             data.intro = "Upload your own gene list here. Only .txt files are allowed, with line containing only one gene name. There should be no extra spaces after the gene"
                             ),
                             introBox(
                             actionButton("go", "Go"),
                             data.step = 13,
                             data.intro = "Click go to start your analysis"
                             ),
                             # Horizontal line ----
                             tags$hr(),
                             introBox(
                             actionButton("test_file", "Use sample file"),
                             data.step = 14,
                             data.intro = "You can also use a sample file by clicking on this button"
                             ),
                           ),
                           mainPanel(
                             br(),
                             plotOutput("heat5", width = "100%"),
                             br(),
                             DT::dataTableOutput("table5_1"),
                             br(),
                             
                           )
                         ),
                         
                ),
                tabPanel("Malaria Cell Atlas",
                         introBox(
                         tabsetPanel(type="pills",
                                     tabPanel("Genewise View ",
                                              fluidRow(
                                                br(),
                                                column(4,
                                                       selectInput("var6_1_1", "1. Select/Type the gene name", choices =c(rownames(mseurat_obj)))),
                                                
                                              ), 
                                              mainPanel( 
                                                br(),
                                                plotOutput("plot6_1_1", click = "plot6_1_1", width = "1000px", height = "350px"),
                                                tags$hr(),
                                                plotOutput("plot6_1_2", click = "plot6_1_2", width = "700px", height = "350px"),
                                                br(),
                                                br(),
                                              )
                                     ),
                                     tabPanel("Upload Gene Set",
                                              br(),
                                              # Sidebar layout with input and output definitions ----
                                              sidebarLayout(
                                                
                                                # Sidebar panel for inputs ----
                                                sidebarPanel(
                                                # Input: Select a file ----
                                                fileInput("file2", "Choose .txt file with one gene name in each line",
                                                          multiple = FALSE,
                                                          accept = c("text/csv",
                                                                     "text/comma-separated-values,text/plain",
                                                                     ".csv")),
                                                
                                                actionButton("go2", "Go"),
                                        
                                                # Horizontal line ----
                                                tags$hr(),
                                                  actionButton("test_file2", "Use sample file"),
                                                ),
                                                mainPanel(
                                                  br(),
                                                  plotOutput("heat6_2_1", width = "100%"),
                                                  br(),
                                                  DT::dataTableOutput("table6_2_2"),
                                                  br(),
                                                  
                                                )
                                              )
                                     )
                         ),
                         data.step = 15,
                         data.intro = "This tab displayes the data from Malaria Cell Atlas"
                         )
                )
                
    )
    
  )
  
)
)  

server <-shinyServer(
  
  function(input, output,session) {
    
    output$plot1_1 <- renderPlot(
      VlnPlot(seurat_obj, features = input$var1_1, split.by = "condition", pt.size = 0, cols=c("#FFFF99", "#FFCCCC"))
    )
    output$plot1_2 <- renderPlot(
      # DimPlot(object = seurat_obj, reduction = "tsne", split.by = "orig.ident", group.by = "seurat_clusters")
      FeaturePlot(seurat_obj,order = T,
                  reduction = "tsne",
                  features = input$var1_1,
                  split.by = "condition",
                  cols = c("gray","orange","red"), 
      )
    )
    output$table1 <- renderDataTable(
      {
        df <- data.frame("Total RNA of gene in Control Cells"= sum(seurat_obj@assays$RNA@counts[input$var1_1,intersect(which(seurat_obj@meta.data[,6]==input$var1_2),which(tt=="_1"))]),
                         "Average RNA of gene per Control Cell" = sum(seurat_obj@assays$RNA@counts[input$var1_1,intersect(which(seurat_obj@meta.data[,6]==input$var1_2),which(tt=="_1"))])/df_frequency_controls[as.integer(input$var1_2) + 1, 2],
                         "Total RNA of gene in Treated Cells"= sum(seurat_obj@assays$RNA@counts[input$var1_1,intersect(which(seurat_obj@meta.data[,6]==input$var1_2),which(tt=="_2"))]),
                         "Average RNA of per Treated Cell" = sum(seurat_obj@assays$RNA@counts[input$var1_1,intersect(which(seurat_obj@meta.data[,6]==input$var1_2),which(tt=="_2"))])/df_frequency_treated[as.integer(input$var1_2) + 1, 2]
        )
        colnames(df) = c("Total RNA of gene in Control Cells","Average RNA of gene per Control Cell",
                         "Total RNA of gene in Treated Cells", "Average RNA of per Treated Cell")
        return(df)
      }, options = list(paging = FALSE, info = FALSE, searching = FALSE)
    )
    output$plot2_1 <- renderPlot(
      DimPlot(object =seurat_obj, reduction = input$var2_1, split.by = "orig.ident", group.by ="seurat_clusters", label = TRUE, label.size = 6 )
      
    )
    output$plot2_2 <- renderPlot(
      VlnPlot(seurat_obj, features = input$var2_2, split.by = "condition", pt.size = 0, cols=c("#FFFF99", "#FFCCCC"))
    )
    # output$table3_1 = renderTable(df_frequency)
    l1 = c("gene_set1","gene_set2", "gene_set3", "gene_set4", "gene_set5", 
           "gene_set6", "gene_set7", "gene_set8", "gene_set9", "gene_set10", 
           "gene_set11", "gene_set12", "gene_set13", "gene_set14")
    l2 = list(gene_set1,gene_set2, gene_set3, gene_set4, gene_set5, 
              gene_set6, gene_set7, gene_set8, gene_set9, gene_set10, 
              gene_set11, gene_set12, gene_set13, gene_set14)
    output$table3_1 = renderDataTable(
      {
        dat = t(na.omit(gene_data[l2[[which(l1==input$var3_1)]],]))
        dat <- ifelse(dat==0, 0.01, dat)
        j=1
        for (i in c(1:nrow(df_frequency)))
        {
          if(j!=17)
          {
            dat[j,] = dat[j,]/as.integer(df_frequency[i,2])
            j=j+1
          }
          dat[j,] = dat[j,]/as.integer(df_frequency[i,3])
          j=j+1
        }
        return(t(dat))
      }
    )
    output$heat3_1 = renderPlot(
      {
        df <- t(na.omit(gene_data[l2[[which(l1==input$var3_1)]],]))
        df = ifelse(df==0,0.01,df)
        j=1
        for (i in c(1:nrow(df_frequency)))
        {
          if(j!=17)
          {
            df[j,] = log10(df[j,]/as.integer(df_frequency[i,2]))
            j=j+1
          }
          df[j,] = log10(df[j,]/as.integer(df_frequency[i,3]))
          j=j+1
        }
        pheatmap(df)
      }, height = 400, width = 1000
    )
    output$plot4_1 = renderPlot(
      plot(data_monocle_controls[,1], data_monocle_controls[,2], 
           col=ifelse(data_monocle_controls[,3]==input$var4, "#FF0000", "#C0C0C0"),
           main="Trajectory: Controls", xlab="UMAP_1", ylab="UMAP_2")
      )
    output$plot4_2 = renderPlot(
      plot(data_monocle_treatment[,1], data_monocle_treatment[,2], 
           col=ifelse(data_monocle_treatment[,3]==input$var4, "#FF0000", "#C0C0C0"),
           main="Trajectory: Treated", xlab="UMAP_1", ylab="UMAP_2")
      )
      
    observeEvent(input$go, {
      output$heat5 = renderPlot(
      {
        req(input$file1)
        tryCatch(
          {
            df <- read.csv(input$file1$datapath, header = FALSE)
            dat <- t(na.omit(gene_data[df[,1],]))
            dat <- ifelse(dat==0, 0.01, dat)
            j=1
            for (i in c(1:nrow(df_frequency)))
            {
              if(j!=17)
              {
                dat[j,] = log10(dat[j,]/as.integer(df_frequency[i,2]))
                j=j+1
              }
              dat[j,] = log10(dat[j,]/as.integer(df_frequency[i,3]))
              j=j+1
            }
          },
          error = function(e) {
            # return a safeError if a parsing error occurs
            stop(safeError(e))
          }
        )
        return(pheatmap(dat))
      }
    )
      output$table5_1 = renderDataTable(
        {
          req(input$file1)
          tryCatch(
            {
              df <- read.csv(input$file1$datapath, header = FALSE)
              dat <- t(na.omit(gene_data[df[,1],]))
              dat <- ifelse(dat==0, 0.01, dat)
              j=1
              for (i in c(1:nrow(df_frequency)))
              {
                if(j!=17)
                {
                  dat[j,] = dat[j,]/as.integer(df_frequency[i,2])
                  j=j+1
                }
                dat[j,] = dat[j,]/as.integer(df_frequency[i,3])
                j=j+1
              }
            },
            error = function(e) {
              # return a safeError if a parsing error occurs
              stop(safeError(e))
            }
          )
          return(t(dat))
        }
      )
    }
    )
    observeEvent(input$test_file, {
      output$table5_1<-renderDataTable(
        {
          dat <- t(na.omit(gene_data[sample_file[,1],]))
          dat <- ifelse(dat==0, 0.01, dat)
          j=1
          for (i in c(1:nrow(df_frequency)))
          {
            if(j!=17)
            {
              dat[j,] = dat[j,]/as.integer(df_frequency[i,2])
              j=j+1
            }
            dat[j,] = dat[j,]/as.integer(df_frequency[i,3])
            j=j+1
          }
          return(t(dat))
        }
        
      )
      output$heat5 <- renderPlot(
        {
          dat <- t(na.omit(gene_data[sample_file[,1],]))
          dat <- ifelse(dat==0, 0.01, dat)
          j=1
          for (i in c(1:nrow(df_frequency)))
          {
            if(j!=17)
            {
              dat[j,] = log10(dat[j,]/as.integer(df_frequency[i,2]))
              j=j+1
            }
            dat[j,] = log10(dat[j,]/as.integer(df_frequency[i,3]))
            j=j+1
          }
          pheatmap(dat)
        }
      )
    })
    
    output$plot6_1_1 <- renderPlot(
      VlnPlot(mseurat_obj, features = input$var6_1_1, pt.size = 0)
    )
    output$plot6_1_2 <- renderPlot(
      DimPlot(object = mseurat_obj, reduction = "umap", label = TRUE)
    )
    
    observeEvent(input$go, {
      output$heat6_2_1 = renderPlot(
        {
          req(input$file2)
          tryCatch(
            {
              df <- read.csv(input$file2$datapath, header = FALSE)
              dat <- t(na.omit(mgene_data[df[,1],]))
              dat <- ifelse(dat==0, 0.01, dat)
              j=
                for (i in c(1:nrow(mdf_frequency)))
                {
                  dat[j,] = log10(dat[j,]/as.integer(mdf_frequency[i,2]))
                  j=j+1
                }
            },
            error = function(e) {
              # return a safeError if a parsing error occurs
              stop(safeError(e))
            }
          )
          return(pheatmap(dat))
        }
      )
      output$table6_2_2 = renderDataTable(
        {
          req(input$file2)
          tryCatch(
            {
              df <- read.csv(input$file2$datapath, header = FALSE)
              dat <- t(na.omit(mgene_data[df[,1],]))
              dat <- ifelse(dat==0, 0.01, dat)
              j=1
              for (i in c(1:nrow(mdf_frequency)))
              {
                dat[j,] = dat[j,]/as.integer(mdf_frequency[i,2])
                j=j+1
              }
            },
            error = function(e) {
              # return a safeError if a parsing error occurs
              stop(safeError(e))
            }
          )
          return(t(dat))
        }
      )
    }
    )
    observeEvent(input$test_file2, {
      output$table6_2_2<-renderDataTable(
        {
          dat <- t(na.omit(mgene_data[sample_file[,1],]))
          dat <- ifelse(dat==0, 0.01, dat)
          j=1
          for (i in c(1:nrow(mdf_frequency)))
          {
            
            dat[j,] = dat[j,]/as.integer(mdf_frequency[i,2])
            j=j+1
            
          }
          return(t(dat))
        }
        
      )
      output$heat6_2_1 <- renderPlot(
        {
          dat <- t(na.omit(mgene_data[sample_file[,1],]))
          dat <- ifelse(dat==0, 0.01, dat)
          j=1
          for (i in c(1:nrow(mdf_frequency)))
          {
            
            dat[j,] = log10(dat[j,]/as.integer(mdf_frequency[i,2]))
            j=j+1
            
          }
          pheatmap(dat)
        }
      )
    })
    introjs(session, options = list("nextLabel"="Next",
                                    "prevLabel"="Back",
                                    "skipLabel"="Skip",
                                    exitOnOverlayClick=F,
                                    exitOnEsc=F),
            events = list(onbeforechange = readCallback("switchTabs"))
            )
    
  }
)

shinyApp(ui, server)

