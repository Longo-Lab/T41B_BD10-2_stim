####################################################################################################
###
###     Shiny app for P75 project - Amira Latif Hernandez
###           Shiny developed by : Patricia Moran Losada pmlosada@stanford.edu
###                                Please if you use this code add my name as a co-authour of your 
###                                research study. 
####################################################################################################


rm(list=ls())
library(shiny)
library(ggplot2)
library(shinythemes)
library(webshot)
library(shinydashboard)
library(DT)
library(plyr)
library(shinyWidgets)
library(forcats)
library(magrittr)
library(dplyr)
library(ggrepel)
library(reshape)
library(Hmisc)
library(corrplot)
library(RColorBrewer)
library(gplots) 
library(markdown)
library(plotly)
library(scales)
library(reshape2)
library(grid)



load("./p75_RSEM_Quant.genes.tpm.RData")
metadata <- read.csv(file="./metadata2.csv",header = TRUE, stringsAsFactors = FALSE,row.names = 1)
gene_list <- read.csv(file="gene_list.csv",header = TRUE, stringsAsFactors = FALSE,row.names = 1)




# Define UI ----
ui <- fluidPage(
  HTML('<meta name="viewport" content="width=1024">'),
  
  theme=shinytheme("sandstone"),
  navbarPage("P75", inverse=FALSE, position = "fixed-top" , 
                  tabPanel("Bulk RNA-seq Visualization",
                      sidebarLayout(
                        sidebarPanel( 
                          tags$style(type='text/css',
                                     ".selectize-dropdown-content{height: 700px;width: 200px;background-color: #b0c4de;
                          }"),
                          width=2,
                          br(),br(),br(),
                          selectInput("Gene", list(h4("Gene:")), choices=c("Trem2"),selected="Trem2",multiple =FALSE ),
                          p("*Type your gene of interest (only first 1000 genes displayed)",align = "left",style = "font-size:10px"),
                          prettyCheckboxGroup("Genotype", label = h4("Genotype"),choices = c("WT","TG"), 
                                              selected = c("WT","TG"),inline   = FALSE, shape = c("curve"), 
                                              thick=FALSE, animation=c("pulse"),fill=FALSE),
                          prettyCheckboxGroup("Treatment", label = h4("Treatment"),choices = c("VEH","C31"), 
                                              selected = c("VEH","C31"),inline   = FALSE, shape = c("curve"), 
                                              thick=FALSE, animation=c("pulse"),fill=FALSE),
                          prettyCheckboxGroup("TBS", label = h4("TBS"),choices = c("TBS", "noTBS"), 
                                              selected = c("TBS","noTBS"),inline   = FALSE, shape = c("curve"), 
                                              thick=FALSE, animation=c("pulse"),fill=FALSE),
                          radioButtons("Graph", list(h4(icon("chart-bar"),"Graph type")),c("Boxplot" = "boxplot","Violin" = "violin")),
                          radioButtons("Scale", list(h4(icon("chart-bar"), "Scale")),c("Linear"="linear","Log"="log"),selected=c("linear"))
                        ),
                        
                        # Main panel for displaying outputs ----
                        mainPanel(
                          br(),br(),br(),br(),
                          
                            fluidRow(br(),
                                   column(5),
                                   column(5,h5(list(h2("  Gene expression"))))),br(),
                                         fluidRow(
                                           column(2),
                                           column(10,
                                                  plotOutput(
                                                    outputId = "plot_bx.var",width ="90%", height = "auto"))),
                                         br(),
                                         fluidRow(column(2,downloadButton(outputId="downloadPlot_bx", label="Download"))),
                                         br(),br(),
                                         DT::dataTableOutput("data"),
                                         br(),br(),br(),br(),br(),br()
                                      
                                      
                          )
                        
                      ))
  ),
  
  navbarPage(inverse=FALSE, position = "fixed-bottom",
             p(a("Please cite:",href=""),align = "justify",style = "font-size:11px"))
)
             

# Define server logic ----
server <- function(input, output,session) {
  
  updateSelectizeInput(session = session, inputId = 'Gene', choices = sort(gene_list$GeneSymbol),selected ="Trem2",  server = TRUE)
  

           data<- reactive({
             
             gene <- gene_list[which(gene_list$GeneSymbol %in% input$Gene),]
             dataSelected <- as.data.frame(tpm[which(rownames(tpm) %in% rownames(gene)),])
             colnames(dataSelected)[1] <- "TPM"
             dataSelected$SampleID <- rownames(dataSelected)
             dataSelected <- dataSelected[rownames(metadata),]
             dataSelected <- cbind(dataSelected,metadata)
             
             metadataSelected <- metadata[which(metadata$Genotype %in% input$Genotype),]
             metadataSelected <- metadataSelected[which(metadataSelected$Treatment %in% input$Treatment),]
             metadataSelected <- metadataSelected[which(metadataSelected$TBS %in% input$TBS),]
             
             validate(
               need(input$Genotype != "", "Please select the Genotype ")
             )
             validate(
               need(input$Treatment != "", "Please select a Treatment")
             )
             validate(
               need(input$TBS != "", "Please select TBS option")
             )
             
             samples <- as.vector(rownames(metadataSelected))
             
             dataSelected <- dataSelected[which(rownames(dataSelected) %in% samples),]
             
             return(dataSelected)
             })
           
           
           dataset_boxplot<- reactive({ 
             mt<-function(x){format(x,nsmall = 0,scientific = FALSE)}
             df <- data()
             
             if(input$Graph== "boxplot"){
               bx<-ggplot(df, aes(x=Group,y=(TPM+0.1))) +geom_boxplot(aes(color=Group),  outlier.shape = NA, lwd=0.8)+
                 geom_jitter(aes(color=Group), position=position_jitter(0.1),size=5,alpha = 0.3)+
                 facet_wrap(.~Genotype, scales = "free_x")+
                 theme(axis.text.x = element_text(angle = 45,hjust = 1))}
             else if(input$Graph == "violin"){
               bx<-ggplot(df, aes(x=Group,y=(TPM+0.1), fill= Group)) +
                 geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")+
                 facet_wrap(.~Genotype, scales = "free_x")+
                 theme(axis.text.x = element_text(angle = 45,hjust = 1))}
             
             bx<- bx + theme(strip.background = element_rect(color="black", fill="gray80", size=1, linetype="solid"), panel.background = element_blank(),
                     panel.border = element_rect(fill = NA, colour = "black")) +
               theme(text = element_text(size=16),axis.title.y = element_text(vjust = 2))  + ylab("TPM+0.1\n")+
               theme(legend.position="none") +xlab("")
             
             if(input$Scale =="linear"){
               bx + scale_y_continuous()+ ylab("TPM+0.1\n")+xlab("")}
             else {
               bx+  scale_y_log10(labels = mt)}
           })
           
           
           
            
           output$plot_bx.var <- renderPlot({
             dataset_boxplot()
           },
           height = function() {
             
             if(length(input$Genotype)>1){
               session$clientData$output_plot_bx.var_width*0.6}else{
                 session$clientData$output_plot_bx.var_width*0.8
               }
             
           })
             
  
           
           output$downloadPlot_bx <- downloadHandler(
             filename = function() {paste("plot.pdf")},
             content = function(file) {
               
                 pdf(file,  width=10, height=8)
                 plot(dataset_boxplot())
                 dev.off()
               }
              
               ) 
           
           output$data <- DT::renderDataTable(DT::datatable({
             data()
           }))
           
}

# Run the app ----
shinyApp(ui = ui, server = server)








