#!/usr/bin/env Rscript --vanilla
# WAJ 2021-11-29
library(shiny)
library(tidyverse)

# Tell R where it can find the goods:
setwd("data")

# Read in files:
fileNames <- list.files()
fileDirs <- paste0(getwd(), "/", fileNames)

# Make the file names human readable; they become figure axis labels later on
names(fileDirs) <- fileNames %>%
  str_remove("results_") %>%
  str_remove("_annotated.csv") %>%
  str_replace_all("_", " ") %>%
  str_replace_all("vs", " vs. ") %>%
  str_replace_all("  ", " ") %>%
  str_replace_all("TG", "Transgenic") %>%
  str_replace_all("WT", "Wild Type") %>%
  str_replace_all("VEH", "Vehicle") %>%
  str_replace_all("NOSTIM", "noSTIM") %>%
  str_replace_all("noSTIM", "No Stimulation") %>%
  str_to_title()

ui <- fluidPage(
  selectInput("choiceX", "Choose an X-Axis",
              choices = fileDirs, width = "450px"),
  selectInput("choiceY", "Choose a Y-Axis",
              choices = fileDirs, width = "450px"),
  plotOutput("plot", width = 600, height = 600)
)

server <- function(input, output, session) {
  output$plot <- renderPlot({
    tibble_x <- read_csv(input$choiceX) %>%
      filter(abs(log2FoldChange) < 15) %>%
      select("GeneSymbol", "log2FoldChange")
    tibble_y <- read_csv(input$choiceY) %>%
      filter(abs(log2FoldChange) < 15) %>%
      select("GeneSymbol", "log2FoldChange")

    results.total <- inner_join(tibble_x, tibble_y,
                                by = "GeneSymbol", suffix = c("_X", "_Y"))

    lreg <- lm(log2FoldChange_Y ~ log2FoldChange_X, data = results.total)

    ggplot(results.total, aes(x = log2FoldChange_X, y = log2FoldChange_Y)) +
      geom_vline(xintercept = 0, color = "black") +
      geom_hline(yintercept = 0, color = "black") +
      geom_point(size = 0.5) +
      geom_smooth(method = "lm", formula = "y ~ x", se = F, size = 0.5, color = "orange") +
      geom_hline(yintercept = c(1,-1), linetype = "dashed", color = "blue") +
      geom_vline(xintercept = c(1,-1), linetype = "dashed", color = "blue") +
      theme_bw() +
      labs(x = "Log(2) Fold Change - Chosen X Axis",
           y = "Log(2) Fold Change - Chosenn Y Axis",
           title = "L2FC Correlation Plots",
           subtitle = paste0("y = ", round(lreg$coefficients[2], 3), "x + ",
                             round(lreg$coefficients[1], 3), " || R2 = ",
                             round(summary(lreg)$adj.r.squared, 3)))
  })
}


# Set options and get it goin'!
options(shiny.launch.browser = TRUE)
shinyApp(ui, server)
