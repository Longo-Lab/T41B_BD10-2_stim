#!/usr/bin/env Rscript --vanilla
# WAJ 2021-11-29

library(shiny)
library(tidyverse)

# Read in datasets:
fileNames <- list.files(pattern = "annotated.csv$", full.names = F, recursive = T)
fileDirs <- paste0(getwd(), "/", fileNames)
files <- lapply(fileDirs, FUN = function(i) {
  read_csv(i, show_col_types = F, col_select = c("GeneSymbol", "log2FoldChange")) %>%
    filter(abs(log2FoldChange) < 15 & !is.na(GeneSymbol))})

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
  str_remove_all(".*/") %>%
  str_to_title()


# Preparing the App ------------------------------------------------------
joinLM <- function(i, j) {
  # Joins two datasets and fits a linear regression to their L2FC's
  results.total <- inner_join(files[[i]], files[[j]], by = "GeneSymbol",
                              suffix = c("_X", "_Y"))

  return(lm(log2FoldChange_Y ~ log2FoldChange_X, data = results.total))
}

tb <- tibble(
  x = character(),
  y = character(),
  Comparison = character(),
  r2Value = numeric()
)
for (i in 1:(length(files) - 1)) {
  for (j in (i + 1):length(files)) {
    lreg <- joinLM(i, j)
    tb <- add_row(tb,
                  x = paste0(i),
                  y = paste0(j),
                  Comparison = paste0(names(fileDirs)[[i]], " || ", names(fileDirs)[[j]]),
                  r2Value = round(summary(lreg)$adj.r.squared, 3))
  }
}

tb <- tb %>%
  arrange(desc(r2Value)) %>%
  mutate(rValue = round(sqrt(r2Value), 3)) %>%
  select("Comparison", "r2Value", "rValue") %>%
  head(10)

# Actual App Parameters --------------------------------------------------

ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "litera"),
  titlePanel("Log(2) Fold Change Correlation Plots"),
  sidebarLayout(
    sidebarPanel(
      selectInput("choiceX", "Choose an X Axis:",
                  choices = fileDirs),
      selectInput("choiceY", "Choose a Y Axis:",
                  choices = fileDirs)
    ),
    mainPanel(
      plotOutput("plot", width = "600px", height = "600px")
    )
  ),
  "Top 10 Correlations:",
  tableOutput("topCorrelations")
)

server <- function(input, output, session) {
  output$topCorrelations <- renderTable({tb}, striped = T, hover = T, bordered = T)

  output$plot <- renderPlot({
    tibble_x <- read_csv(input$choiceX, show_col_types = F, col_select = c("GeneSymbol", "log2FoldChange")) %>%
      filter(abs(log2FoldChange) < 15 & !is.na(GeneSymbol))
    tibble_y <- read_csv(input$choiceY, show_col_types = F, col_select = c("GeneSymbol", "log2FoldChange")) %>%
      filter(abs(log2FoldChange) < 15 & !is.na(GeneSymbol))
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

  session$onSessionEnded(function() {
    stopApp()
  })
}

options(shiny.launch.broswer = TRUE)
shinyApp(ui, server)
