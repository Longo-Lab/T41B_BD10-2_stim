library(tidyverse)
library(data.table)
library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(gprofiler2)
library(waiter)
library(R.utils)


# ----------------
# Define settings
# ----------------

round_num <- 'R1'

# Project selection
projs <- c(
  'Stim effect by group' = 'T41B_stim'
)


img_dir <- 'imgs'
addResourcePath(img_dir, '.')

# Define functions
add_count <- function(x, s1 = ',', s2 = '|') {
  str_c(str_pad(str_count(x, s1) + 1, width = 2, side = 'left', pad = '0'), x, sep = s2)
}

get_tab_box <- function(typeout, cluster, l2fc_correlations, gprofilers) {
  print(str_c('get_tab_box() for ', cluster, '...'))
  
  # Transcriptomics modules
  imgs <- str_c(
    file.path(img_dir, typeout, round_num), 
    cluster,
    c('TREAT-AD', 'Mostafavi_etal', 'Milind_etal', 'Wan_etal'), 
    'Modules_Up-Down.full.logfdr.png',
    sep = '.'
  )
  names(imgs) <- c('Treat-AD', 'Mostafavi, et al.', 'Milind, et al.', 'Wan, et al.')
  
  transcriptomics_modules <- list(tabPanel(
    'Transcriptomics enrichment', 
    do.call(div, c(lapply(names(imgs), function(i) { a(`data-value` = imgs[[i]], class = ifelse(i == 'Treat-AD', 'active', ''), i) }), class = 'modules')),
    img(src = imgs[['Treat-AD']])
  ))
  
  # Proteomics modules
  imgs <- str_c(
    file.path(img_dir, typeout, round_num), 
    cluster,
    c('TMT-AD', 'BA6', 'BA37'), 
    'Modules_Up-Down.full.logfdr.png',
    sep = '.'
  )
  names(imgs) <- c('Tmt-AD', 'BA6 resilience', 'BA37 resilience')
  
  proteomics_modules <- list(tabPanel(
    'Proteomics enrichment', 
    do.call(div, c(lapply(names(imgs), function(i) { a(`data-value` = imgs[[i]], class = ifelse(i == 'Tmt-AD', 'active', ''), i) }), class = 'modules')),
    img(src = imgs[['Tmt-AD']])
  ))
  
  # TREAT-AD correlation
  img_corr <- str_c(
    file.path(img_dir, typeout, round_num), 
    cluster,
    'TREAT-AD.correlations.png',
    sep = '.'
  )
  
  biodomain_correlation <- list(tabPanel('Treat-AD correlation', img(src = img_corr)))
  
  # L2FC correlation
  l2fc_corr <- list(tabPanel(
    'L2FC correlation',
    radioButtons('l2fc_comparison', 'Comparison', c('WtV vs. (TgV - WtV)' = 'WtV', 'TgD vs. (TgV - TgD)' = 'TgD')),
    uiOutput('l2fc_corr')
  ))
  
  # gProfiler2
  gp_names <- c('LTP-dependent', 'Shared')

  missing_data <- div(style = 'height:300px;', 'No data to show.')
  
  gprofilers <- lapply(seq(gp_names), function(i) {
    g <- gprofilers[[gp_names[[i]]]][['g']]
    g_tbl <- gprofilers[[gp_names[[i]]]][['g_tbl']]
    
    g_title <- str_c(gp_names[[i]], ' gain')
    l_title <- str_c(gp_names[[i]], ' loss')
    
    g_length <- 0
    l_length <- 0
    
    if (length(g) > 1) {
      g_length <- length(g$meta$query_metadata$queries$g)
      l_length <- length(g$meta$query_metadata$queries$l)
    }
    
    sketch <- htmltools::withTags(table(
      class = 'display',
      thead(
        tr(
          th(rowspan = 2, 'Source'),
          th(rowspan = 2, 'Term ID'),
          th(rowspan = 2, 'Term name'),
          th(rowspan = 2, 'Term size'),
          th(colspan = 2, style = 'background-color:#e1e8e3;border-bottom:none;text-align: center;', str_c(g_title, ' (N=', g_length, ')')),
          th(colspan = 2, style = 'background-color:#e8e1e1;border-bottom:none;text-align: center;', str_c(l_title, ' (N=', l_length, ')'))
        ),
        tr(
          lapply(rep(c('Intersection', 'Pval'), 2), th)
        )
      )
    ))
    
    render_sort <- "
    function(data, type, row, meta) {
      if (type === 'sort' && data === null) {
        return meta.settings.aaSorting[0][1] === 'asc' ? 999 : -999;
      }
      return data === null ? data : data.toPrecision(3);
    }
    "
    
    render_tooltip <- "
    function(row, data) {
      // 5th & 7th columns (intersections)
      for (let i = 4; i < 7; i+=2) {
        let genes = data[i] ? data[i].split('|')[1].split(',') : [];
        $('td:eq(' + i + ')', row).
          html('<span class=\"tbl-info\" title=\"' + genes.join('&#010;') + '\">' + genes.length + '</span>');
      }
    }
    "
    
    out_plot <- missing_data
    if (length(g) > 1) out_plot <- gostplot(g, capped = F, interactive = T)
    
    out_tbl <- missing_data
    if (length(g_tbl) > 1) {
      out_tbl <- g_tbl %>% 
        datatable(
          height = '300px',
          rownames = F,
          filter = 'top',
          container = sketch,
          class = 'compact',
          escape = F,
          options = list(
            dom = 'tip',
            scrollX = T,
            scrollY = T,
            rowCallback = htmlwidgets::JS(render_tooltip),
            columnDefs = list(
              list(targets = c(4, 6), className = 'dt-right'),
              list(targets = c(5, 7), render = htmlwidgets::JS(render_sort))
            )
          )
        ) %>%
        formatRound('term_size', 0)
    }
    
    list(
      tabPanel(str_c('gProfiler plot (', ifelse(gp_names[[i]] == 'Shared', 'shared', gp_names[[i]]), ')'), out_plot),
      tabPanel(str_c('gProfiler table (', ifelse(gp_names[[i]] == 'Shared', 'shared', gp_names[[i]]), ')'), out_tbl)
    )
  }) %>% flatten()
  
  # GO terms enrichment
  go_terms <- list(tabPanel(
    'GO terms enrichment',
    div(
      selectInput('go_term', label = 'GO terms category', choices = c('GO:BP', 'GO:MF', 'GO:CC')),
      numericInput('top_n_terms', 'Top N terms', value = 5, min = 5, max = 20),
      sliderInput('term_size_range', label = 'Term size filter', value = c(0, 2000), min = 0, max = 2000),
      textOutput('term_size_text') %>% tagAppendAttributes(style = 'font-weight: 600;')
    ) %>% tagAppendAttributes(style = 'width: 300px;'),
    plotOutput('go_terms_plot', height = 800) %>% tagAppendAttributes(style = 'margin-left: 100px;')
  ))
  
  # tabBox object
  m <- do.call(tabBox, c(transcriptomics_modules, proteomics_modules, biodomain_correlation, l2fc_corr, gprofilers, go_terms))
  m$attribs$class <- 'col-sm-12'
  
  m
}


# -----------------
# Create dashboard
# -----------------

ui <- dashboardPage(
  skin = 'black',
  header = dashboardHeader(
    tags$li('Dashboard by: Robert R Butler III & Crystal Han', br(), 'Please cite:', a('TBD', href = '#'), class = 'dropdown'),
    title = 'T41B + BD10-2 + stim'
  ),
  sidebar = dashboardSidebar(
    div(
      selectInput('proj', label = 'Select:', choices = projs),
    ),
    uiOutput('page_tab'),
    uiOutput('page_legend')
  ),
  body = dashboardBody(
    tags$head(
      tags$link(rel = 'stylesheet', type = 'text/css', href = 'style.css'),
      tags$script(src = 'script.js')
    ),
    useWaiter(),
    uiOutput('page_content')
  )
)

server <- function(input, output, session) {
  
  # Initialize loader
  w <- Waiter$new(
    id = 'page_content',
    html = spin_throbber(), 
    color = transparent(0.5)
  )
  
  # Update w/ cell type selection
  page_data <- reactive({
    data_file <- file.path(input$proj, str_c(round_num, input$proj, 'dashboard_files', 'rdata', sep = '.'))
    
    print(str_c('Loading data at ', data_file, '...'))
    load(data_file)
    
    list(
      'results' = results,
      'l2fc_correlations' = l2fc_correlations, 
      'gprofilers' = gprofilers,
      'go_terms' = go_terms,
      'meta' = meta
    )
  })
  
  # Render sidebar tabs
  output$page_tab <- renderUI({
    w$show()
    
    cls <- names(page_data()[['results']])
    subclass <- cls[[1]]
    
    pages <- lapply(
      cls,
      function(cl) {
        m <- menuItem(
          cl, 
          tabName = str_c('page_', cl), 
          icon = icon(if_else(cl == subclass, 'brain', 'dna'))
        )
        
        if (cl == subclass) {
          m <- m %>% tagAppendAttributes(class = 'active')
        }
        
        m
      }
    )
    
    do.call(sidebarMenu, pages)
  })
  
  # Render sidebar text
  output$page_legend <- renderUI({
    geno <- page_data()[['meta']][['geno']]
    drug <- page_data()[['meta']][['drug']]
    footnote <- page_data()[['meta']][['footnote']]
    
    div(
      p(str_c('WT_VEH_TBS vs. WT_VEH'), br(), icon('arrow-right'),
              span(str_c(' Stimulation effect in WT group'))),
      p(str_c(geno, '_VEH_TBS vs. ', geno, '_VEH'), br(), icon('arrow-right'),
              span(str_c(' Stimulation effect in ', geno, ' group'))),
      p(str_c(geno, '_', drug, '_TBS vs. ', geno, '_', drug), br(), icon('arrow-right'),
              span(str_c(' Stimulation effect in ', geno, '_', drug, ' group'))),
      # Wt drug group, uncomment if present
      # p(str_c('WT_', drug, '_TBS vs. WT_', drug), br(), icon('arrow-right'),
              # span(str_c(' Stimulation effect in WT_', drug, ' group'))),
      br(),
      footnote,
      br(),
      p(paste('Copyright \U00A9 Longo Lab', format(Sys.Date(), '%Y')))
    )
  })
  
  # Render tabs content
  output$page_content <- renderUI({
    cls <- names(page_data()[['results']])
    subclass <- cls[[1]]
    
    de_names <- page_data()[['meta']][['de_names']]
    geno <- page_data()[['meta']][['geno']]
    drug <- page_data()[['meta']][['drug']]
    footnote <- page_data()[['meta']][['footnote']]
    
    analyses_cols <- list(
      htmltools::tags$th(colspan = length(de_names), style = 'background-color:#f3f7eb;border-bottom:none;text-align: center;', title = 'Stimulation effect in WT group', 'WT_VEH_TBS vs WT_VEH'),
      htmltools::tags$th(colspan = length(de_names), style = 'background-color:#fdf2f1;border-bottom:none;text-align: center;', title = str_c('Stimulation effect in ', geno, ' group'), str_c(geno, '_VEH_TBS vs ', geno, '_VEH')),
      htmltools::tags$th(colspan = length(de_names), style = 'background-color:#eef8f9;border-bottom:none;text-align: center;padding: 4px 15px;', title = str_c('Stimulation effect in ', geno, '_', drug, ' group'), str_c(geno, '_', drug, '_TBS vs ', geno, '_', drug))
    )
    
    sketch <- htmltools::withTags(table(
      class = 'display',
      thead(
        tr(
          th(rowspan = 2, 'Ensembl ID'),
          th(rowspan = 2, 'Symbol'),
          th(rowspan = 2, 'Category'),
          analyses_cols,
          th(colspan = 5, style = 'background-color:#eeeeee;border-bottom:none;text-align: center;', 'Modules'),
          th(colspan = 2, style = 'background-color:#f9f3d7;border-bottom:none;text-align: center;', 'TF')
        ),
        tr(
          lapply(rep(de_names, length(analyses_cols)), th),
          th('TreatAD'),
          th('TmtAD'),
          th('Mostafavi'),
          th('Milind'),
          th('Wan'),
          th('AnimalTFDB'),
          th('ChEMBL')
        )
      )
    ))
    
    modules_col <- 3 + length(analyses_cols) * length(de_names)
    
    render_tooltip <- sprintf("
    function(row, data) {
      // 2nd column (Ensembl gene info)
      var info = data[1].split('|'),
          title = 'Chromosome: ' + info[2] + '&#010;Biotype: ' + info[1];
      
      $('td:eq(1)', row)
        .html('<span class=\"tbl-info\" title=\"' + title + '\">' + info[0] + '</span>');
        
      // Module columns
      var m = %s;
      for (var i = m; i < m + 5; i++) {
        var modules = data[i] ? data[i].split('|')[1].split(',') : [];
        $('td:eq(' + i + ')', row)
          .html('<span class=\"tbl-info\" title=\"' + modules.join('&#010;') + '\">' + modules.length + '</span>');
      }
      
      // TF columns
      var c = '';
      if (data[m + 6] !== null) {
        var info = data[m + 6].split('~'),
            num = parseInt(info[0]),
            labels = info[1].split('|'),
            source = info[2];
      
        c = '<span class=\"tbl-info\" title=\"Source: ' + source + '&#010;&#010;' + labels.join('&#010;') + '\">' + num + '</span>';
      }
        
      $('td:eq(' + (m + 6) + ')', row)
        .html(c);
    }
    ", modules_col)
    
    pages <- lapply(
      cls,
      function(cl) {
        results <- page_data()[['results']][[cl]]
        l2fc_correlations <- page_data()[['l2fc_correlations']][[cl]]
        gprofilers <- page_data()[['gprofilers']][[cl]]

        tabItem(
          tabName = str_c('page_', cl),
          class = if_else(cl == subclass, 'active', ''),
          fluidRow(
            get_tab_box(subclass, cl, l2fc_correlations, gprofilers)
          ),
          fluidRow(
            tabBox(
              width = 12,
              tabPanel(
                'Genes list',
                results %>% 
                  mutate(
                    ensembl_gene_id = str_c('<a href="https://www.ncbi.nlm.nih.gov/gene/?term=', ensembl_gene_id, '">', ensembl_gene_id, '</a>'),
                    external_gene_name = str_c(external_gene_name, gene_biotype, chromosome_name, sep = '|'),
                    category = factor(category),
                    treatAD = add_count(treatAD),
                    tmtAD = add_count(tmtAD),
                    mostafavi = add_count(mostafavi),
                    milind = add_count(milind),
                    wan = add_count(wan),
                    chembl = str_c(add_count(chembl_protein_labels, '\\|', '~'), if_else(is_mouse_chembl, 'Mouse', 'Human'), sep = '~')
                  ) %>% 
                  dplyr::select(-gene_biotype, -chromosome_name, -is_mouse_chembl, -chembl_protein_labels) %>% 
                  datatable(
                    selection = 'none',
                    rownames = F, 
                    filter = 'top',
                    container = sketch,
                    class = 'compact',
                    escape = F,
                    callback = htmlwidgets::JS("$(\"[data-type='number'] input[type='search']\").attr('placeholder','');"),
                    options = list(
                      dom = 'tip',
                      scrollX = T,
                      scrollY = T,
                      rowCallback = htmlwidgets::JS(render_tooltip),
                      columnDefs = list(
                        list(targets = modules_col:(modules_col + 6), className = 'dt-right')
                      ),
                      search = list(regex = T)
                    )
                  ) %>% 
                  formatRound(grep('log2|L2FC', names(results), value = T), 3) %>%
                  formatSignif(grep('adj|val', names(results), value = T), 3) %>% 
                  DT::renderDataTable()
              ),
              tabPanel(
                'Summary table',
                results[order(category), .N, by = c('category')] %>%
                  datatable(
                    width = '400px',
                    selection = 'none',
                    rownames = F,
                    colnames = c('Category', 'N'),
                    class = 'compact display',
                    options = list(
                      dom = 't'
                    )
                  ),
                footnote
              )
            )
          )
        )
      }
    )
    
    do.call(tabItems, pages)
  })
  
  # Render L2FC correlation
  output$l2fc_corr <- renderUI({
    l2fc_correlations <- page_data()[['l2fc_correlations']][[1]]
    l2fc_correlations[[input$l2fc_comparison]]
  })
  
  # Render GO terms enrichment
  output$term_size_text <- renderText({
    paste0('Showing top ', input$top_n_terms, ' up and down-regulated ', input$go_term, ' terms in each experimental group. Including only terms whose size is between ', input$term_size_range[[1]], ' and ', input$term_size_range[[2]], '.')
  })
  
  output$go_terms_plot <- renderPlot({
    go_terms <- page_data()[['go_terms']][[1]]
    geno <- page_data()[['meta']][['geno']]
    drug <- page_data()[['meta']][['drug']]

    go_terms <- go_terms[input$term_size_range[[1]] < term_size & term_size < input$term_size_range[[2]]]

    go_terms <- go_terms[source == input$go_term]  # subset terms for specific category

    top_terms <- go_terms[, head(.SD, input$top_n_terms), by = c('query', 'Tiss')][, term_name]

    go_terms <- go_terms[term_name %in% top_terms] %>%
      complete(Tiss, term_name, fill = list(LogBH = 0)) %>%
      data.table()

    term_order <- go_terms[go_terms[, .I[which.max(abs(LogBH))], by=term_name]$V1]
    term_down <- term_order[order(LogBH)][LogBH < 0]$term_name
    term_up <- term_order[order(-LogBH)][LogBH > 0]$term_name
    term_order <- c(term_down, term_up)

    go_terms$term_name = factor(go_terms$term_name, levels=rev(term_order))
    go_terms$Tiss = factor(go_terms$Tiss, levels=rev(c('Wt', geno, str_c(geno, ' + ', drug))))

    ggplot(go_terms, aes(x = LogBH, y = term_name, fill = Tiss)) +
      geom_bar(stat = 'identity', position = 'dodge') +
      guides(fill = guide_legend(reverse = T)) +
      geom_vline(xintercept = 0) +
      scale_fill_brewer(type = 'qual', palette = 'Set2') +
      labs(fill = 'TBS exp. group') +
      xlab(expression('Downregulated' %<-% '-Log10 FDR' %->% 'Upregulated')) +
      geom_text(
        aes(x = 0, label = term_name),
        hjust = 0,
        nudge_x = 0.3,
        size = 3.5,
        color = ifelse(go_terms$term_name %in% term_up, '#f0f0f0', 'black')
      ) +
      theme(
        panel.background = element_rect(fill = 'white'),
        panel.grid.major.x = element_line(color = '#A8BAC4', size = 0.3),
        axis.ticks.length = unit(0, 'mm'),
        axis.title.y = element_blank(),
        axis.line.y.left = element_line(color = '#A8BAC4'),
        axis.text.y = element_blank()
      )
  })
}

shinyApp(ui, server)
