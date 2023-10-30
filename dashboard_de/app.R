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
is_sc <- F

# Project selection
projs <- c(
  'Drug effect (stim)' = 'T41B_BD10-2_stim',
  'Drug effect (unstim)' = 'T41B_BD10-2_unstim',
  'Stim effect' = 'T41B_TBS'
)


img_dir <- 'imgs'
addResourcePath(img_dir, '.')

# Define functions
add_count <- function(x, s1 = ',', s2 = '|') {
  str_c(str_pad(str_count(x, s1) + 1, width = 2, side = 'left', pad = '0'), x, sep = s2)
}

get_tab_box <- function(typeout, cluster, l2fc_correlations, gprofilers, gseas) {
  print(str_c('get_tab_box() for ', cluster, '...'))
  
  # Transcriptomics modules
  imgs <- str_c(
    file.path(img_dir, typeout, round_num), 
    ifelse(is_sc, str_c(nameset, '.', cluster), cluster),
    c('Mostafavi_etal', 'Milind_etal', 'Wan_etal'), 
    'Modules_Up-Down.full.logfdr.png',
    sep = '.'
  )
  names(imgs) <- c('Mostafavi, et al.', 'Milind, et al.', 'Wan, et al.')
  
  transcriptomics_modules <- list(tabPanel(
    'Transcriptomics enrichment', 
    do.call(div, c(lapply(names(imgs), function(i) { a(`data-value` = imgs[[i]], class = ifelse(i == 'Mostafavi, et al.', 'active', ''), i) }), class = 'modules')),
    img(src = imgs[['Mostafavi, et al.']])
  ))
  
  # Proteomics modules
  imgs <- str_c(
    file.path(img_dir, typeout, round_num), 
    ifelse(is_sc, str_c(nameset, '.', cluster), cluster),
    c('TMT-AD', 'BA6'), 
    'Modules_Up-Down.full.logfdr.png',
    sep = '.'
  )
  names(imgs) <- c('Tmt-AD', 'BA6/BA37 resilience')
  
  proteomics_modules <- list(tabPanel(
    'Proteomics enrichment', 
    do.call(div, c(lapply(names(imgs), function(i) { a(`data-value` = imgs[[i]], class = ifelse(i == 'Tmt-AD', 'active', ''), i) }), class = 'modules')),
    img(src = imgs[['Tmt-AD']])
  ))
  
  # TREAT-AD correlation
  img_corr <- str_c(
    file.path(img_dir, typeout, round_num), 
    ifelse(is_sc, str_c(nameset, '.', cluster), cluster),
    'TREAT-AD.correlations.png',
    sep = '.'
  )
  
  biodomain_correlation <- list(tabPanel('Treat-AD correlation', img(src = img_corr)))
  
  # L2FC correlation
  l2fc_corr <- list(tabPanel('L2FC correlation', l2fc_correlations))
  
  # gProfiler2
  gp_names <- c('all', 'direct', 'compensatory')
  gp_abbrevs <- c('indir + dir', 'dir', 'comp')
  
  missing_data <- div(style = 'height:300px;', 'No data to show.')
  
  gprofiler_tabs <- lapply(seq(gp_names), function(i) {
    g <- gprofilers[[gp_names[[i]]]][['g']]
    g_tbl <- gprofilers[[gp_names[[i]]]][['g_tbl']]
    
    gp_n <- ifelse(gp_names[[i]] == 'all', '', gp_names[[i]])

    e_title <- str_to_sentence(str_trim(str_c(gp_n, ' enhancement')))
    s_title <- str_to_sentence(str_trim(str_c(gp_n, ' suppression')))
    
    e_length <- 0
    s_length <- 0
    
    if (length(g) > 1) {
      e_length <- length(g$meta$query_metadata$queries$e)
      s_length <- length(g$meta$query_metadata$queries$s)
    }
    
    gprofiler_sketch <- htmltools::withTags(table(
      class = 'display',
      thead(
        tr(
          th(rowspan = 2, 'Source'),
          th(rowspan = 2, 'Term ID'),
          th(rowspan = 2, 'Term name'),
          th(rowspan = 2, 'Term size'),
          th(colspan = 2, style = 'background-color:#e1e8e3;border-bottom:none;text-align: center;', str_c(e_title, ' (N=', e_length, ')')),
          th(colspan = 2, style = 'background-color:#e8e1e1;border-bottom:none;text-align: center;', str_c(s_title, ' (N=', s_length, ')'))
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
        mutate(
          term_id = if_else(str_detect(term_id, '^GO:'), str_c('<a href="https://www.ebi.ac.uk/QuickGO/term/', term_id, '">', term_id, '</a>'), term_id)
        ) %>% 
        datatable(
          selection = 'none',
          height = '300px',
          rownames = F,
          filter = 'top',
          container = gprofiler_sketch,
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
      tabPanel(str_c('gProfiler plot (', gp_abbrevs[[i]], ')'), out_plot),
      tabPanel(str_c('gProfiler table (', gp_abbrevs[[i]], ')'), out_tbl)
    )
  }) %>% flatten()
  
  gprofiler_tabs <- c(gprofiler_tabs[c(1, 3, 5)], gprofiler_tabs[c(2, 4, 6)])
  
  # GSEA
  gsea_sketch <- htmltools::withTags(table(
    class = 'display',
    thead(
      tr(
        th('Biodomain'),
        th('Source'),
        th('Term ID'),
        th('Term name'),
        th('Term size'),
        th('Pval (adj)'),
        th('NES')
      )
    )
  ))
  
  gsea_tabs <- lapply(c('geno', 'drug', 'geno_drug'), function(category) {
    
    gsea_tbl <- gseas[[category]][['enr_tbl']] %>%
      dplyr::select(Biodomain, ONTOLOGY, ID, Description, setSize, p.adjust, NES) %>%
      mutate(
        ID = str_c('<a href="https://www.ebi.ac.uk/QuickGO/term/', ID, '">', ID, '</a>')
      ) %>% 
      datatable(
        selection = 'none',
        height = '300px',
        rownames = F,
        filter = 'top',
        container = gsea_sketch,
        class = 'compact',
        escape = F,
        options = list(
          dom = 'tip',
          scrollX = T,
          scrollY = T
        )
      ) %>%
      formatSignif(c('p.adjust', 'NES'), 3)
    
    list(
      tabPanel(str_c('GSEA plot (', category, ')'), gseas[[category]][['enr_plt']]),
      tabPanel(str_c('GSEA table (', category, ')'), gsea_tbl)
    )
  }) %>% flatten()
  
  gsea_tabs <- c(gsea_tabs[c(1, 3, 5)], gsea_tabs[c(2, 4, 6)])
  
  # GSEA correlation
  gsea_correlation <- list(tabPanel('GSEA correlation', gseas$corr))
  
  # tabBox object
  m <- do.call(tabBox, c(
    transcriptomics_modules, proteomics_modules, biodomain_correlation, l2fc_corr, 
    list(tabPanel('gProfiler', '')), gprofiler_tabs, 
    list(tabPanel('GSEA', '')), gsea_correlation, gsea_tabs)
  )
  m$attribs$class <- 'col-sm-12'
  
  # remove default click event on grouping tabs
  for (idx in c(5, 12)) {
    m$children[[1]]$children[[1]]$children[[idx]]$children[[1]]$attribs$`data-toggle` = ''
    m$children[[1]]$children[[1]]$children[[idx]]$children[[1]]$attribs$`data-bs-toggle` = ''
  }
  
  # add class to collapsing tabs for easier selection
  for (idx in 5:19) {
    cl <- (m$children[[1]]$children[[1]]$children[[idx]]$children[[1]]$children[[2]] %>% str_split(' '))[[1]][[1]] %>% str_to_lower()
    
    if (idx %in% c(5, 12)) {
      cl <- c(cl, 'tab-group')
    }
    
    m$children[[1]]$children[[1]]$children[[idx]]$attribs$class = cl
  }
  
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
    fn <- ifelse(is_sc, str_c(nameset, '.', input$proj), input$proj)
    data_file <- file.path(input$proj, str_c(round_num, fn, 'dashboard_files', 'rdata', sep = '.'))
    
    print(str_c('Loading data at ', data_file, '...'))
    load(data_file)
    
    list(
      'results' = results,
      'l2fc_correlations' = l2fc_correlations, 
      'gprofilers' = gprofilers,
      'gseas' = gseas,
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
    
    veh <- ifelse(drug == 'TBS', 'VEH_', '')
    
    div(
      p(str_c(geno, '_VEH vs. WT_VEH'), br(), icon('arrow-right'),
              span(str_c(' ', geno, ' effect'))),
      p(str_c(geno, '_', veh, drug, ' vs. WT_VEH'), br(), icon('arrow-right'),
              span(str_c(' ', geno, '_', drug, ' effect'))),
      p(str_c(geno, '_', veh, drug, ' vs. ', geno, '_VEH'), br(), icon('arrow-right'), 
              span(str_c(' ', drug, ' effect'))),
      # Wt drug effect group, uncomment if present
      # p(str_c('WT_', veh, drug, ' vs. WT_VEH'), br(), icon('arrow-right'),
      #         span(str_c(' ', drug, ' effect (in WT)'))),
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
    
    analyses <- page_data()[['meta']][['analyses']]
    de_names <- page_data()[['meta']][['de_names']]
    geno <- page_data()[['meta']][['geno']]
    drug <- page_data()[['meta']][['drug']]

    veh <- ifelse(drug == 'TBS', 'VEH_', '')
    
    include_wt <- ifelse(length(analyses) == 4, T, F)
    
    analyses_cols <- list(
      htmltools::tags$th(colspan = length(de_names), style = 'background-color:#f3f7eb;border-bottom:none;text-align: center;', title = str_c(geno, ' effect'), str_c(geno, '_VEH vs. WT_VEH')),
      htmltools::tags$th(colspan = length(de_names), style = 'background-color:#fdf2f1;border-bottom:none;text-align: center;', title = str_c(geno, ' + ', drug, ' effect'), str_c(geno, '_', veh, drug, ' vs. WT_VEH')),
      htmltools::tags$th(colspan = length(de_names), style = 'background-color:#eef8f9;border-bottom:none;text-align: center;', title = str_c(drug, ' effect (in ', geno, ')'), str_c(geno, '_', veh, drug, ' vs. ', geno, '_VEH'))
    )
    
    if (include_wt) {
      analyses_cols[[4]] <- htmltools::tags$th(colspan = length(de_names), style = 'background-color:#f9f2ff;border-bottom:none;text-align: center;', title = str_c(drug, ' effect (in WT)'), str_c('WT_', veh, drug, ' vs. WT_VEH'))
    }
    
    sketch <- htmltools::withTags(table(
      class = 'display',
      thead(
        tr(
          th(rowspan = 2, 'Ensembl ID'),
          th(rowspan = 2, 'Symbol'),
          th(rowspan = 2, 'Category'),
          th(rowspan = 2, 'Direct'),
          analyses_cols,
          th(colspan = 6, style = 'background-color:#eeeeee;border-bottom:none;text-align: center;', 'Modules'),
          th(colspan = 2, style = 'background-color:#f9f3d7;border-bottom:none;text-align: center;', 'TF')
        ),
        tr(
          lapply(rep(de_names, length(analyses_cols)), th),
          th('TreatAD'),
          th('Mostafavi'),
          th('Milind'),
          th('Wan'),
          th('TmtAD'),
          th('Resilience'),
          th('AnimalTFDB'),
          th('ChEMBL')
        )
      )
    ))
    
    modules_col <- 4 + length(analyses) * length(de_names)
    
    render_tooltip <- sprintf("
    function(row, data) {
      // 2nd column (Ensembl gene info)
      var info = data[1].split('|'),
          title = 'Chromosome: ' + info[2] + '&#010;Biotype: ' + info[1];
      
      $('td:eq(1)', row)
        .html('<span class=\"tbl-info\" title=\"' + title + '\">' + info[0] + '</span>');
        
      // Module columns
      var m = %s;
      for (var i = m; i < m + 6; i++) {
        var modules = data[i] ? data[i].split('|')[1].split(',') : [];
        $('td:eq(' + i + ')', row)
          .html('<span class=\"tbl-info\" title=\"' + modules.join('&#010;') + '\">' + modules.length + '</span>');
      }
      
      // TF columns
      var c = '';
      if (data[m + 7] !== null) {
        var info = data[m + 7].split('~'),
            num = parseInt(info[0]),
            labels = info[1].split('|'),
            source = info[2];
      
        c = '<span class=\"tbl-info\" title=\"Source: ' + source + '&#010;&#010;' + labels.join('&#010;') + '\">' + num + '</span>';
      }
        
      $('td:eq(' + (m + 7) + ')', row)
        .html(c);
    }
    ", modules_col)
    
    pages <- lapply(
      cls,
      function(cl) {
        results <- page_data()[['results']][[cl]]
        l2fc_correlations <- page_data()[['l2fc_correlations']][[cl]]
        gprofilers <- page_data()[['gprofilers']][[cl]]
        gseas <- page_data()[['gseas']][[cl]]

        tabItem(
          tabName = str_c('page_', cl),
          class = if_else(cl == subclass, 'active', ''),
          fluidRow(
            get_tab_box(subclass, cl, l2fc_correlations, gprofilers, gseas)
          ),
          fluidRow(
            tabBox(
              width = 12,
              tabPanel(
                'Genes list',
                class = 'genes-list',
                results %>% 
                  mutate(
                    ensembl_gene_id = str_c('<a href="https://www.ncbi.nlm.nih.gov/gene/?term=', ensembl_gene_id, '">', ensembl_gene_id, '</a>'),
                    external_gene_name = str_c(external_gene_name, gene_biotype, chromosome_name, sep = '|'),
                    category = factor(category),
                    treatAD = add_count(treatAD),
                    mostafavi = add_count(mostafavi),
                    milind = add_count(milind),
                    wan = add_count(wan),
                    tmtAD = add_count(tmtAD),
                    resilience = add_count(resilience),
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
                        list(targets = modules_col:(modules_col + 7), className = 'dt-right')
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
                results[order(-category, -direct), .N, by = c('category', 'direct')] %>%
                  datatable(
                    width = '400px',
                    selection = 'none',
                    rownames = F,
                    colnames = c('Category', 'Direct', 'N'),
                    class = 'compact display',
                    options = list(
                      dom = 't'
                    )
                  )
              )
            )
          )
        )
      }
    )
    
    do.call(tabItems, pages)
  })
  
  # Render plot
  output$gsea_corr <- renderText({input$biodomains})
}

shinyApp(ui, server)
