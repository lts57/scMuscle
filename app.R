# scMuscle Explorer
#   written by Leo Song
# Load libraries ----
# library(shiny)
library(Seurat)
library(ggplot2)
library(shinythemes)
library(pals)
library(dplyr)
library(patchwork)

library(shinycssloaders)
library(shades) 

# Load data ----
load("scMuscle_slim_v3.RData")
load("myo_slim_v2.RData") #TODO: don't load this until someone clicks on the myogenic tab
# User interface ----
ui <- fluidPage(
  theme = shinytheme("paper"),
  navbarPage(
    "scMuscle Explorer",
    # Each tabPanel makes a tab on the upper navbar
    # About tab
    tabPanel(
      "About",
      # Title of the app----
      fluidRow(
        style = 'border-bottom: 1px solid grey; border-top: 1px solid grey',
        br(),
        h1("Single-Cell Muscle Explorer", align = "center"),
        br()
      ),
      br(),
      br(),
      # What you'll find here----
      fluidRow(
        style = 'border-bottom: 1px solid grey',
        column(
          6,
          br(),
          h3(
            "This web app is an interactive tool to explore
102 single-cell RNA-sequencing mouse skeletal muscle datasets."
            )
        ),
        column(
          6,
          img(src = "img1.png", width = 700, height = 600),
          br()
        )
      ),
      # How to use the app ----
      fluidRow(
        style = 'border-bottom: 1px solid grey',
        column(
          5,
          br(),
          img(src = "img2.png", width = 500, height = 500),
          br()
        ),
        column(
          7,
          h2("How to use the app", align = "center"),
          br(),
          h3("User-selected inputs appear on the left-hand sidebar panel
to control the figures plotted on the main panel of the app.
The main panel includes four different tabs to showcase
UMAPs, Single Violin plots, Split Violin plots, and Dot Plots.")
        )
      ),
      # Who we are----
      fluidRow(
        column(
          6,
          br(), br(), br(),
          div(img(src = "img3.png", width = 500, height = 100), align = "center")
        ),
        column(
          6,
          h2("Who we are", align = "center"),
          br(),
          h3(p("This app was developed by the",
               a(href = "http://devlaminck.bme.cornell.edu/",
                 "De Vlaminck Lab"), "and the",
               a(href = "https://cosgrovelab.bme.cornell.edu/",
                 "Cosgrove Lab"),
               "of Cornell University Biomedical
                Engineering. Here is the",
               a(href = "https://github.com/lts57/scMuscle",
                 "Github"), "and",
               a(href = "https://www.biorxiv.org/content/10.1101/2020.12.01.407460v1",
                 "preprint", "."))),
          br(), br()
        )
      )
    ),
    # Skeletal Muscle tab----
    tabPanel(
      "All Cell Types",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          # UMAP Panel Inputs----
          conditionalPanel(
            condition = "input.tabselected==1",
            # inputs for reduction to be shown
            br(),
            helpText("Generate UMAPs with the integration of 102 datasets."),

            selectInput(
              "reduction1",
              label = "Choose a dimensional reduction to display",
              choices = c("Harmony", "BBKNN", "Scanorama"),
              selected = "Harmony"
            ),

            # inputs gene for feature plot
            br(),
            selectInput(
              "gene1",
              label = "Select gene to examine:",
              choices = rownames(scMuscle.seurat),
              selected = "Sox17"
            )
          ),
          # Single Violin Panel Inputs----
          conditionalPanel(
            condition = "input.tabselected==2",
            # selects gene to be examined on violin plot
            br(),
            selectInput(
              "gene2",
              label = "Select gene to display",
              choices = rownames(scMuscle.seurat),
              selected = "Sox17"
            )
          ),
          # Split Violin Panel Inputs----
          conditionalPanel(
            condition = "input.tabselected==3",
            # selects gene to be examined on violin plots
            br(),
            selectInput(
              "gene3",
              label = "Select gene to display:",
              choices = rownames(scMuscle.seurat),
              selected = "Sox17"
            ),
            # selects cell types to label on x axis
            br(),
            selectInput(
              "splitviolincelltype",
              label = "Select cell type IDs:",
              choices = c("Harmony", "BBKNN", "Scanorama"),
              selected = "Harmony"
            )
          ),
          # Dot Panel Inputs----
          conditionalPanel(
            condition = "input.tabselected==4",
            # dot plot action button
            br(),
            helpText("Click to apply input changes below"),
            actionButton(
              "action1", label = "Generate Dot Plot",
              style="color: #B31B1B; background-color: #F7F7F7; border-color: #B31B1B"
            ),
            # inputs for dot plot
            br(), br(), br(),
            selectizeInput(
              "dot",
              label = "Select gene(s) to display:",
              choices = rownames(scMuscle.seurat),
              selected = "Sox17",
              multiple = TRUE
            )
          ),
          # Inputs present in all tabs----
          # metadata features
          br(),
          helpText("Group cells by metadata features..."),

          selectInput(
            "variables",
            label = "Metadata Features:",
            choices = c(
              "source", "Harmony cell types", "BBKNN cell types", "Scanorama cell types",
              "sample", "chemistry", "injury days", "injury agent",
              "age", "type", "tissue", "mouse strain", "sex", "mice per sample", "sequencing instrument"
            ),
            selected = "source"
          ),
          # downloadable plot type
          br(),
          helpText("Download Specifications"),
          selectInput(
            "downloadable",
            label = "file type:",
            choices = c("pdf", "png", "eps"),
            selected = "pdf"
          ),
          numericInput(
            "plotsizex",
            label = "horizontal dimension (px)",
            value = 1000
          ),
          numericInput(
            "plotsizey",
            label = "vertical dimension (px)",
            value = 1000
          )
        ),

        # Establishes spaces for plots in the main panel ----
        mainPanel(
          tabsetPanel(
            # umap panel----
            tabPanel(
              "UMAP", value = 1,
              # umap grouped by cell types----
              br(),
              downloadButton("down1", label = "Download"),
              br(), br(),
              plotOutput("umap") %>% withSpinner(type = 1, color = "#B31B1B"),
              br(),
              # umap grouped by metadata features----
              downloadButton("down2", label = "Download"),
              br(), br(),
              plotOutput("feature") %>% withSpinner(type = 1, color = "#B31B1B"),
              br(),
              # feature plot----
              downloadButton("down3", label = "Download"),
              br(), br(),
              plotOutput("grouping") %>% withSpinner(type = 1, color = "#B31B1B"), 
              br() 
            ),

            # single violinplot panel----
            tabPanel(
              "Single Violin", value = 2,
              # by different metadata variables----
              br(),
              downloadButton("down4", label = "Download"),
              br(), br(),
              plotOutput("violin1") %>% withSpinner(type = 1, color = "#B31B1B"),
              br()
            ),
            # split violinplot panel----
            tabPanel(
              "Split Violin", value = 3,
              # makes different violin plots for each unique instance of a metadata variable----
              # grouped by cell types IDs of different reductions
              br(),
              downloadButton("down5", label = "Download"),
              br(), br(),
              imageOutput("violin2") %>% withSpinner(type = 1, color = "#B31B1B"),
              br()
            ),
            # Dot Plot panel----
            tabPanel(
              "Dot Plot", value = 4,
              # DotPlot----
              br(),
              downloadButton("down6", label = "Download"),
              br(), br(),
              plotOutput("dotplot") %>% withSpinner(type = 1, color = "#B31B1B")
            ),
            id = "tabselected"
          )
        )
      )
    ),
    # Myogenic Muscle tab----
    tabPanel(
      "Myogenic Cells",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          # PHATE Panel Inputs----
          # inputs for reduction to be shown
          br(),
          helpText("Visualize gene expression in myogenic cells alone..."),
          selectInput(
            "reduction3",
            label = "Choose a dimensional reduction to display:",
            choices = c("PHATE + Harmony", "PHATE + Scanorama"),
            selected = "PHATE + Harmony"
          ),
          br(),
          helpText("Group cells by metadata features..."),

          selectInput(
            "bins",
            label = "Metadata Features:",
            choices = c("PHATE bins"),
            selected = "PHATE bins"),
          # violin plot action button
          br(),
          helpText("Click to generate new violin plots"),
          actionButton(
            "action2", label = "Generate",
            style="color: #B31B1B; background-color: #F7F7F7; border-color: #B31B1B"
          ),
          # inputs for violin plots
          br(), br(), br(),
          selectizeInput(
            "gene4",
            label = "Select gene(s) to display:",
            choices = rownames(myo.slim.seurat),
            selected = "Sox17",
            multiple = TRUE
          ),
          # downloadable plot type
          br(),
          helpText("Download Specifications"),
          selectInput(
            "downloadable",
            label = "file type:",
            choices = c("pdf", "png", "eps"),
            selected = "pdf"
          ),
          numericInput(
            "plotsizex",
            label = "horizontal dimension (px)",
            value = 1000
          ),
          numericInput(
            "plotsizey",
            label = "vertical dimension (px)",
            value = 1000
          )
        ),
        # Establishes spaces for plots in the main panel ----
        mainPanel(
          # PHATE grouped by variables----
          br(),
          downloadButton("down7", label = "Download"),
          br(), br(),
          plotOutput("PHATE") %>% withSpinner(type = 1, color = "#B31B1B"),
          br(),
          downloadButton("down8", label = "Download"),
          br(), br(),
          imageOutput("phateviolin") %>% withSpinner(type = 1, color = "#B31B1B"),
          br()
        )
      )
    )
  )
)

# Server logic ----
server <- function(input, output) {
  # Plot themes and colors----
  # Figure settings
  small.font = reactive(10)
  big.font = reactive(12)
  line.width = reactive(0.8)
  pt.size = reactive(0.01)
  pt.stroke = reactive(0.3)
  label.size = reactive(4)
  # Color Palette to cell delegation
  colors1 <- reactive({
    as.vector(polychrome())[c(3:9,17,11:13,21,14:16,32,18:20,28,10,30,29,31,33,34,35,36,22:27)]
    })

  # Plots
  umap.theme <- reactive({theme(
    axis.line = element_line(color = "black", size = line.width()),
    axis.title = element_text(face='bold',size = small.font(), hjust = 0.5, vjust = 1),
    axis.text = element_text(size=small.font(),color="black"),
    axis.ticks = element_line(color = "black", size = line.width()),
    legend.text = element_text(size=small.font(),color="black"),
    legend.title = element_text(size=big.font(),color="black")
  )})
  dot.theme <- reactive({theme(
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    legend.text = element_text(size=small.font(), color="black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x=element_text(angle = 90,hjust = 1,vjust= 0.5),
    panel.grid.major = element_line(colour = "gray", size = 0.5)
  )})
  vln.theme <- reactive({theme(
    panel.background = element_blank(),
    axis.line = element_line(color = "black", size = line.width()),
    legend.text = element_text(size=small.font(), color="black"),
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y=element_text(size=big.font(), face="bold",color="black"),
    axis.text = element_text(color="black",size=small.font()),
    axis.ticks.x = element_line(color="black")
  )})

  # Passes UMAP Panel Inputs----

  # passes the reduction type to be plotted
  reduction1 <- reactive({switch(input$reduction1,
                                 "Harmony" = "umap_harmony",
                                 "BBKNN" = "umap_bbknn",
                                 "Scanorama" = "umap_scanorama")})
  umapxlabel <- reactive(paste(reduction1(), "1", sep = "_"))
  umapylabel <- reactive(paste(reduction1(), "2", sep = "_"))

  # creates variable to group by for UMAPs by having the user just decide the reduction
  reduction2 <- reactive(if (reduction1()== "umap_harmony"){"harmony_factorIDs"}
                         else if(reduction1()=="umap_bbknn"){"bbknn_factorIDs"}
                         else if(reduction1()=="umap_scanorama"){"scanorama_factorIDs"}
  )


  # passes gene to feature
  gene1 <- reactive({input$gene1})

  # Passes Single Violin Panel Inputs----

  # passes gene for Violin Plot
  gene2 <- reactive({input$gene2})

  # Passes Split Violin Panel Inputs----

  # passes gene for Multiple Violin Plots
  gene3 <- reactive({input$gene3})

  # allows violin plots to be grouped by chosen cell types
  splitviolincelltype <- reactive(if (input$splitviolincelltype == "Harmony"){"harmony_factorIDs"}
                                  else if(input$splitviolincelltype =="BBKNN"){"bbknn_factorIDs"}
                                  else if(input$splitviolincelltype == "Scanorama"){"scanorama_factorIDs"}
  )

  # Passes Dot Plot Panel Inputs----

  # reverses vector to have dotplot plot genes list left to right
  dot <- reactive(rev(input$dot))

  # Passes Inputs common to all panels----
  # passes different variables to group by
  variables <- reactive({
    switch(
      input$variables,
      "source" = "source.label",
      "Harmony cell types" = "harmony_factorIDs",
      "BBKNN cell types" = "bbknn_factorIDs",
      "Scanorama cell types" = "scanorama_factorIDs",
      "sample" = "sample",
      "chemistry" = "chemistry",
      "injury days" = "injury",
      "injury agent" = "injury.agent",
      "age" = "age",
      "type" = "type",
      "tissue" = "tissue",
      "mouse strain" = "mouse.strain",
      "sex" = "sex",
      "mice per sample" = "mice.per.sample",
      "sequencing instrument" = "Sequencing.Instrument"
    )
  })
  # Passes PHATE Tab Inputs----
  # passes the reduction type to be plotted
  reduction3 <- reactive({
    switch(
      input$reduction3,
      "PHATE + Harmony" = "phate_harmony",
      "PHATE + Scanorama" = "phate_scanorama"
    )
  })
  # passes variable to group by
  bins <- reactive({switch(input$bins,
                           "PHATE bins" = "phate1.bins")})
  # passes gene for Multiple Violin Plots
  gene4 <- reactive({input$gene4})

  # generates first plot (All Cells - Cell Type UMAP) ----
  output$umap <- renderPlot({
    DimPlot(
      scMuscle.seurat,
      cells = sample(colnames(scMuscle.seurat)), #plot cells in random order
      reduction=reduction1(),
      group.by=reduction2(),
      cols=colors1(), #adds colors for just the cell types present in this clustering
      na.value = NA, # removes noisy cells from plot
      pt.size = pt.size(), # see value above
      label.size = label.size(), # see value above
      repel = T,label= TRUE
    ) +
      NoLegend() +
      aes(stroke=pt.stroke())+
      umap.theme()
  })

  # generates second plot (All Cells - metadata UMAP) ----
  output$grouping <- renderPlot({
    DimPlot(
      scMuscle.seurat,
      cells = sample(colnames(scMuscle.seurat)), #plot cells in random order
      reduction=reduction1(),
      group.by=variables(),
      cols=colors1(), #adds colors for just the cell types present in this clustering
      na.value = NA, # removes noisy cells from plot
      pt.size = pt.size(), # see value above
      label.size = label.size(), # see value above
      repel = T,label= TRUE
    ) +
      aes(stroke=pt.stroke())+
      umap.theme()
  })

  # generates third plot (All Cells - FeaturePlot/UMAP) ----
  output$feature <- renderPlot({
    FeaturePlot(
      scMuscle.seurat,
      cells = sample(colnames(scMuscle.seurat)), #plot cells in random order
      features = gene1(),
      reduction = reduction1()
    )+
      scale_colour_viridis_c()+
      labs(color="Log-Normalized\nExpression")+
      umap.theme()
  })


  # generates fourth plot (All Cells - Violin Plot)----
  output$violin1 <- renderPlot({
    VlnPlot( #TODO - add multiple gene plotting
      scMuscle.seurat,
      features = gene2(),
      group.by = variables(),
      cols = colors1(),
      combine = F,
      pt.size = 0
    ) %>% lapply(
      FUN=function(X) X +
        NoLegend() +
        scale_y_continuous(expand=c(0,0)) +
        scale_colour_viridis_c() +
        vln.theme()
    ) %>% wrap_plots(ncol=1)
  })

  # renders image of fifth plot (All Cells - Split Violin Plot)----

  # scales the height of image (300 px of height given per violin plot)
  scaler <- reactive({200*length(unique(scMuscle.seurat@meta.data[[variables()]]))})

  # generates image of split violin plot
  output$violin2 <- renderImage({

    # A temp file to save the output. It will be deleted after renderImage
    # sends it, because deleteFile=TRUE.
    outfile <- tempfile(fileext='.png')

    # Generate a png (1100 px width)
    png(outfile, width = 1100, height = scaler())
    print(
      VlnPlot(
        scMuscle.seurat,
        features = gene3(),
        group.by = splitviolincelltype(),
        pt.size = 0
      ) +
        NoLegend() +
        scale_y_continuous(expand=c(0,0.5))+
        facet_grid(rows = vars(scMuscle.seurat@meta.data[[variables()]])) +
        scale_colour_viridis_c() +
        vln.theme()
    )
    dev.off()

    # Return a list
    list(src = outfile,
         contentType = 'image/png',
         alt = "This is alternate text")
  }, deleteFile = TRUE)


  # generates sixth plot (All Cells - DotPlot)----
  output$dotplot <- renderPlot({
    input$action1
    isolate(
      DotPlot(
        scMuscle.seurat,
        features = dot(),
        group.by = variables()
      )+
        scale_colour_viridis_c() +
        dot.theme() +
        labs(title = "Metadata Features")
    )
  })

  # generates seventh plot (Myogenic Cells - PHATE) ----
  output$PHATE <- renderPlot({
    DimPlot(
      myo.slim.seurat,
      cells = sample(colnames(myo.slim.seurat)), #plot cells in random order
      reduction=reduction3(),
      group.by=bins(),
      cols=rainbow(25)%>%saturation(values=0.75) %>% brightness(values=0.9) %>% as.vector(), #adds colors for just the cell types present in this clustering
      na.value = NA, # removes noisy cells from plot
      pt.size = pt.size(), # see value above
      label.size = label.size(), # see value above
      repel = T,label= TRUE
    ) +
      NoLegend() +
      aes(stroke=pt.stroke())+
      xlab("PHATE_Harmony_1") +
      ylab("PHATE_Harmony_2") +
      umap.theme()
  })
  # generates eigth plot (Myogenic Cells - PHATE violins)----
  # scales the height of image (300 px of height given per violin plot)
  scaler2 <- reactive({200*length(gene4())})

  # generates image of PHATE violin plot
  output$phateviolin <- renderImage({
    input$action2
    isolate({
      # A temp file to save the output. It will be deleted after renderImage
      # sends it, because deleteFile=TRUE.
      outfile <- tempfile(fileext='.png')

      # Generate a png (1100 px width)
      png(outfile, width = 1100, height = scaler2())
      print(
        VlnPlot(
          myo.slim.seurat,
          features = gene4(),
          group.by = bins(),
          cols=rainbow(25)%>%saturation(values=0.75) %>% brightness(values=0.9) %>% as.vector(),
          combine=F,
          pt.size = 0
          ) %>% lapply(
            FUN = function(X) X +
              NoLegend() +
              scale_y_continuous(expand=c(0,.5))+
              scale_colour_viridis_c() +
              vln.theme()
            ) %>% wrap_plots(ncol=1)
      )
      dev.off()

      # Return a list
      list(src = outfile,
           contentType = 'image/png',
           alt = "This is alternate text")
    })}, deleteFile = TRUE)


  # DownloadHandler----
  # Allows plots to be donwloaded in specificed file type
  # All Cells - Cell Type UMAP----
  output$down1 <- downloadHandler(
    # specify file name
    filename = function() {
      paste("name", input$downloadable, sep = ".")
    },
    # creates the plot
    content = function(file) {
      if (input$downloadable == "pdf")
        pdf(file, width = 11, height = 8.5)
      else
        png(file, width = 11, height = 8.5, units = "in", res = 300)
      print(DimPlot(
        scMuscle.seurat,
        cells = sample(colnames(scMuscle.seurat)), #plot cells in random order
        reduction=reduction1(),
        group.by=reduction2(),
        cols=colors1(), #adds colors for just the cell types present in this clustering
        na.value = NA, # removes noisy cells from plot
        pt.size = pt.size(), # see value above
        label.size = label.size(), # see value above
        repel = T,label= TRUE
      ) +
        NoLegend() +
        aes(stroke=pt.stroke())+
        xlab(umapxlabel()) +
        ylab(umapylabel()) +
        umap.theme())
      dev.off()
    }
  )

  # All Cells - metadata UMAP----
  output$down2 <- downloadHandler(
    # specify file name
    filename = function() {
      paste("name", input$downloadable, sep = ".")
    },
    # creates the plot
    content = function(file) {
      if (input$downloadable == "pdf")
        pdf(file, width = 11, height = 8.5)
      else
        png(file, width = 11, height = 8.5, units = "in", res = 300)
      print(DimPlot(
        scMuscle.seurat,
        cells = sample(colnames(scMuscle.seurat)), #plot cells in random order
        reduction=reduction1(),
        group.by=variables(),
        cols=colors1(), #adds colors for just the cell types present in this clustering
        na.value = NA, # removes noisy cells from plot
        pt.size = pt.size(), # see value above
        label.size = label.size(), # see value above
        repel = T,label= TRUE
      ) +
        NoLegend() +
        aes(stroke=pt.stroke())+
        xlab(umapxlabel()) +
        ylab(umapylabel()) +
        umap.theme())
      dev.off()
    }
  )

  # All Cells - FeaturePlot/UMAP----
  output$down3 <- downloadHandler(
    # specify file name
    filename = function() {
      paste("name", input$downloadable, sep = ".")
    },
    # creates the plot
    content = function(file) {
      if (input$downloadable == "pdf")
        pdf(file, width = 11, height = 8.5)
      else
        png(file, width = 11, height = 8.5, units = "in", res = 300)
      print(FeaturePlot(scMuscle.seurat, features = gene1(), reduction = reduction1())+scale_colour_viridis_c()+umap.theme())
      dev.off()
    }
  )

  # All Cells - Violin Plot----
  output$down4 <- downloadHandler(
    # specify file name
    filename = function() {
      paste("name", input$downloadable, sep = ".")
    },
    # creates the plot
    content = function(file) {
      if (input$downloadable == "pdf")
        pdf(file, width = 11, height = 8.5)
      else
        png(file, width = 11, height = 8.5, units = "in", res = 300)
      print(VlnPlot(scMuscle.seurat, features = gene2(), group.by = variables(), pt.size = 0) +
              NoLegend() + scale_y_continuous(expand=c(0,0))+scale_colour_viridis_c() + vln.theme())
      dev.off()
    }
  )

  # All Cells - Split Violin Plot----
  output$down5 <- downloadHandler(
    # specify file name
    filename = function() {
      paste("name", input$downloadable, sep = ".")
    },
    # creates the plot
    content = function(file) {
      if (input$downloadable == "pdf")
        pdf(file, width = 11, height = scaler()/100)
      else
        png(file, width = 11, height = scaler()/100, units = "in", res = 300)
      print(VlnPlot(scMuscle.seurat, features = gene3(), group.by = splitviolincelltype(), pt.size = 0) + NoLegend() +
              scale_y_continuous(expand=c(0,.5))+facet_grid(rows = vars(scMuscle.seurat@meta.data[[variables()]]))
            +scale_colour_viridis_c()+ vln.theme())
      dev.off()
    }
  )
  # All Cells - DotPlot----
  output$down6 <- downloadHandler(
    # specify file name
    filename = function() {
      paste("name", input$downloadable, sep = ".")
    },
    # creates the plot
    content = function(file) {
      if (input$downloadable == "pdf")
        pdf(file, width = 11, height = 8.5)
      else
        png(file, width = 11, height = 8.5, units = "in", res = 300)
      print(DotPlot(scMuscle.seurat, features = dot(), group.by = variables())
            +scale_colour_viridis_c()+ dot.theme())
      dev.off()
    }
  )

  
  # 
  # Myogenic Cells - PHATE----
  output$down7 <- downloadHandler(
    # specify file name
    filename = function() {
      paste("name", input$downloadable, sep = ".")
    },
    # creates the plot
    content = function(file) {
      if (input$downloadable == "pdf")
        pdf(file, width = 11, height = 8.5)
      else
        png(file, width = 11, height = 8.5, units = "in", res = 300)
      print(DimPlot(
        myo.slim.seurat,
        cells = sample(colnames(myo.slim.seurat)), #plot cells in random order
        reduction=reduction3(),
        group.by=bins(),
        cols=rainbow(25)%>%saturation(values=0.75) %>% brightness(values=0.9) %>% as.vector(), #adds colors for just the cell types present in this clustering
        na.value = NA, # removes noisy cells from plot
        pt.size = pt.size(), # see value above
        label.size = label.size(), # see value above
        repel = T,label= TRUE
      ) +
        NoLegend() +
        aes(stroke=pt.stroke())+
        xlab("PHATE_Harmony_1") +
        ylab("PHATE_Harmony_2") +
        umap.theme())
      dev.off()
    }
  )
  
  # Myogenic Cells - PHATE violins----
  output$down8 <- downloadHandler(
    # specify file name
    filename = function() {
      paste("name", input$downloadable, sep = ".")
    },
    # creates the plot
    content = function(file) {
      if (input$downloadable == "pdf")
        pdf(file, width = 11, height = 8.5)
      else
        png(file, width = 1100, height = scaler2(), units = "px", res = 300)
      print(VlnPlot(
        myo.slim.seurat,
        features = gene4(),
        group.by = bins(),
        cols=rainbow(25)%>%saturation(values=0.75) %>% brightness(values=0.9) %>% as.vector(),
        combine=F,
        pt.size = 0
      ) %>% lapply(
        FUN = function(X) X +
          NoLegend() +
          scale_y_continuous(expand=c(0,.5))+
          scale_colour_viridis_c() +
          vln.theme()
      ) %>% wrap_plots(ncol=1))
      dev.off()
    }
  )
}

# Run app ----
shinyApp(ui, server)
