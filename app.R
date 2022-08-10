library(shiny)
library(DT)
#library(dplyr)


ui <- fluidPage(

    HTML(r"(<p style="text-align:center;">Joey Mays - Updated 2022-08-10</p>)"),
    
    titlePanel("NGS Library Pooling Calculator"),

    sidebarLayout(
        sidebarPanel(width = 3,
            textInput("sample.name", "Sample Names", placeholder = "paste here"),
            textInput("input.conc", "Concentration (ng/uL)", placeholder = "paste here"),
            textInput("n.req.reads", "Reads Requested (M)", placeholder = "paste here"),
            h5("Set all average library sizes:"),
            fluidRow(
                column(6,
                       actionButton("avg.size.320", "320 bp", class = "btn-primary")),
                column(6, 
                       actionButton("avg.size.1000", "1000 bp", class = "btn-primary"))
            ),
            p(" "),
            numericInput("final.vol", "Pool Volume (uL)", value = 50),
            downloadButton("TableDownload", label = "Download Results", class = "btn-success"),
            hr(),
            h3("Instructions"),
            p("1. Copy sample names from spreadsheet and paste in box."),
            p("2. Copy and paste concentrations and requested reads."),
            p("3. Click buttons to set avg library size."),
            p("4. Edit cells as necessary by double clicking. Click outside the table to save.")
        ),

        mainPanel(
           DT::dataTableOutput("input.data.table"),
           textOutput("n.reads"),
           textOutput("pooling.conc")
        )
    )
)

server <- function(input, output) {

    input.data.table <- reactiveVal(
        data.frame(matrix(data = NA, ncol = 7, 
                          dimnames = list(NULL, c("sample.name", "input.conc", "n.req.reads", "relative.representation", "avg.lib.size", "calc.conc", "amt.to.pool.prelim"))))
    )
    
    pooling.conc <- reactiveVal(1)
    
    total.reads <- reactive({
        sum(input.data.table()$n.req.reads)
    })
    
    toListen <- reactive({
        list(input$final.vol, input.data.table())
    })
    
    observeEvent(toListen(), {
        req(input$sample.name)
        print(input.data.table())
        temp.input <- input.data.table()

        if(sum(temp.input$avg.lib.size) != 0){

        temp.input$calc.conc <- temp.input$input.conc/1000/660*10^6/temp.input$avg.lib.size*10^6
        }
        print(sum(temp.input$calc.conc))
        if(sum(temp.input$calc.conc) != 0){
            amt.to.pool.prelim <- (temp.input$relative.representation*input$final.vol)/(temp.input$calc.conc)
            amt.to.pool.prelim.sum <- sum(amt.to.pool.prelim) / 100
            print(paste0("amt to pool prelim sum ", amt.to.pool.prelim.sum))
            pooling.conc(input$final.vol / amt.to.pool.prelim.sum)
            temp.input$amt.to.pool.prelim <- amt.to.pool.prelim * pooling.conc() * 0.01
        }
        input.data.table(temp.input)
    })

    output$n.reads <- renderText(paste0("Total Number of Reads: ", total.reads(), "M"))
    
    output$pooling.conc <- renderText(paste0("Pooling Concentration (pM): ", format(round(pooling.conc(),0), nsmall = 0, big.mark=",")))
    
    observeEvent(input$sample.name, {
        req(input$sample.name)
        raw.paste <- input$sample.name
        raw.paste <- strsplit(raw.paste, "\\s")[[1]]
        input.data.table(data.frame(sample.name = raw.paste, input.conc = 0, n.req.reads = 0, relative.representation = 0, avg.lib.size = 0, calc.conc = 0, amt.to.pool.prelim = 0))
    })

    observeEvent(input$input.conc, {
        req(input$sample.name)
        raw.paste <- input$input.conc
        raw.paste <- strsplit(raw.paste, "\\s")[[1]]
        temp.input <- input.data.table()
        temp.input$input.conc <- as.numeric(raw.paste)
        input.data.table(temp.input)
    })
    
    observeEvent(input$n.req.reads, {
        req(input$sample.name)
        raw.paste <- input$n.req.reads
        raw.paste <- strsplit(raw.paste, "\\s")[[1]]
        temp.input <- input.data.table()
        temp.input$n.req.reads <- as.numeric(raw.paste)
        temp.input$relative.representation <- round(100 * temp.input$n.req.reads / sum(temp.input$n.req.reads),2)
        input.data.table(temp.input)
    })
    
    observeEvent(input$avg.size.320, {
        req(input$sample.name)
        temp.input <- input.data.table()
        temp.input$avg.lib.size <- 320
        input.data.table(temp.input)
    })
    
    observeEvent(input$avg.size.1000, {
        req(input$sample.name)
        temp.input <- input.data.table()
        temp.input$avg.lib.size <- 1000
        input.data.table(temp.input)
    })
    
    output$input.data.table <- DT::renderDataTable(
        datatable(input.data.table(), 
                  options = list(ordering=F, dom="t"), 
                  editable = T, rownames = F, 
                  colnames = c("Sample Name", "Concentration (ng/uL)", "Reads Requested (M)", "Relative Representation (%)","Avg Library Size (bp)", "Concentration (pM)", "Amount to Pool (ul)")) %>% 
            formatRound(columns = c(2,3,4,6,7), digits = c(2, 0, 2, 0, 2)))
    
    observeEvent(input$input.data.table_cell_edit, {
        input.data.table(editData(input.data.table(), input$input.data.table_cell_edit, rownames = F)) 
    })
    
    
    
    output$TableDownload <- downloadHandler(
        filename = function() {
            paste0("ngs-pooling_",as.character(Sys.Date()),".txt")
        },
        content = function(file) {
            write.table(input.data.table(), file, quote = F, row.names = F, col.names = c("Sample Name", "Concentration (ng/uL)", "Reads Requested (M)", "Relative Representation (%)","Avg Library Size (bp)", "Concentration (pM)", "Amount to Pool (ul)"), sep = '\t')
        }
    )
    
}

# Run the application 
shinyApp(ui = ui, server = server)
