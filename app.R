library(shiny)
library(deSolve)
library(ggplot2)
library(reshape2)

# =========================
# Lotka–Volterra model
# =========================
phytoplankton_model <- function(t, state, parameters) {
  
  P1 <- state[1]
  P2 <- state[2]
  
  with(as.list(parameters), {
    
    # species 2 introduction
    if (t < t_intro2) {
      P2 <- 0
    }
    
    dP1_dt <- r1 * P1 * (1 - (alpha11*P1 + alpha12 * P2))  # Logistic growth with competition for species 1
    dP2_dt <- r2 * P2 * (1 - (alpha22*P2 + alpha21 * P1))  # Logistic growth with competition for species 2

    # Alternative formulation with carrying capacities:
    # dP1_dt <- r1 * P1 * (1 - (P1 + alpha12 * P2) / K1)
    # dP2_dt <- r2 * P2 * (1 - (P2 + alpha21 * P1) / K2)
    
    list(c(dP1_dt, dP2_dt))
  })
}

# =========================
# UI
# =========================
ui <- fluidPage(
  
  titlePanel("Lotka–Volterra Competition: First-Come–First-Served"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      h4("Initial conditions"),
      sliderInput("P1_0", "Initial biomass P1", min = 0.1, max = 50, value = 10),
      sliderInput("P2_0", "Initial biomass P2", min = 0.1, max = 50, value = 10),
      
      h4("Growth rates"),
      sliderInput("r1", "r1", min = 0.1, max = 2, value = 0.8, step = 0.05),
      sliderInput("r2", "r2", min = 0.1, max = 2, value = 0.8, step = 0.05),
      
      h4("intraspecific Competition coefficients"),
      sliderInput("alpha11", "α11 (effect of P1 on P1)", min = 0, max = 3, value = 1, step = 0.05),
      sliderInput("alpha22", "α22 (effect of P2 on P2)", min = 0, max = 3, value = 1, step = 0.05),

      # h4("Carrying capacities"),
      # sliderInput("K1", "K1", min = 10, max = 300, value = 100),
      # sliderInput("K2", "K2", min = 10, max = 300, value = 90),
      
      h4("interspecific Competition coefficients"),
      sliderInput("alpha12", "α12 (effect of P2 on P1)", min = 0, max = 3, value = 1, step = 0.05),
      sliderInput("alpha21", "α21 (effect of P1 on P2)", min = 0, max = 3, value = 1, step = 0.05),
      
      h4("Species 2 introduction"),
      sliderInput("t_intro2", "Introduction time of species 2", min = 0, max = 50, value = 0, step = 0.5),
      
      h4("Simulation"),
      sliderInput("tmax", "Simulation length", min = 100, max = 2000, value = 1000, step = 100)
    ),
    
    mainPanel(
      plotOutput("timeseries", height = "400px"),
      plotOutput("dominance", height = "250px"),
      verbatimTextOutput("summary")
    )
  )
)

# =========================
# Server
# =========================
server <- function(input, output) {
  
  simulation <- reactive({
    
    state <- c(
      P1 = input$P1_0,
      P2 = input$P2_0
    )
    
    parameters <- c(
      r1 = input$r1,
      r2 = input$r2,
      alpha12 = input$alpha12,
      alpha21 = input$alpha21,
      alpha11 = input$alpha11,
      alpha22 = input$alpha22,
      
      #alternative formulation with carrying capacities
      # K1 = input$K1,
      # K2 = input$K2,
      t_intro2 = input$t_intro2
    )
    
    times <- seq(0, input$tmax, by = 1)
    
    out <- ode(
      y = state,
      times = times,
      func = phytoplankton_model,
      parms = parameters
    )
    
    as.data.frame(out)
  })
  
  output$timeseries <- renderPlot({
    
    df <- simulation()
    df_long <- melt(df, id.vars = "time",
                    variable.name = "Species",
                    value.name = "Biomass")
    
    ggplot(df_long, aes(time, Biomass, color = Species)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = c("P1" = "red", "P2" = "darkgreen")) +
      theme_bw() +
      labs(y = "Biomass", x = "Time")
  })
  
  output$dominance <- renderPlot({
    
    df <- simulation()
    
    df$P1_frac <- df$P1 / (df$P1 + df$P2)
    df$P2_frac <- df$P2 / (df$P1 + df$P2)
    
    ggplot(df, aes(time, P1_frac)) +
      geom_line(data=df, aes(time, P1_frac), color = "red", linewidth = 1) +
      geom_line(data=df, aes(time, P2_frac), color = "darkgreen", linewidth = 1) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
      ylim(0, 1) +
      theme_bw() +
      labs(y = "Fraction of species 1", x = "Time")
  })
  
  output$summary <- renderPrint({
    
    df <- simulation()
    
    P1_end <- tail(df$P1, 1)
    P2_end <- tail(df$P2, 1)
    
    cat("Final biomasses:\n")
    cat("P1 =", round(P1_end, 2), "\n")
    cat("P2 =", round(P2_end, 2), "\n\n")
    
    cat("Final dominance (% P1):\n")
    cat(round(100 * P1_end / (P1_end + P2_end), 1), "%\n")
  })
}

# =========================
# Run app
# =========================
shinyApp(ui = ui, server = server)
