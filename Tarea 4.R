install.packages("biostrings")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
library("msa")
library(BiocGenerics)
library(BiocManager)
library(Biostrings)
library("msa")

install.packages("deSolve")
install.packages("ggplot2")
install.packages("shiny")
library(shiny)
library(deSolve)
library(ggplot2)
library(bslib)

# Define the models
si_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I
    dI <- beta * S * I
    return(list(c(dS, dI)))
  })
}

sis_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I + gamma * I
    dI <- beta * S * I - gamma * I
    return(list(c(dS, dI)))
  })
}

sir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I
    dI <- beta * S * I - gamma * I
    dR <- gamma * I
    return(list(c(dS, dI, dR)))
  })
}

# UI
ui <- page_navbar(
  title = "Modelos básicos por compartimentos",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  nav_panel("SI", 
            layout_sidebar(
              sidebar = sidebar(
                sliderInput("si_beta", "Transmission rate (β)", 0.001, 1, 0.3, step = 0.001),
                sliderInput("si_I0", "Proporción inicial infectada", 0.001, 0.1, 0.01, step = 0.001),
                sliderInput("si_tmax", "Rango de tiempo", 1, 100, 50)
              ),
              plotOutput("si_plot")
            )
  ),
  
  nav_panel("SIS ",
            layout_sidebar(
              sidebar = sidebar(
                sliderInput("sis_beta", "Transmission rate (β)", 0.001, 1, 0.3, step = 0.001),
                sliderInput("sis_gamma", "Tasa de recuperación (γ)", 0.001, 1, 0.1, step = 0.001),
                sliderInput("sis_I0", "Proporción inicial infectada", 0.001, 0.1, 0.01, step = 0.001),
                sliderInput("sis_tmax", "Rango de tiempo", 1, 100, 50)
              ),
              plotOutput("sis_plot")
            )
  ),
  
  nav_panel("SIR ",
            layout_sidebar(
              sidebar = sidebar(
                sliderInput("sir_beta", "Transmission rate (β)", 0.001, 1, 0.3, step = 0.001),
                sliderInput("sir_gamma", "Tasa de recuperación (γ)", 0.001, 1, 0.1, step = 0.001),
                sliderInput("sir_I0", "Proporción inicial infectada", 0.001, 0.1, 0.01, step = 0.001),
                sliderInput("sir_tmax", "Rango de tiempo", 1, 100, 50)
              ),
              plotOutput("sir_plot")
            )
  ),
  
  nav_panel("Ecuaciones diferenciales",
            card(
              card_header("Ecuaciones de los modelos"),
              card_body(
                h4("Modelo SI :"),
                withMathJax("$$\\frac{dS}{dt} = -\\beta S I$$"),
                withMathJax("$$\\frac{dI}{dt} = \\beta S I$$"),
                hr(),
                h4("Model SIS:"),
                withMathJax("$$\\frac{dS}{dt} = -\\beta S I + \\gamma I$$"),
                withMathJax("$$\\frac{dI}{dt} = \\beta S I - \\gamma I$$"),
                hr(),
                h4("Model SIR :"),
                withMathJax("$$\\frac{dS}{dt} = -\\beta S I$$"),
                withMathJax("$$\\frac{dI}{dt} = \\beta S I - \\gamma I$$"),
                withMathJax("$$\\frac{dR}{dt} = \\gamma I$$")
              )
            )
  )
)

# Server
server <- function(input, output) {
  
  output$si_plot <- renderPlot({
    parameters <- c(beta = input$si_beta)
    initial_state <- c(S = 1 - input$si_I0, I = input$si_I0)
    times <- seq(0, input$si_tmax, by = 0.1)
    
    out <- ode(y = initial_state, times = times, func = si_model, parms = parameters)
    df <- as.data.frame(out)
    
    ggplot(df, aes(x = time)) +
      geom_line(aes(y = S, color = "Susceptibles")) +
      geom_line(aes(y = I, color = "Infectados")) +
      labs(title = "SI Model", x = "Tiempo", y = "Proporción", color = "Compartimento") +
      theme_minimal() +
      scale_color_manual(values = c("Susceptibles" = "blue", "Infectados" = "red"))
  })
  
  output$sis_plot <- renderPlot({
    parameters <- c(beta = input$sis_beta, gamma = input$sis_gamma)
    initial_state <- c(S = 1 - input$sis_I0, I = input$sis_I0)
    times <- seq(0, input$sis_tmax, by = 0.1)
    
    out <- ode(y = initial_state, times = times, func = sis_model, parms = parameters)
    df <- as.data.frame(out)
    
    ggplot(df, aes(x = time)) +
      geom_line(aes(y = S, color = "Susceptibles")) +
      geom_line(aes(y = I, color = "Infectados")) +
      labs(title = "SIS ", x = "Tiempo", y = "Proporción", color = "Compartimento") +
      theme_minimal() +
      scale_color_manual(values = c("Susceptibles" = "blue", "Infectados" = "red"))
  })
  
  output$sir_plot <- renderPlot({
    parameters <- c(beta = input$sir_beta, gamma = input$sir_gamma)
    initial_state <- c(S = 1 - input$sir_I0, I = input$sir_I0, R = 0)
    times <- seq(0, input$sir_tmax, by = 0.1)
    
    out <- ode(y = initial_state, times = times, func = sir_model, parms = parameters)
    df <- as.data.frame(out)
    
    ggplot(df, aes(x = time)) +
      geom_line(aes(y = S, color = "Susceptibles")) +
      geom_line(aes(y = I, color = "Infectados")) +
      geom_line(aes(y = R, color = "Recuperados")) +
      labs(title = "SIR Model", x = "Tiempo", y = "Proporción", color = "Compartimento") +
      theme_minimal() +
      scale_color_manual(values = c("Susceptibles" = "blue", "Infectados" = "red", "Recuperados" = "green"))
  })
}

# Run the app
shinyApp(ui = ui, server = server)


library(tidyverse)
library(patchwork)
library(RColorBrewer)


a <- data.frame(
  x = c("a", "b", "c", "d", "e"),
  y = c(10, 20, 30, 450, 500)
)

lower <- a %>% 
  ggplot(aes(x = x, y = y)) +
  geom_bar(stat = "identity", aes(fill = x),
           width = 0.7, color = "black") +
  geom_text(aes(label = y), vjust = -0.2) +
  scale_fill_manual(values = brewer.pal(8, "Set2")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.y = element_text(hjust =1 )) +
  coord_cartesian(ylim = c(0, 50))

upper1 <- a %>% 
  ggplot(aes(x = x, y = y)) +
  geom_bar(stat = "identity", aes(fill = x),
           width = 0.7, color = "black") +
  geom_text(aes(label = y), vjust = -0.2) +
  scale_fill_manual(values = brewer.pal(8, "Set2")) +
  labs(x = NULL,
       y = NULL) +
  theme_classic() +
  theme(legend.position = "none",
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  coord_cartesian(ylim = c(400, 550))

upper2 <- a %>% 
  ggplot(aes(x = x, y = y)) +
  geom_bar(stat = "identity", aes(fill = x),
           width = 0.7, color = "black") +
  geom_text(aes(label = y), vjust = -0.2) +
  scale_fill_manual(values = brewer.pal(8, "Set2")) +
  labs(x = NULL,
       y = NULL) +
  theme_classic() +
  theme(legend.position = "none",
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  coord_cartesian(ylim = c(450, 520))


broken1 <- wrap_plots(upper1, lower, nrow = 2) &
  theme(axis.text = element_text(color = "black"),
        text = element_text(size = 14))

broken2 <- wrap_plots(upper2, lower, nrow = 2) &
  theme(axis.text = element_text(color = "black"),
        text = element_text(size = 14))

wrap_plots(broken1, broken2, nrow = 1)

ggsave("../Results/Broken_axis.svg", height = 4, width = 6, bg = "white")
ggsave("../Results/Broken_axis.png", height = 4, width = 6, bg = "white")




left <- a %>% 
  ggplot(aes(x = y, y = x)) +
  geom_bar(stat = "identity", aes(fill = x),
           width = 0.7, color = "black") +
  geom_text(aes(label = y), hjust = -0.2) +
  scale_fill_manual(values = brewer.pal(8, "Set2")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_text(hjust =1 )) +
  coord_cartesian(xlim = c(0, 50))


right <- a %>% 
  ggplot(aes(x = y, y = x)) +
  geom_bar(stat = "identity", aes(fill = x),
           width = 0.7, color = "black") +
  geom_text(aes(label = y),hjust = -0.2) +
  scale_fill_manual(values = brewer.pal(8, "Set2")) +
  labs(x = NULL,
       y = NULL) +
  theme_classic() +
  theme(legend.position = "none",
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  coord_cartesian(xlim = c(400, 600))

wrap_plots(left, right, nrow = 1) &
  theme(axis.text = element_text(color = "black"),
        text = element_text(size = 14))




Corynebacterium<-readAAStringSet("C:/Users/nat_g/OneDrive/Escritorio/Bioinfo/Archivo de secuencias2.txt")
Corynebacterium
letterFrequency(Corynebacterium,"A")
letterFrequency(Corynebacterium, "D")
letterFrequency(Corynobacterium, "E")
letterFrequency(Corynebacterium, "E")
letterFrequency(Corynebacterium, "R")
letterFrequency(Corynebacterium, "K")
letterFrequency(Corynebacterium, "N")
letterFrequency(Corynebacterium, "H")
letterFrequency(Corynebacterium, "Q")
letterFrequency(Corynebacterium, "S")
letterFrequency(Corynebacterium, "T")
letterFrequency(Corynebacterium, "G")
letterFrequency(Corynebacterium, "V")
letterFrequency(Corynebacterium, "P")
letterFrequency(Corynebacterium, "L")
letterFrequency(Corynebacterium, "F")
letterFrequency(Corynebacterium, "Y")
letterFrequency(Corynebacterium, "I")
letterFrequency(Corynebacterium, "M")
letterFrequency(Corynebacterium, "W")
letterFrequency(Corynebacterium, "C")

A<-letterFrequency(Corynebacterium,"A")
A
D<-letterFrequency(Corynebacterium, "D")
E<-letterFrequency(Corynebacterium, "E")
R<-letterFrequency(Corynebacterium, "R")
K<-letterFrequency(Corynebacterium, "K")
N<-letterFrequency(Corynebacterium, "N")
H<-letterFrequency(Corynebacterium, "H")
Q<-letterFrequency(Corynebacterium, "Q")
S<-letterFrequency(Corynebacterium, "S")
Tr<-letterFrequency(Corynebacterium, "T")
G<-letterFrequency(Corynebacterium, "G")
V<-letterFrequency(Corynebacterium, "V")
P<-letterFrequency(Corynebacterium, "P")
L<-letterFrequency(Corynebacterium, "L")
Fe<-letterFrequency(Corynebacterium, "F")
Y<-letterFrequency(Corynebacterium, "Y")
I<-letterFrequency(Corynebacterium, "I")
M<-letterFrequency(Corynebacterium, "M")
W<-letterFrequency(Corynebacterium, "W")
C<-letterFrequency(Corynebacterium, "C")
library("plyr")
install.packages("plyr")
barplot("A")
barplot(A)
A
barplot(A, main = "Frequencia de Alanina",
        col = rainbow(3))
barplot(A, main = "Frequencia de Alanina",
        col = aquamarine)
hist(A, main = "Frequencia de Alanina",
     col = aquamarine)
hist(A, main = "Frequencia de Alanina",
     col = rainbow(3))
A
boxplot.default(A, main = "Frequencia de Alanina",
                col = rainbow(3))
boxplot(A, main = "Frequencia de Alanina",
        col = rainbow(3))
plot(A, main = "Frequencia de Alanina",
     col = rainbow(3))
Sq1<-as.character(Corynebacterium[[1]])
sq2<-as.character(Corynebacterium[[2]])
PairwiseAlignments(pattern = sq2, subject = Sq1)

library(outbreaks)
library(incidence)
library(epiR)
library(ggplot2)
#con error relativo, usamos la función epi.sssimpleestb
relativo<-epi.sssimpleestb(N = NA, Py = 0.35,epsilon = 0.035, error = "relative",se = 1,sp = 1,conf.level = 0.95)#el valor de epsilon lo cambiamos a 0.035 cuando queremos hacerlo con epsilon de la proporcion de 10% de 35%
relativo
### ahora con error absoluto
absoluto<-epi.sssimpleestb(N = NA, Py = 0.35,epsilon = 0.035, error = "absolute",se = 1,sp = 1,conf.level = 0.95)
absoluto
#

###prueba de detección
relativo2<-epi.sssimpleestb(N=NA,Py = .35,epsilon = .1,error = "relative", se = 0.96,sp = 0.85,conf.level = 0.95)
absoluto2<-epi.sssimpleestb(N=NA,Py = .35,epsilon = .1,error = "absolute", se = 0.96,sp = 0.85,conf.level = 0.95)
relativo2
absoluto2
#población de riesgo
riesgo3k<-epi.sssimpleestb(N = 3000,Py = 0.35,epsilon = .1,error = "absolute",se = 1, sp = 1,conf.level = .95)
riesgo300<-epi.sssimpleestb(N = 300,Py = 0.35,epsilon = .1,error = "absolute",se = 1, sp = 1,conf.level = .95)
riesgo3k
riesgo300

########ejercicio2
tamaño<-epi.sscc(OR = 2, p1 = .3,p0 = .3,n = NA,power = 0.8,conf.level = 0.95, sided.test = 2)
tamaño

#ejercicio 3
#parte 1
epi.sscohortt(FT = 6,irexp1 = 0.005,irexp0 = 0.007,conf.level = 0.95,n = NA,r = 1,power = 0.8)
#Ejercicio 3 500 gatos
epi.sscohortt(FT = 6,irexp1 = 0.005,irexp0 = 0.007,conf.level = 0.95,n = 500,r = 1,power = NA)

