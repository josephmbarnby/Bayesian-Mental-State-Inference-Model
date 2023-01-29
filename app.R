
#install.packages(c('ggplot2','dplyr', 'shiny'))
library(shiny)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(tidyquant)
source('gen_ut.R')


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel(h1("Shallow DoM Attribution Model")),
    titlePanel(h4(HTML(paste(
                    "This is an interactive environment to simulate attributional outcomes of a social observer during an modified repeated Dictator game.",
                  "<br/>",
                  sep="<br/>")))),

    titlePanel(h4("This agent uses the following equations:")),
    br(),
    withMathJax(),
      tabPanel(
      title = "Diagnostics",
      h4(textOutput("diagTitle")),
      uiOutput("formula")
      ),
    br(),
    titlePanel(h4(HTML(paste("Move the sliders to change the task structure, and the agent's and the subject's policy.",
                  "<br/>",
                  "The red line is simulated harmful intent attributions emitted over all trials,
                  and the blue line is the simulated self interest attributions emmited over all trials.",
                  "</br/>",
                  "Click 'Select a new agent' to start a new agent from scratch using the same settings",
                  sep="<br/>")))),
    br(),
    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            actionButton("setseed", "Select a new agent"),
            br(),
            br(),
            sliderInput("pHI0",
                        HTML(paste("Prior magnitude over HI: pHI0")),
                        min = 0.1,
                        max = 1,
                        value = 0.5,
                        step = 0.1),
            br(),
            #sliderInput("uHI0",
            #            HTML(paste("Prior uncertainty over HI: uHI0")),
            #            min = 0.01,
            #            max = 5,
            #            value = 1,
            #            step = 0.1),
            #br(),
            sliderInput("pSI0",
                        HTML(paste("Prior magnitude over SI: pSI0")),
                        min = 0.1,
                        max = 1,
                        value = 0.5,
                        step = 0.1),
            br(),
            sliderInput("uSIHI0",
                        HTML(paste("Uncertainty over priors: uPri")),
                        min = 0.1,
                        max = 5,
                        value = 1,
                        step = 0.1),
            br(),
            sliderInput("upi",
                        HTML(paste("Posterior Integraton Noise: u&pi;")),
                        min = 0.1,
                        max = 10,
                        value = 1,
                        step = 0.1),
            br(),
            #sliderInput("w0",
            #            HTML(paste("Likelihood Intercept: w0")),
            #            min = -5,
            #            max = 5,
            #            value = 0,
            #            step = 0.5),
            #br(),
            sliderInput("wHI",
                        HTML(paste("Ratio of HI and SI updating: wHI:wSI")),
                        min = 0.5,
                        max = 2,
                        value = 1,
                        step = 0.1),
            br(),
            #sliderInput("wSI",
            #            HTML(paste("Sensitivity to update SI: wSI")),
            #            min = 0.1,
            #            max = 5,
            #            value = 0.5,
            #            step = 0.1),
            #br(),
            #sliderInput("eta",
            #            HTML(paste("Transfer of policies between blocks: &eta;")),
            #            min = 0.1,
            #            max = 1,
            #            value = 1,
            #            step = 0.1),
            #br(),
            sliderInput("trials",
                        "Task Length:",
                        min = 10,
                        max = 100,
                        value = 50,
                        step = 10),
            br(),
            #sliderInput("partner",
            #            "Block Length:",
            #            min = 1,
            #            max = 4,
            #            value = 1,
            #            step = 1),
            #br(),
            sliderInput("partner_decisions",
                        "p(Decision = Unfair):",
                        min = 0.1,
                        max = 1,
                        value = 0.5,
                        step = 0.1),
            br(),
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot", height = '900px', width = 'auto')
        )
    ),

    titlePanel(h5(HTML(paste("CC <a href='https://www.joebarnby.com/'>Dr Joe M Barnby</a> 2022"))))

)

# Define server logic required to draw a histogram
server <- function(input, output) {

    #output$formula <- renderUI({
    #withMathJax(paste0("
    #                   $$Q^{t}_{c} = Q^{t-1}_{c} * \\lambda + ({Reward - Q^{t-1}_{c}})$$
    #                   $$p(\\hat{c} = c) = \\frac{e^{\\frac{Q^{t}_{c}}{\\tau}}}{\\sum_{c'\\in(c_1, c_2)} e^{\\frac{Q^{t}_{c'}}{\\tau}}}$$
#
    #                   $$\\text{Therefore, } Q^{t}_{c} = \\text{the internal beliefs the agent holds about the value of each card at each trial}$$
    #                   "))
    #})

    output$distPlot <- renderPlot({

        phase   <- 1#input$partner;  # number or partners or phases
        tn      <- input$trials
        tot_t   <- phase*tn

        Nb = 9;
        Na = 2;
        err= 0.02/(Nb*Nb)

        PSI0 <- noisyBino(input$pSI0,
                          input$uSIHI0,
                          Nb);
        PHI0 <- noisyBino(input$pHI0,
                          input$uSIHI0,
                          Nb);

        PSIHI0 <- PSI0 %*% t(PHI0);

        # generate bins based on input$bins from ui.R

        upi     <- input$upi;
        w0      <- 0#input$w0;
        whi     <- input$wHI;
        wsi     <- 1/input$wHI;
        #eta     <- input$eta;
        seed    <- input$setseed;

        pi      <- array(NA,c(Nb,Nb,Na));

        offs = (Nb+1)/2
        for (SI in 1:Nb){
          for (HI in 1:Nb){
            pi[SI,HI,1]  = invlogit(w0  + (wsi*(SI-offs)) + (whi*(HI-offs)))
            pi[SI,HI,2]  = 1 - pi[SI,HI,1]
          }
        }

        pri0 <- PSIHI0;
        post <- pri0; # this is the belief that will be updated with each trial

        # rows to be processed and convenience copies of other's actions,
        # malice and greed attributions:

        llout               <- list()
        llout[[1]]          <- matrix(NA,tot_t+1,6);
        llout[[1]][,1]      <- c(0,1:(tn*phase))
        colnames(llout[[1]])<- c('trial','ret','HImode','SImode','HIsim','SIsim')
        pol0                <- pri0^(1/upi)
        pol0                <- pol0 / sum(as.vector(pol0))
        pol0                <- (pol0+err)/(1+err*length(pol0))
        llout[[2]]          <- array(dim=c(dim(pri0),(1+tn*phase)))
        llout[[2]][,,1]     <- pri0;
        names(llout)        <- c('evo','policy')

        for (t in 1:tot_t){  # loop

          #if (t == tn+1 | t == (tn*2)+1 | t == (tn*3)+1 | t == (tn*4)+1){
          #  post = (pri0 * (1-eta)) + (post * eta);
          #}

          outcome            <- sample(c(0, 0.5), 1, prob = c(input$partner_decisions, 1-input$partner_decisions))
          llout[[1]][(t+1),2]<- outcome
          aind               <- round((Na-1)*outcome+1)

          pri                <- post;
          post               <- pi[,,aind] * pri
          post               <- post / sum(as.vector(post))

          # Now the probability of the response
          pol                <- post^(1/upi);
          pol                <- pol/sum(as.vector(pol)); #renormalise
          pol                <- (pol+err)/(1+err*length(pol)) # add noise floor

          # find mode of the policy
          c=max.col(pol);  # this finds the max. col. pos. for each row
          m=rep(NA, Nb);
          for (x in 1:Nb){m[x]=pol[x,c[x]]};  # retrieve the values ...
          r=which.max(m);  # ... so was to find the row with the mode of the distro.
          llout$evo[t+1,c('HImode','SImode')] <- c(c[r]-0.5,r-0.5)/Nb
          # Now sample randomly
          hisim <- rmdf(1,colSums(pol)) # joint=marginal*conditional, so sample first dim from the marginal ...
          sisim <- rmdf(1,pol[,hisim]/sum(pol[,hisim]))  # ... and the phase from the corresponding conditional.
          llout$evo[(t+1),c('HIsim','SIsim')] <- c((hisim-0.5),(sisim-0.5))/Nb
          llout$policy[,,(t+1)] <- pol

        }


        #Set colours

        shift_df <- data.frame(HIstart = colSums(llout$policy[,,1]),
                  SIstart = rowSums(llout$policy[,,1]),
                  HIend = colSums(llout$policy[,,(tn*phase)+1]),
                  SIend = rowSums(llout$policy[,,(tn*phase)+1]),
                  res = 1:9)
        policy_change <- ggplot(shift_df %>%
                                  pivot_longer(HIstart:SIend, 'policy', values_to = 'probability') %>%
                                  mutate(facet = c(rep(c(1,2,1,2), 9))),
                                aes(res/10, probability, fill = policy))+
          geom_col(position = position_dodge())+
          #geom_density(aes(after_stat(density)))+
          labs(x = 'Attribution', y = 'Probability', title = 'Policy Change from Trial 1 to N')+
          scale_fill_manual(values = c('#A31621', '#DB222A', '#1F7A8C', '#BFDBF7'), name = "")+
          facet_wrap(~facet)+
          theme_tq()+
          theme(plot.background = element_blank(),
                legend.position = 'bottom',
                legend.direction = 'horizontal',
                strip.background = element_blank())
        initial_policy <- pi[,,1] %>%
          as.data.frame() %>%
          mutate(SI = c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85)) %>%
          pivot_longer(1:9, 'HI', values_to = 'Probability') %>%
          ggplot(aes(HI, SI, fill = Probability)) +
          geom_tile() +
          scale_fill_gradient2(low = '#98C1D9', mid = 'white', high = '#C81D25', midpoint = 0.5,
                               labels = seq(0, 1, 0.2), breaks = seq(0, 1, 0.2))+
          scale_x_discrete(labels = c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85))+
          labs(title = 'Subject Policy for Updating')+
          theme_tq()+
          theme(legend.key.width = unit(1, 'cm'))
        plotobject <- llout$evo %>%
          as.data.frame() %>%
          na.omit()
        plotobject2 <- plotobject  %>%
          pivot_longer(HIsim:SIsim, 'Attribution', values_to = 'values') %>%
          mutate(ret = ifelse(ret == 0.5, 1, 0))
        trial_wise <- ggplot( plotobject2, aes(trial, values, color = Attribution))+
            geom_line()+
            geom_point(aes(trial, ret), color = ifelse( plotobject2$ret == 1, 'black', 'dark red'))+
            scale_color_brewer(palette = 'Set1')+
            labs(x = 'Trial', y = 'Simulated Observation', title = 'Simulated Agent')+
            theme_tq()
        (trial_wise | (initial_policy/policy_change)) +
          plot_annotation(tag_levels = 'A') &
          theme(text = element_text(size = 20))
    })
}

# Run the application
shinyApp(ui = ui, server = server)
