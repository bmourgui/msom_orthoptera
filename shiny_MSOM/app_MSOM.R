

setwd(dir = here::here("shiny_MSOM"))
load(file=here::here("Data", "Data_analyse", "data_model_article.RData"))
load(here::here("Data", "Data_analyse", "out_articlev2.RData"))
out.model <- out_articlev2
points_pnm <- readr::read_csv(here::here("shiny_MSOM", "points_pnm.csv"))

library(shiny)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(RGraphics)
library(plyr)
library(reshape2)
library(tidyr)
library(leaflet)







#### Application ####

ui <- fluidPage(
    titlePanel("Supplementary materials"),
    sidebarPanel(
      selectInput(inputId = "nb1",
                  h3("Which species?"),
                  choices = list("Anonconotus occidentalis"=1, "Anonconotus ghilianii"=2, "Anonconotus mercantouri"=3, "Antaxius pedestris"=4,
                                 "Arcyptera fusca"=5, "Bicolorana bicolor"=6, "Calliptamus italicus"=7, "Calliptamus siciliae"=8,
                                 "Chorthippus apricarius"=9, "Chorthippus biguttulus"=10, "Chorthippus brunneus brunneus"=11, "Chorthippus dorsatus"=12,
                                 "Chorthippus saulcyi daimei"=13, "Chorthippus vagans vagans"=14, "Decticus verrucivorus verrucivorus"=15, "Tetrix depressa"=16, 
                                 "Ephippiger terrestris"=17, "Euchorthippus declivus"=18, "Euchorthippus elegantulus"=19, "Eupholidoptera chabrieri"=20, 
                                 "Euthystira brachyptera"=21, "Gomphocerus sibiricus sibiricus"=22,  "Gryllus campestris"=23, "Leptophyes punctatissima"=24, 
                                 "Metrioptera saussuriana"=25, "Myrmeleotettix maculatus"=26, "Nemobius sylvestris"=27, "Oecanthus pellucens"=28, 
                                 "Oedipoda caerulescens caerulescens"=29, "Oedipoda germanica"=30, "Omocestus haemorrhoidalis"=31, "Omocestus raymondi raymondi"=32, 
                                 "Omocestus rufipes"=33, "Omocestus viridulus"=34, "Pezotettix giornae"=35, "Pholidoptera aptera"=36, 
                                 "Pholidoptera fallax"=37, "Pholidoptera griseoaptera"=38, "Platycleis albopunctata"=39, "Podisma dechambrei"=40, 
                                 "Podisma pedestris"=41, "Polysarcus denticauda"=42, "Pseudochorthippus parallelus"=43, "Psophus stridulus"=44, 
                                 "Roeseliana roeselii"=45, "Sepiana sepium"=46, "Stauroderus scalaris"=47, "Stenobothrus cotticus"=48, 
                                 "Stenobothrus lineatus"=49, "Stenobothrus nigromaculatus"=50,  "Stenobothrus rubicundulus"=51, "Tessellana tessellata"=52, 
                                 "Tettigonia cantans"=53, "Tettigonia viridissima"=54, "Tylopsis lilifolia"=55, "Yersinella beybienkoi"=56),
                  selected = 1)),
    
    mainPanel(
      h2("Overview"),
      br(),
      p("The aim is to provide supplementary information to readers of the paper 'Multi-species occupancy models: an effective and flexible framework for studies of insect communities' by B. Mourguiart, T. Couturier, Y. Braud, J. Mansons, D. Combrisson and A. Besnard."),
      br(),
      
      h2("Detection map"),
      leafletOutput(outputId = "map"),
      uiOutput(outputId = "figS1"),
      br(),
      
      h2("Altitudinal distributions"),
      plotOutput(outputId = "occu"),
      uiOutput(outputId = "figS2"),
      br(),
      
      h2("Detection probability by technique"),
      plotOutput(outputId = "det_technique"),
      uiOutput(outputId = "figS3"),
      br(),
      
      h2("Grass height effect on sighting detection"),
      plotOutput(outputId = "det_grass"),
      uiOutput(outputId = "figS4"),
      br(),
      br()
      )
    )






server <- function(input, output) {
  
  
    
    
    i <- reactive({as.numeric(input$nb1)})
    
    spnames <- c("Anonconotus occidentalis", "Anonconotus ghilianii", "Anonconotus mercantouri", "Antaxius pedestris",
                "Arcyptera fusca", "Bicolorana bicolor", "Calliptamus italicus", "Calliptamus siciliae",
                "Chorthippus apricarius", "Chorthippus biguttulus", "Chorthippus brunneus brunneus", "Chorthippus dorsatus",
                "Chorthippus saulcyi daimei", "Chorthippus vagans vagans", "Decticus verrucivorus verrucivorus", "Tetrix depressa", 
                "Ephippiger terrestris", "Euchorthippus declivus", "Euchorthippus elegantulus", "Eupholidoptera chabrieri", 
                "Euthystira brachyptera", "Gomphocerus sibiricus sibiricus",  "Gryllus campestris", "Leptophyes punctatissima", 
                "Metrioptera saussuriana", "Myrmeleotettix maculatus", "Nemobius sylvestris", "Oecanthus pellucens", 
                "Oedipoda caerulescens caerulescens", "Oedipoda germanica", "Omocestus haemorrhoidalis", "Omocestus raymondi raymondi", 
                "Omocestus rufipes", "Omocestus viridulus", "Pezotettix giornae", "Pholidoptera aptera", 
                "Pholidoptera fallax", "Pholidoptera griseoaptera", "Platycleis albopunctata", "Podisma dechambrei", 
                "Podisma pedestris", "Polysarcus denticauda","Pseudochorthippus parallelus", "Psophus stridulus", 
                "Roeseliana roeselii", "Sepiana sepium", "Stauroderus scalaris", "Stenobothrus cotticus", 
                "Stenobothrus lineatus", "Stenobothrus nigromaculatus",  "Stenobothrus rubicundulus", "Tessellana tessellata", 
                "Tettigonia cantans", "Tettigonia viridissima", "Tylopsis lilifolia", "Yersinella beybienkoi")
    
    #output$species <- renderImage({
      #filename <- normalizePath(file.path('./www',
                                          #paste('sp', input$nb1, '.jpg', sep='')))
      
      # Return a list containing the filename and alt text
      #list(src = filename,
           #alt = dimnames(X)[[3]][[i()]])
      
    #}, deleteFile = FALSE)
    
    
    output$map <- renderLeaflet({
      
      data.frame("site"=cov.sites$station,
                 "alti"=cov.sites$altitude,
                 "det"=apply(X[,,i()], 1, function(x)sum(sum(x)>0))) %>% 
        mutate("psi"=plogis(out.model$mean$a0[i()] + out.model$mean$a1[i()]*cov.sites_standard$altitude + out.model$mean$a2[i()]*cov.sites_standard$altitude^2))%>%
        mutate("psi"=round(psi,2)) %>%
        inner_join(points_pnm, by=c('site'='id_releve')) -> points_pnm2
      
      points <- data.frame(longitudes = points_pnm2$X,
                           latitudes = points_pnm2$Y,
                           det = as.factor(points_pnm2$det),
                           occu=paste("mean(occupancy)=", points_pnm2$psi))
      
      pal <- colorFactor(  palette = c("red","darkgreen"),  domain = points$det)
      leaflet(points) %>%
        addTiles() %>%
        setView(lng = 7.14, lat = 44.103, zoom =10)%>%
        addCircles(lng = ~longitudes, lat = ~latitudes, 
                   weight = 10,
                   radius = 10, 
                   popup = paste0("Estimated occupancy probability: ",
                                  points_pnm2$psi,
                                  "<br>",
                                  "Site altitude: ",
                                  round(points_pnm2$alti,2),
                                  "</div>"),
                   color = ~pal(det), 
                   fillOpacity = 0.9) %>%
        addLegend("bottomleft", 
                  values = points$det,
                  pal = pal,
                  title = "Detection",
                  opacity = 1,
                  layerId = "legend")
      
    })
    
    output$figS1 <- renderUI({
      HTML(
        paste("Figure S1: Localization of the 81 mountain grasslands in the Mercantour National Park where Orthopera communities were sampled between July and October 2018. Green points correspond to sites where ", em(spnames[i()]), " was present and red points to sites where it was not detected. (Click on the points to have information on estimated occupancy probability and site altitude.)")
      )
    })
    
    
    output$occu <- renderPlot({
      
      alti.min <- floor(min(cov.sites$altitude)/100)*100
      alti.max <- ceiling(max(cov.sites$altitude)/100)*100
      o.ele <- seq(alti.min,alti.max,,500)
      mean.ele <- cov.sites_mean[1,1]
      sd.ele <- cov.sites_sd[1,1]
      ele.pred <- (o.ele - mean.ele)/sd.ele
      
      sim <- out.model$sims.list
      
      nsamp <- 3000
      predEsp <- array(0, dim=c(500, nsamp))
      for (m in 1: nsamp){
        predEsp[,m] <- plogis(sim$a0[m,i()] + sim$a1[m,i()]*ele.pred + sim$a2[m,i()]*ele.pred^2)
      }
      
      pmC <- apply(predEsp, 1, mean)
      criC <- apply(predEsp, 1, function(x)quantile(x, prob=c(0.025,0.975)))
      
      pred <- data.frame("o.ele"=o.ele,
                         "inf.occ"=criC[1,],
                         "occ"=pmC,
                         "sup.occ"=criC[2,])
      
      pres <- data.frame("alti"=cov.sites$altitude,
                         "data"=apply(X[,,i()], 1, function(x)sum(sum(x)>0)))
      
      pred %>%
        ggplot(aes(x=o.ele, y=occ)) +
        geom_line(size=1.5) +
        geom_ribbon(aes(ymin=inf.occ, ymax=sup.occ), alpha=0.2) +
        theme_classic() +
        ggtitle("") +
        ylim(c(0,1)) +
        geom_point(data = pres, mapping = aes(x=alti, y=data), color="black", size=2, shape=3) +
        ylab("Occupancy probability") +
        xlab("Elevation") +
        theme(axis.title = element_text(size=16),
              axis.text = element_text(size=14))
        
    
        
        
        
    })
    
    output$figS2 <- renderUI({
      HTML(paste(sep="", "Figure S2: Effect of the altitude on the occupancy probability for ", em(spnames[i()]),". The black line represents the posterior mean, and the grey surface corresponds to the 95% credible interval. "))
    })
    
    
    output$det_technique <- renderPlot({
      
      data_summary <- function(x) {
        m <- mean(x)
        ymin <- quantile(x=x, 0.025)
        ymax <- quantile(x=x, 0.975)
        return(c(y=m,ymin=ymin[[1]],ymax=ymax[[1]]))
      }
      
      sim <- out.model$sims.list
      
      nsamp <- 3000
      det <- data.frame("Technique"=rep(c("Sighting", "Listening", "Netting","Sighting", "Listening", "Netting"), each=nsamp), 
                        "Level"=rep(c("Plot-level","Site-level"), each=nsamp*3),
                        "Detection"=c(plogis(sim$b0[,i()]), plogis(sim$b0[,i()]+sim$b1[,i()]), plogis(sim$b0[,i()]+sim$b2[,i()]),
                                      1-(1-plogis(sim$b0[,i()]))^5, 1-(1-plogis(sim$b0[,i()]+sim$b1[,i()]))^5, 1-(1-plogis(sim$b0[,i()]+sim$b2[,i()]))^5))
      
      det %>%
        ggplot(aes(x=Technique, y=Detection)) +
        geom_violin() +
        facet_grid(.~Level) +
        stat_summary(fun.data=data_summary, color="black", position = position_dodge(0.8)) +
        theme_classic() +
        ggtitle("") +
        ylim(c(0,1)) +
        ylab("Detection probability") +
        xlab("") +
        theme(axis.title = element_text(size=16),
              axis.text = element_text(size=14),
              strip.text.x = element_text(size = 14))
      
      
      
      
      
    })
    
    
    output$figS3 <- renderUI({
      HTML(
        paste(sep="", "Figure S3: Posterior distribution of detection probabilities of ", em(spnames[i()])," for each sampling technique at plot- and site-level, i.e. after 5 plots. The black points are the means and the segments are the 95% credible intervals.")
      )
      
    })
    
    
    
    
    output$det_grass <- renderPlot({
      
      
      hauteur.min <- min(cov.detection$hauteur)
      hauteur.max <- max(cov.detection$hauteur)
      o.hauteur <- seq(hauteur.min, hauteur.max,,500)
      mean.hauteur <- cov.detection_mean[1,2]
      sd.hauteur <- cov.detection_sd[1,2]
      hauteur.pred <- (o.hauteur - mean.hauteur)/sd.hauteur
      
      sim <- out.model$sims.list
      
      nsamp <- 3000
      predEsp <- array(0, dim=c(500, nsamp))
      for (m in 1: nsamp){
        predEsp[,m] <- plogis(sim$b0[m,i()] + sim$b3[m,i()]*hauteur.pred)
      }
      
      pmC <- apply(predEsp, 1, mean)
      criC <- apply(predEsp, 1, function(x)quantile(x, prob=c(0.025,0.975)))
      
      pred <- data.frame("o.haut"=o.hauteur,
                         "inf.det"=criC[1,],
                         "det"=pmC,
                         "sup.det"=criC[2,])
      
      pred %>%
        ggplot(aes(x=o.haut, y=det)) +
        geom_line(size=1.5) +
        geom_ribbon(aes(ymin=inf.det, ymax=sup.det), alpha=0.2) +
        theme_classic() +
        ggtitle("") +
        ylim(c(0,1)) +
        #geom_point(data = pres, mapping = aes(x=alti, y=data), color="black", size=2, shape=3) +
        ylab("Detection probability") +
        xlab("Grass height") +
        theme(axis.title = element_text(size=16),
              axis.text = element_text(size=14))
      
      
      
      
      
    })
    
    output$figS4 <- renderUI({
      HTML(paste(sep="", "Figure S4: Effect of the grass height on the sighting detection probability of ", em(spnames[i()]),". The black line represents the posterior mean, and the grey surface corresponds to the 95% credible interval. "))
      })
    
    
    output$map <- renderLeaflet({
      
      data.frame("site"=cov.sites$station,
                 "alti"=cov.sites$altitude,
                 "det"=apply(X[,,i()], 1, function(x)sum(sum(x)>0))) %>% 
        mutate("psi"=plogis(out.model$mean$a0[i()] + out.model$mean$a1[i()]*cov.sites_standard$altitude + out.model$mean$a2[i()]*cov.sites_standard$altitude^2))%>%
        mutate("psi"=round(psi,2)) %>%
        inner_join(points_pnm, by=c('site'='id_releve')) -> points_pnm2
      
      points <- data.frame(longitudes = points_pnm2$X,
                           latitudes = points_pnm2$Y,
                           det = as.factor(points_pnm2$det),
                           occu=paste("mean(occupancy)=", points_pnm2$psi))
      
      pal <- colorFactor(  palette = c("red","darkgreen"),  domain = points$det)
      leaflet(points) %>%
        addTiles() %>%
        setView(lng = 7.14, lat = 44.103, zoom =10)%>%
        addCircles(lng = ~longitudes, lat = ~latitudes, 
                   weight = 10,
                   radius = 10, 
                   popup = paste0("Estimated occupancy probability: ",
                                   points_pnm2$psi,
                                   "<br>",
                                   "Site altitude: ",
                                   round(points_pnm2$alti,2),
                                   "</div>"),
                   color = ~pal(det), 
                   fillOpacity = 0.9) %>%
        addLegend("bottomleft", 
                  values = points$det,
                  pal = pal,
                  title = "Detection",
                  opacity = 1,
                  layerId = "legend")
      
    })
    
}



shinyApp(ui, server)

