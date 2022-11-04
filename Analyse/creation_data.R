#########################################################################%
#### Creation of datasets for the MSOM used in the published article ####%
#########################################################################%

# Last modification: 24 oct 2022
library(magrittr)

data <- readxl::read_excel(here::here("Data/Data_analyse/data_PNM_ortho_2018.xlsx"))
data_station <- readxl::read_excel(here::here("Data/Data_analyse/data_PNM_ortho_2018.xlsx"), sheet = 2)

data %>%
  dplyr::filter(taxon != "Gomphocerinae sp") %>%
  dplyr::filter(taxon != "Calliptamus sp") %>%
  dplyr::filter(taxon != "Oedipoda sp") %>%
  dplyr::filter(taxon != "Mantis religiosa") %>%
  dplyr::filter(taxon != "Pijnackeria masettii") %>%
  dplyr::mutate("taxon"=as.factor(taxon)) -> data
levels(data$taxon) <- c("Anonconotus occidentalis", #rename species according to actual nomenclature
                        "Anonconotus ghiliani",
                        "Anonconotus mercantouri",
                        "Antaxius pedestris",
                        "Arcyptera fusca",
                        "Bicolorana bicolor",
                        "Calliptamus italicus",
                        "Calliptamus siciliae",
                        "Gomphocerippus apricarius",
                        "Gomphocerippus biguttulus",
                        "Gomphocerippus brunneus brunneus",
                        "Chorthippus dorsatus",
                        "Gomphocerippus saulcyi daimei",
                        "Gomphocerippus vagans vagans",
                        "Decticus verrucivorus verrucivorus",
                        "Depressotetrix depressa",
                        "Ephippiger terrestris",
                        "Euchorthippus declivus",
                        "Euchorthippus elegantulus",
                        "Eupholidoptera chabrieri",
                        "Euthystira brachyptera",
                        "Gomphocerus sibiricus sibiricus",
                        "Gryllus campestris",
                        "Leptophyes punctatissima",
                        "Metrioptera saussuriana",
                        "Myrmeleotettix maculatus",
                        "Nemobius sylvestris",
                        "Oecanthus pellucens",
                        "Oedipoda caerulescens caerulescens",
                        "Oedipoda germanica",
                        "Omocestus haemorrhoidalis",
                        "Omocestus raymondi raymondi",
                        "Omocestus rufipes",
                        "Omocestus viridulus",
                        "Pezotettix giornae",
                        "Pholidoptera aptera",
                        "Pholidoptera fallax",
                        "Pholidoptera griseoaptera",
                        "Platycleis albopunctata",
                        "Podisma dechambrei",
                        "Podisma pedestris",
                        "Polysarcus denticauda",
                        "Pseudochorthippus parallelus",
                        "Psophus stridulus",
                        "Roeseliana roeselii",
                        "Sepiana sepium",
                        "Stauroderus scalaris",
                        "Stenobothrus cotticus",
                        "Stenobothrus lineatus",
                        "Stenobothrus nigromaculatus",
                        "Stenobothrus rubicundulus",
                        "Tessellana tessellata",
                        "Tettigonia cantans",
                        "Tettigonia viridissima",
                        "Tylopsis lilifolia",
                        "Yersinella beybienkoi")


data  %>%
  dplyr::select(c(station, releve, taxon)) %>%
  dplyr::mutate(pres=1) %>%
  tidyr::spread(taxon, pres) -> tab_contingence

tab_contingence[is.na(tab_contingence)] <- 0


tab_contingence %>%
  dplyr::select(-2) -> t

t$rep <- 1
t <- t[,c(1,ncol(t),2:(ncol(t)-1))]
for (j in 2:nrow(t)){
  if(t[j-1,1]==t[j,1])(t[j,"rep"] <- t[j-1,"rep"]+1)
}

t %>%
  tidyr::gather("species","Occ", -c(1:2)) %>%
  reshape2::melt(id.var=c("species", "station", "rep"), measure.var="Occ") %>%
  reshape2::acast(station ~ rep ~ species) -> X 

X[is.na(X)] <- 0

data  %>%
  dplyr::select(c(station, releve, taxon , entendu)) %>%
  tidyr::spread(taxon, entendu) -> contingence_entendu

contingence_entendu[is.na(contingence_entendu)] <- 0


data  %>%
  dplyr::select(c(station, releve, taxon , vu)) %>%
  tidyr::spread(taxon, vu) -> contingence_vu

contingence_vu[is.na(contingence_vu)] <- 0


data  %>%
  dplyr::select(c(station, releve, taxon , fauche)) %>%
  tidyr::spread(taxon, fauche) -> contingence_fauche

contingence_fauche[is.na(contingence_fauche)] <- 0

contingence_methode <- rbind(contingence_vu, contingence_entendu, contingence_fauche)

contingence_methode %>%
  dplyr::select(-2) -> t2

t2$met <- rep(c("aVU","E","F"), each=403)
t2$rep <- 1
t2 <- t2[,c(1,(ncol(t2)-1):ncol(t2),2:(ncol(t2)-2))]
for (j in 2:nrow(t2)){
  if(t2[j-1,1]==t2[j,1])(t2[j,"rep"] <- t2[j-1,"rep"]+1)
}

t2 %>%
  dplyr::mutate(rep=paste(met, rep, sep="-")) %>%
  dplyr::select(-2) %>%
  tidyr::gather("species","Occ", -c(1:2)) %>%
  reshape2::melt(id.var=c("species", "station", "rep"), measure.var="Occ") %>%
  reshape2::acast(station ~ rep ~ species) -> X2 #matrice detection/non-detection 3D [j,k*m,i] : j = site, k = replicat, m = methode, i = especes

X2[is.na(X2)] <- 0

J <- dim(X2)[1] #number of sites
n <- dim(X2)[3] #species number
K <- dim(X2)[2] #visits number (5 visits per method)


#Data augmenation
nzeroes <- 50 #number of species non detected added
X.zero2 = matrix(0, nrow=J, ncol=K)
X.zero2[is.na(X2[,,1])] <- NA

#Xaug is the augmented version of X.  The first n species were actually observed
#and the n+1 through nzeroes species are all zero encounter histories  
X <- array(0, dim=c(dim(X2)[1],dim(X2)[2],dim(X2)[3]+nzeroes))
X[,,(dim(X2)[3]+1):dim(X)[3]] <- rep(X.zero2, nzeroes)
X[,,1:dim(X2)[3]] <-  X2
dimnames(X) <- list("station"=dimnames(X2)[[1]], "rep"=dimnames(X2)[[2]], "species"=c(dimnames(X2)[[3]],paste("sp",1:nzeroes)))


#Covariables de detection
methodeV <- matrix(0, nrow=dim(X2)[1], ncol=dim(X2)[2])
methodeV[,1:5] <- 1
methodeE <- matrix(0, nrow=dim(X2)[1], ncol=dim(X2)[2])
methodeE[,6:10] <- 1
methodeF <- matrix(0, nrow=dim(X2)[1], ncol=dim(X2)[2])
methodeF[,11:15] <- 1

tab_contingence %>%
  dplyr::inner_join(data, by=c("station","releve")) %>%
  dplyr::filter(!duplicated(releve)) %>%
  dplyr::select(c(station, releve, heure,hauteur)) %>%
  dplyr::mutate("heure"=lubridate::hour(heure)+ lubridate::minute(heure)/60) %>%
  # manually add two replicates where no orthoptera were detected and no noted in the dataframe 'data'
  dplyr::add_row(station = "108",
                 releve = 896,
                 heure = 16 + 13 / 60,
                 hauteur = 3) %>%
  dplyr::add_row(station = "336",
                 releve = 763,
                 heure = 11 + 51 / 60,
                 hauteur = 4) %>% 
  dplyr::arrange(station) -> cov.detection

cov.detection_mean <- matrix(data = colMeans(cov.detection[,-c(1:2)]), nrow=nrow(cov.detection), ncol=ncol(cov.detection)-2, byrow = T)
cov.detection_sd <- matrix(data = apply(cov.detection[,-c(1:2)],2, FUN=sd), nrow=nrow(cov.detection), ncol=ncol(cov.detection)-2, byrow = T)
cov.detection_standard <- cbind(cov.detection[,c(1:2)], ((cov.detection[,-c(1:2)]-cov.detection_mean)/cov.detection_sd))


cov.detection_standard$releve <- 1
cov.detection_standard <- cov.detection_standard[,c(1,ncol(cov.detection_standard),2:(ncol(cov.detection_standard)-1))]
for (j in 2:nrow(cov.detection_standard)){
  if(cov.detection_standard[j-1,1]==cov.detection_standard[j,1])(cov.detection_standard[j,"releve"] <- cov.detection_standard[j-1,"releve"]+1)
}

cov.detection_standard %>%
  tidyr::gather("Cov","Value", -c(1,3)) %>%
  reshape2::melt(id.var=c("Cov", "station", "releve"), measure.var="Value") %>%
  reshape2::acast(station ~ releve ~ Cov) %>%
  apply(c(2,3), FUN=as.numeric) -> cov.detection_standard


cov.detection_hauteur <- matrix(rep(cov.detection_standard[,,"hauteur"], 3), ncol=15)
cov.detection_heure <- matrix(rep(cov.detection_standard[,,"heure"], 3), ncol=15)



#Covariables de sites potentielles
data_station %>%
  dplyr::inner_join(data, by="station") %>%
  dplyr::select(station, exposition, "altitude"=altitude.y, hauteur) %>%
  dplyr::group_by(station) %>%
  dplyr::mutate_at(-1, as.numeric) %>%
  dplyr::summarise_if(is.numeric, .funs=mean) %>%
  dplyr::right_join(data.frame(station=rownames(X2)), by="station") -> cov.sites



cov.sites_mean <- matrix(data = colMeans(cov.sites[,-c(1:2)]), nrow=nrow(cov.sites), 
                         ncol=ncol(cov.sites)-2, byrow = T)
cov.sites_sd <- matrix(data = apply(cov.sites[,-c(1:2)],2, FUN=sd), nrow=nrow(cov.sites), 
                       ncol=ncol(cov.sites)-2, byrow = T)
cov.sites_standard <- (cov.sites[,-c(1:2)]-cov.sites_mean)/cov.sites_sd



cov.sites$exposition %>%
  Hmisc::cut2(cuts = seq(0, 360, 45)) -> exposition
levels(exposition) <- c(0,1,1,1,1,0,0,0)
exposition <- as.numeric(exposition)-1



cov.sites_standard <- cbind("station"=cov.sites$station, exposition, cov.sites_standard)



save(tab_contingence, cov.detection, cov.detection_mean, cov.detection_sd, 
     cov.detection_hauteur, cov.detection_heure, 
     methodeV, methodeE, methodeF, cov.sites, cov.sites_standard, cov.sites_mean, cov.sites_sd,
     X, n,J,K,nzeroes, file=here::here("/Data/Data_analyse/data_model_article.RData"))
