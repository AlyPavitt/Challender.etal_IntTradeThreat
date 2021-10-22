
rm(list = ls())

#####################################################################################################################################################
#####################################################################################################################################################

# Packages ----------

library("here")
library("reshape2")


#####################################################################################################################################################
#####################################################################################################################################################

# Full automated assessment from IUCN raw data
# STEP 1: database query - code up the initial IUCN database query (ies?) - actually might be able to get away without doing this, since the selection should exclude anything that doesnt fit
# STEP 2: first step automated assessment 
# STEP 3: second step automated assessment following decision tree process

#####################################################################################################################################################
#####################################################################################################################################################


## prepare data ##
# code works on raw downloads of the following files from https://www.iucnredlist.org/ : "assessment", "taxonomy", "threats" and "usetrade"


#####################################################################################################################################################


# read in RL data files -----------
# note that if multiple downloads are taken from the IUCN Red List website (e.g. for different taxonomic subsets), you will need to rbind() them before running the script below
# e.g. to rbind multiple csv files saved in the folder "Assessments": 
    #IUCNassessment <- do.call(rbind,lapply(paste0(paste0(here(),"/Input/Assessments"),sep="/",dir(paste0(here(),"/Input/Assessments"))),read.csv))

# assessments (core data - species, free text for threat, rationale, conservation, useTrade)
IUCNassessment <- read.csv(paste0(here(),"/Input/Assessment.csv"))

# taxonomy
IUCNTaxonomy <- read.csv(paste0(here(),"/Input/Taxonomy.csv"))

# threats (species threats)
IUCNThreats <- read.csv(paste0(here(),"/Input/Threats.csv"))

# usetrade (species end uses)
IUCNUseTrade <- read.csv(paste0(here(),"/Input/useTrade.csv"))


#####################################################################################################################################################


# Universal vectors for threat codes
# full descriptions of threat classifications can be found here: https://www.iucnredlist.org/resources/threat-classification-scheme 

threatcodes <- c("5.1.1","5.1.4","5.4.1","5.4.2","5.4.4","5.4.6","5.2.1","5.2.4","5.3.1","5.3.2","5.3.5")
threatcodes.animals <- c("5.1.1","5.1.4","5.4.1","5.4.2","5.4.4","5.4.6")
threatcodes.plants <- c("5.2.1","5.2.4","5.3.1","5.3.2","5.3.5")


#####################################################################################################################################################


# add higher taxonomy ----------
traderaw <- merge(IUCNTaxonomy[c("internalTaxonId","kingdomName","phylumName","className","orderName","familyName")],
                  IUCNassessment[c("assessmentId","internalTaxonId","redlistCategory","scientificName","rationale","threats","useTrade")], by="internalTaxonId")


# add threat codes  ----------

# IS species where threat$internationalTrade == "YES"
# exclude instances where the threat timing is "Past, Unlikely to Return"
# create separate vector because multiple threats might have different responses
Species_internationalTrade <- c(IUCNThreats$internalTaxonId[IUCNThreats$internationalTrade == "Yes" & !IUCNThreats$timing == "Past, Unlikely to Return"])

IUCNThreats <- IUCNThreats[order(IUCNThreats$code),]
species.threats <- reshape2::dcast(IUCNThreats[which(IUCNThreats$code %in% c(threatcodes) & !IUCNThreats$timing=="Past, Unlikely to Return"),],internalTaxonId~"ThreatCodes",
                                   value.var="code",fun.aggregate=function(x) paste(x,collapse="; ")) # existing threats

species.threats_Past <- reshape2::dcast(IUCNThreats[which(IUCNThreats$code %in% c(threatcodes) & IUCNThreats$timing=="Past, Unlikely to Return"),],internalTaxonId~"ThreatCodes_PastUnlikely",
                                        value.var="code",fun.aggregate=function(x) paste(x,collapse="; "))

species.threats$internationalTrade [species.threats$internalTaxonId %in% Species_internationalTrade] <- "Yes"

traderaw <- merge(traderaw,merge(species.threats, species.threats_Past,by="internalTaxonId",all=T),by="internalTaxonId",all.x=T) # add to trade


# add end uses ----------
# details of the use and trade classifications can be found here: https://www.iucnredlist.org/resources/general-use-trade-classification-scheme 
# each use and trade classification can also be classified at the scale of use: international, national and/or subsistence
# this is not a compulsory field in the Red List assessment and so may not always be completed

# number of end uses under "intenational", "national" and "subsistence"
enduse <- data.frame(merge(table(IUCNUseTrade$internalTaxonId[IUCNUseTrade$international=="true"]),
                           merge(table(IUCNUseTrade$internalTaxonId[IUCNUseTrade$national=="true"]),
                                 table(IUCNUseTrade$internalTaxonId[IUCNUseTrade$subsistence=="true"]), by="Var1", all=T), by="Var1", all=T))
names(enduse) <- c("internalTaxonId","International","National","Subsistence")


# concatenate end uses for "international"
enduse.int <- reshape2::dcast(IUCNUseTrade[which(IUCNUseTrade$international=="true"),],internalTaxonId~"International_EndUses",
                              value.var="name",fun.aggregate=function(x) paste(x,collapse="; "))

enduse.intONLY <- reshape2::dcast(IUCNUseTrade[which(IUCNUseTrade$international=="true" & IUCNUseTrade$subsistence=="" & IUCNUseTrade$national==""),],
                                  internalTaxonId~"InternationalOnly_EndUses",
                                  value.var="name",fun.aggregate=function(x) paste(x,collapse="; "))

# combine number of end uses and lists of international end uses with previous df
multi_merge_function <- function(df1, df2){
  merge(df1, df2, by="internalTaxonId", all=T) 
} # "full outer join", by internalTaxonId

traderaw <- Reduce(multi_merge_function, list(traderaw,enduse,enduse.int,enduse.intONLY)) 
traderaw <- unique(traderaw)


#####################################################################################################################################################


# clean free text -----------

# key column names
cn <- c("threats","rationale","useTrade")
usetrade <- c("Subsistence","National","International")

# standardise case and remove unwanted characters
traderaw[,c(cn)] <- sapply(traderaw[,c(cn)],function(x) tolower(x)) # to lower
traderaw[,c(cn)] <- sapply(traderaw[,c(cn)],function(x) gsub("<[^>]+>","",x)) # exclude everything between < and > 
traderaw[,c(cn)] <- sapply(traderaw[,c(cn)],function(x) gsub("[^a-zA-Z|^ ]", "", x)) # remove all characters except alpha and space
traderaw[,c(cn)] <- sapply(traderaw[,c(cn)],function(x) gsub("   |  ", " ", x)) # remove double and triple spaces

# standardise abbreviations
traderaw[,c(cn)] <- sapply(traderaw[,c(cn)],function(x) gsub("cannot|can not", "cant", x)) 
traderaw[,c(cn)] <- sapply(traderaw[,c(cn)],function(x) gsub("is not", "isnt", x)) 
traderaw[,c(cn)] <- sapply(traderaw[,c(cn)],function(x) gsub("will not", "wont", x))
traderaw[,c(cn)] <- sapply(traderaw[,c(cn)],function(x) gsub("does not", "doesnt", x))
traderaw[,c(cn)] <- sapply(traderaw[,c(cn)],function(x) gsub("has not", "hasnt", x))
traderaw[,c(cn)] <- sapply(traderaw[,c(cn)],function(x) gsub("have not", "havent", x))

# standardise other keywords
traderaw[,c(cn)] <- sapply(traderaw[,c(cn)],function(x) gsub("utiliz", "utilis", x)) 

# remove stopwords
stop <- c("a","about","also","an","and","any","are","as","be","been","being","but","by","can","for","from","has","have","if","in","includ.?.?.?","into","is","it","its","much","of","on","or","so","still","such","than","that","the","there","therefore","these","this","to","were","where","which","with","within","would", # standard stops
          "species","wild") # specific stops
stop <- paste0(" ",c(stop)," ")


# if traderaw[,c(cn)] contains stopwords, replace with space and repeat until no more stopwords
repeat{
  trade.stop <- unique(traderaw[sapply(traderaw[c(cn)], grepl, pattern = paste(stop,collapse="|")),])
  
  if(nrow(trade.stop)=="0"){
    break
  }
  
  if(nrow(trade.stop)>0){
    traderaw[,c(cn)] <- sapply(traderaw[,c(cn)],function(x) gsub(paste(stop,collapse="|"), " ", x))
  }
}


#####################################################################################################################################################
#####################################################################################################################################################


## First queries to identify shortlist of species ##

# Defining parameters ----------

# #1. species categorised as Critically Endangered (CR), Endangered (EN), Vulnerable (VU), Near Threatened (NT), 
# Low Risk/near threatened (LR/nt) or Low Risk/conservation dependent (LR/cd); and either 
# (a) >= one of 55 text strings in rationale, threats, and/or use and trade; or 
# (b) >= one of 11  biological resource use threat codes (5.1.1, 5.1.4, 5.2.1, 5.2.4, 5.3.1, 5.3.2, 5.3.5, 5.4.1, 5.4.2, 5.4.4, 5.4.6)
# #2. "international" use under end uses = TRUE (identify any additional species to Q1)
# #3. "International trade is a significant driver of threat" = TRUE (identify any additional species to Q1 and Q2)


#####################################################################################################################################################


# vectors ----------
threatcategory <- c("Critically Endangered", "Endangered", "Vulnerable", "Near Threatened", "Lower Risk/near threatened", "Lower Risk/conservation dependent")
threatstrings.major <- c("trade", "export", "collect", "enthusiast", "harvest", "exploit", "demand", "market", "cites")
threatstrings.minor1 <- c("nternational", "ransboundary", "ntercontinental", "border", "pet", "aquari", "horticultur", "timber", "commercial", "world", "ornamental",
                          "cage.?bird", "curio", "medicin", "fisher", "regional", "pharmaceutical", "unsustainabl")
threatstrings.minor2 <- c(" use", "utili")


# 1. text string | BRU threat code ----------
# text string major = keyword alone
traderaw$Shortlisted <- NA
traderaw$Shortlisted [grepl(paste(threatcodes,collapse="|"),traderaw$ThreatCodes) | grepl(paste(threatcodes,collapse="|"),traderaw$ThreatCodes_PastUnlikely) | 
                        !!rowSums(sapply(traderaw[c("threats", "rationale", "useTrade")], grepl, pattern = paste(threatstrings.major,collapse="|")))] <- "yes.1a"

# text string minor = must include at least one from each of minor1 and minor2
traderaw$Shortlisted [grepl(paste(threatstrings.minor1,collapse="|"),traderaw$threats) & grepl(paste(threatstrings.minor2,collapse="|"),traderaw$threats)] <- "yes.1b"
traderaw$Shortlisted [grepl(paste(threatstrings.minor1,collapse="|"),traderaw$rationale) & grepl(paste(threatstrings.minor2,collapse="|"),traderaw$rationale)] <- "yes.1b"
traderaw$Shortlisted [grepl(paste(threatstrings.minor1,collapse="|"),traderaw$useTrade) & grepl(paste(threatstrings.minor2,collapse="|"),traderaw$useTrade)] <- "yes.1b"


# 2. international end use ----------
traderaw$Shortlisted [is.na(traderaw$Shortlisted) & traderaw$International>0] <- "yes.2"


# 3. international trade is a signficiant driver of threat ----------
# include all timings, including "past, unlikely to return" or not
traderaw$Shortlisted [is.na(traderaw$Shortlisted) & traderaw$internalTaxonId %in% c(IUCNThreats$internalTaxonId[IUCNThreats$internationalTrade == "Yes"])] <- "yes.3"


# anything not in accepted red list categories = NA ----------
traderaw$Shortlisted [!traderaw$redlistCategory %in% c(threatcategory)] <- NA

# take shortlisted taxa forwards onto first automation step

trade <- traderaw[grepl("yes", traderaw$Shortlisted),]


#####################################################################################################################################################
#####################################################################################################################################################


## General assessment [ first automation step ] ##

# Defining parameters of scoring ----------

# Three basis for automated scoring step 1
# #1: "Use and Trade" field contains keywords indicating that use and trade are (a) unlikely to be a threat (score "U") or (b) there is insufficient information (score "I") 
# #2: "International trade is a significant threat" field is selected (score "L")
# #3: An animal was selected in the initial IUCN RL query on the basis of plant-related threat keywords or threat scores only (score "U")


#####################################################################################################################################################


# 1. use and trade keywords ----------

# trade-related keywords -----------
TradeKeywords1 <- c("trad.?.?.?","harvest.?.?.?","export.?.?.?.?.?.?","exploit.?.?.?.?.?.?","utilis.?.?.?.?.?.?","trapp.?.?.?.?","collect.?.?.?.?")

################################################################

###### CODE U ######

# keywords indicating that there is no trade/ international trade or it is not a threat (U) -----------

NoTrade_Prefix1 <- c("no","not","not known","not found","no.? current.?.?","current.?.? no.?","no evidence","no evidence current","no.? direct.?.?","no existing",
                     "neither","appears not","no information suggests","not subject", "not thought use","not thought","no.? document.?.?","no.? apparent",
                     "no known specific","no.? known current.?.?","no.? know.?.?.?.?.?.?","no documentation","isnt","no known",
                     "no.? report.?.?","no.? record.?.?","no record.?.?.?.?","no report.?.?.?.?","no.? re.?.?.?.?ed.?.?.?","no report.?.?.?.? know.?","never re.?.?.?.?ed",
                     "not believed","no.? thought","doesnt appear","doesnt appear use.?","unlikely","no.? significant.?.? use.?")

NoTrade_Keywords1 <- c(paste(sapply(NoTrade_Prefix1,paste,c(TradeKeywords1,"international trad","use"))),"not targeted specifically","none known",
                       "no known human use","no use.? know","no.? know.? use","no report.?.?.?.? use","no know.? record.?.?.?.? use","no know.? report.?.?.?.? use",
                       "no known removal trade","not found commercial pet trade","no way use.? trad...?","trade isnt apparent",
                       "trade unrecorded") 

###### CODE I ######

# keywords indicating that  there is insufficient information (I) to know whether trade is a threat -----------

UnknownTrade_Keywords1 <- c("no use.? information","no trad.?.?.? information","no use.?trade information","no.? known whether use",
                            "no information use.?","no information trade","no information regarding use.?","no information regarding trade","no information regarding.?.?.?.? trade",
                            "no information found use.?","no information found trade","no information known about use","no information suggest.? trade",
                            "no specific use.?trade information","no specific trade information","no information.?.?.?.? available","no information about .?.?.?.?use",
                            "no current information use.?","no current information trade","no data particular use.?","no data particular trad.?",
                            "no known information use","no known information trade","no known information about use","no known information about trade",
                            "more information needed","more information required","no available data","no information whether harvest",
                            "information use.?trade not available","information isnt available","no specific data available","no.? available use.?trade",
                            "use.?trad.?.?.? information not available","trade status unknown","trade.?use unknown","use.?trade unknown",
                            "unknown use.?trade","unknown whether use","unknown whether trade","use.? unknown","trad.?.?.? unknown","unknown trade","information.?.? unknown",
                            "trade information unknown","use.? specifically unknown","trad.?.?.? specifically unknown",
                            "information regarding trade.?use isnt known","information.?.? missing",
                            "no report.?.? information","no report.?.? use.? information","no record.?.? information use","no.? re.?.?.?.?.?.? information trade",
                            "precaution warranted")

###########################################################

# score based on keywords in use and trade AND international end uses that ----------
# U = indicate it is not in trade, or
# I = indicate that it is unknown whether it is in trade due to insufficient information
trade$Score <- ""
trade$Score[which(grepl(paste(c(NoTrade_Keywords1),collapse="|"),trade$useTrade)& is.na(trade$International))] <- "U"
trade$Score[which(grepl(paste(c(UnknownTrade_Keywords1),collapse="|"),trade$useTrade))] <- "I"

# indicate basis of the score
trade$BasisOfScore <- ""
trade$BasisOfScore[which(grepl(paste(c(NoTrade_Keywords1),collapse="|"),trade$useTrade) & is.na(trade$International))] <- "TradeAndUse_keywords"
trade$BasisOfScore[which(grepl(paste(c(UnknownTrade_Keywords1),collapse="|"),trade$useTrade))] <- "TradeAndUse_keywords"


#####################################################################################################################################################


# 2. "International trade is a significant threat" field is selected (score "L")
# excluding timing = "Past, unlikely to return"
# "International trade is a significant threat" is not a compulsory field in the Red List assessment and so may not always be completed

trade$Score[trade$internationalTrade=="Yes"] <- "L"
trade$BasisOfScore[trade$internationalTrade=="Yes"] <- "IntTradeThreat"


#####################################################################################################################################################


# 3. Animal selected under plant threats ----------
# An animal was selected in the initial IUCN RL query on the basis of plant-related threat keywords or threat scores only (score "U")

# identify instances where animal species 
# (a) DO include plant-related keywords and do NOT include animal-related keywords AND
# (b) DO include plant-related threat scores and do NOT include animal-related threat scores AND
# (c) "International trade is a significant threat" is NOT selected

# (a) identify where plant-related keywords are included ----------
trade$Score4 [!!rowSums(sapply(trade[c("threats", "rationale", "useTrade")], grepl, pattern = paste(c("horticultur","timber"),collapse="|"))) & 
                !!rowSums(sapply(trade[c("threats", "rationale", "useTrade")], grepl, pattern = paste(threatstrings.minor2,collapse="|")))] <- "plantonly"

# ignore any instances where animal and non-taxon-specific keywords were also included
trade$Score4 [!!rowSums(sapply(trade[c("threats", "rationale", "useTrade")], grepl, pattern = paste(c(threatstrings.major,"pet","cage.?bird","fisher"),collapse="|")))] <- NA


# (b) ignore where animal-related threat scores OR (c) "International trade is a significant threat" are used ----------
trade$Score4 [grepl(paste(threatcodes.animals, collapse="|"), trade$ThreatCodes) | grepl(paste(threatcodes.animals, collapse="|"), trade$ThreatCodes_PastUnlikely) | 
                trade$internationalTrade =="yes"] <- NA


# identify species that meet either the initial query keywords OR the initial query threats ----------
# with additional caveat of "timber" as a keyword and five threat codes relating to plants/timber only refering to plants
trade$Score[trade$kingdomName %in% c("ANIMALIA","CHROMISTA","FUNGI") & trade$Score4=="plantonly"] <- "U"

trade$BasisOfScore[trade$kingdomName %in% c("ANIMALIA","CHROMISTA","FUNGI") & trade$Score4=="plantonly"] <- "Animal_PlantThreat"

trade$Shortlisted <- trade$AnimalThreat4 <- trade$Score4 <- NULL #clean unnecessary datafields


#####################################################################################################################################################
#####################################################################################################################################################


## General assessment [ second automation step ] ##

# Following the rules of the decision tree -----------
# note that the first step of the decision tree "is international trade a sig driver of trade" = YES 
# already forms part of the first automated scoring

# 1. Is there evidence that use and/or trade takes place?
# (a) No/past/future/potential/possiblr
# (b) Yes/probable 
# if yes/probably keywords = 1(b) ELSE 1(a)

# 2a. IF 1(a) THEN is there evidence that use and/or trade is a potential future threat?
# (a) No == U1
# (b) Yes == I3

# 2b. IF 1(b) THEN is there evidence that use and/or trade is NOT international
# (a) No
# (b) Yes == U2
# if 2b(b) ELSE 2b(a)
# 2b(b) = international is not ticked but national/subsistance is,  use/trade NOT international keywords

# 3. IF 2b(a) THEN is there evidence that use/trade is NOT a threat?
# (a) No 
# (b) Yes == U3
# if 3(b) ELSE 3(a)
# 3(b) = threat codes OTHER THAN use threat codes ticked, use/trade NO threat keywords

# 4. IF 3(a) THEN is there evidence that use/trade IS a threat?
# (a) No/past/future/potential/possible = I2
# (b) Yes/probable
# if 4(b) ELSE 4(a)
# 4(b) = use-related threat code selected, use/trade as threat keywords

# 5. IF 4(b) THEN is there evidence that use/trade is international?
# (a) yes/probably = L2
# (b) No/past/future/potential/possible = I1


#####################################################################################################################################################


# for the purposes of proof of concept, focus on amphibia and actinopterygii
trade <- trade[trade$className %in% c("AMPHIBIA","ACTINOPTERYGII"),]
#trade <- trade[trade$kingdomName == "ANIMALIA",]

#####################################################################################################################################################
#####################################################################################################################################################


## Set keyword vectors for all steps of the decision tree ##


#####################################################################################################################################################


### DECISION TREE 1: IS THERE EVIDENCE THAT USE AND/OR TRADE TAKES PLACE? ###

# TRADE 
Trade_core <- c("trad.?.?.?","harvest.?.?.?","export.?.?.?.?.?.?","utilis.?.?.?.?.?.?","trapp.?.?.?.?","collect.?.?.?.?") #"exploit.?.?.?.?.?.?"

Trade_prefix <- c("medicin.?.?","food","cosmetic.?.?.?.?","curio.?","ornamental.?.?","animal.?","pet.?","zoo.?",
                  "fish.?.?.?.?.?","aquari.?.?.?",
                  "bird.?.?.?","falcon.?.?.?.?.?") 

# trade suffix needs some more strings to account for addition text after trade_core
Trade_suffix <- c("medicin.?.?","food","cosmetic.?.?.?.?","curio.?","ornamental.?.?","animal.?","liv.?.?.? animal.?","pet.?","zoo.?",
                  "fish.?.?.?.?.?","line.?fish.?.?.?.?.?","sport.? fish.?.?.?.?.?","aquari.?.?.?","aquari.?.?.? pet.?",
                  "bird.?.?.?","cage.?bird.?","falcon.?.?.?.?.?") 

#Trade_prefixPLANT <- c("plant.?","horticulur.?.?","liv.?.?.? plant.?","timber") # need to incorporate plant keywords later

TradeKeywords <- c(paste(sapply(Trade_prefix,paste,Trade_core)),
                   paste(sapply(Trade_core,paste,Trade_suffix)),
                   paste(sapply(Trade_core,paste,c("trade"))))

##########################

# UNKNOWN TRADE

UnknownTrade_Keywords <- c("no information utilisation", "no information use trade", "no information potential use trade", "no report.?.? stud.?.?.? use trade")

##########################

# NO TRADE
# including unlikely and past trade

NoTrade_Prefix <- c("not","not known","not found","not found commerical","not subject","neither","no interest","unlikely","no known removal",
                    "no existing","no.? current.?.?","current.?.? no.?","no evidence","no evidence .?.?.?.?.?.?.?","not believed",
                    "no.? direct.?.?","appears not","no information suggest.?.?.?","no documentation","doesnt appear",
                    "not thought use.?","not thought","no.? document.?.?","no.? apparent.?.?",
                    "no.? know.?.?.?.?","no.? known current.?.?","no known specific","no.? significant.?.?",
                    "no re.?or.?.?.?","not re.?or.?.?.?","no known re.?or.?.?.?","no re.?or.?.?.? known","never re.?or.?.?.?",
                    "historical") 

NoTrade_Keywords <- c(paste(sapply(NoTrade_Prefix,paste,c(Trade_prefix,Trade_core,"international trad","use"))),
                      "not targeted specifically","no known human use","no use.? know","no.? know.? use",
                      "no re.?or.?.?.?.?.? use","no know.? re.?or.?.?.?.?.? use","trade unre.?or.?ed",
                      "no way use.??trad...?","trade isnt apparent","trade not issue",
                      "past some commercial") 


################################################################


### DECISION TREE 2a: IS THERE EVIDENCE THAT USE/TRADE IS A POTENTIAL FUTURE THREAT? ###

# is genus in trade?
# possible the tree is utilised
# likely attract the attention of  specialist collectors if recollected

Future_Suffix <- c("may","future","potential","possib.?l.?.?.?","rarely issue")
FutureKeywords <- c(paste(sapply(c(Trade_core,Trade_prefix),paste,Future_Suffix)),
                    paste(sapply(Future_Suffix,paste,c(Trade_core,Trade_prefix))),
                    "collector.? item")


################################################################


### DECISION TREE 2b & 5: IS THERE EVIDENCE THAT USE/TRADE (2b) IS NOT or (5) IS INTERNATIONAL? ###
#

# INTERNATIONAL
International_prefix <- c("international.?.?.?","trans.?boundar.?.?.?","global.?.?","world","export.?.?.?.?.?.?")

Int_TradeKeywords <- unique(c(paste(sapply(International_prefix,paste,Trade_core)),
                              paste(sapply(International_prefix,paste,Trade_prefix)),
                              paste(sapply(International_prefix,paste,c("market","demand","touris.?.?"))),
                              paste(sapply(Trade_core,paste,International_prefix)),
                              #"ebay","internet trad.?.?.?","internet collect.?.?.?.?","online demand.?","online market.?","webmarket.?","on.?.?.?.? internet",
                              "commercial.?.? overexploit.?.?.?.?.?","international.?.? sought after","smuggl.?.?","out country"))

##########################

# NOT INTERNATIONAL
No_prefix <- c("no.?","no know.?")
NoInt_TradeKeywords <- unique(c(paste(sapply(No_prefix,paste,Int_TradeKeywords)),
                                paste(sapply(Int_TradeKeywords,paste,No_prefix)),
                                "only local.?.?","domestic use only","domestic trad.?.?.? only","use.? local.?.?","consume.? local.?.?","subsistence only"))


################################################################


### DECISION TREE 3 & 4: IS THERE EVIDENCE THAT USE/TRADE (3) IS NOT or (4) IS A THREAT? ###

# IS A THREAT
Threat_prefix <- c("unsustainab.?.?.?.?.?","threat.?.?.?.?","exploit.?.?.?.?.?")

Threat_suffixterminal <- c("threat","problem","issue") #"concern"

Threat_suffix <- c(paste(sapply(c("considered","thought","likely","probabl.?"),paste,Threat_suffixterminal)))

# split out to allow for more strings to be searched
Threat_TradeKeywordsPref <- c(paste(sapply(Threat_prefix,paste,c(TradeKeywords,"captur.?.?.?","collect.?.?.?.?","over.?collect",
                                                                 "over.?exploit","over.?harvest","over.?utilis","over.?fish"))))
Threat_TradeKeywordsSuff <- c(paste(sapply(TradeKeywords,paste,Threat_suffix)))

##########################

# IS NOT A THREAT
NoThreat_prefix <- c(paste0(sapply(c("no known","no.?"),paste,c("threat.?.?.?.?"))))

NoThreat_TradeKeywords <- c(paste(sapply(NoThreat_prefix,paste,c(Trade_core,"collect.?.?.?.?"))),
                            paste(sapply(Trade_core,paste,NoTrade_Prefix)))

##########################

# MAY BE A THREAT
PotentialThreat_prefix <- c("potential.?.?","possibl.?","future","uncertain","may.?be")
PotentialThreat_prefix <- c(PotentialThreat_prefix,paste(sapply(c(PotentialThreat_prefix,"may"),paste,c("pose.?"))))

#PotentialThreat_prefix <- c(paste(sapply(PotentialThreat_prefix, paste, c("over.?harvest", "threat", "risk", "problem"))))

PotentialThreat_suffix <- c("potential","possible","future","uncertain")
PotentialThreat_suffix <- c(paste(sapply(c("may pose", "may.?be","may become"),paste,PotentialThreat_suffix)),
                            "may.?be","could become","may pose",PotentialThreat_suffix)

PotentialThreat_suffix <- c(paste(sapply(c(PotentialThreat_suffix,"may", "might", "might suffer"), paste, 
                                         c("risk","threat.?","additional threat.?","problem","impact.? population"))))


PotentialThreat_TradeKeywords <- c(paste(sapply(Trade_core,paste,PotentialThreat_suffix)),
                                   paste(sapply(c(PotentialThreat_suffix,PotentialThreat_prefix), paste, c(Trade_core,"over.?harvest", "over.?exploit"))))


#####################################################################################################################################################
#####################################################################################################################################################


## Score taxa based on the steps of the decision tree ##

#####################################################################################################################################################


# 1. Is there evidence that use and/or trade takes place? -----------
# (a) No/past/future/potential/possiblr
# (b) Yes/probable 
# if yes/probably keywords = 1(b) ELSE 1(a)

# evidence of use/trade
# i. subsistence, national or international > 0
# ii. keywords for use/trade in rationale, threat or use/trade
# iii. threat from BRU (excluding Past, unlikely to return)


trade$DT1 <-"1a.no" # default = no

trade$DT1<-apply(trade[,c(usetrade)],MARGIN=1, FUN=function(x) {ifelse(any(!is.na(x)),"1b.yes",trade$DT1)}) # (i)

trade$DT1[grepl(paste(c(TradeKeywords,UnknownTrade_Keywords),collapse="|"),trade$useTrade) & !grepl(paste(NoTrade_Keywords,collapse="|"),trade$useTrade)]<-"1b.yes" # (ii)
trade$DT1[grepl(paste(c(TradeKeywords,UnknownTrade_Keywords),collapse="|"),trade$rationale) & !grepl(paste(NoTrade_Keywords,collapse="|"),trade$rationale)]<-"1b.yes" # (ii)
trade$DT1[grepl(paste(c(TradeKeywords,UnknownTrade_Keywords),collapse="|"),trade$threats) & !grepl(paste(NoTrade_Keywords,collapse="|"),trade$threats)]<-"1b.yes" # (ii)

trade$DT1[(trade$kingdomName %in% c("ANIMALIA","CHROMISTA","FUNGI") & grepl(paste(threatcodes.animals,collapse="|"),trade$ThreatCodes)) | (trade$kingdomName=="PLANTAE" & grepl(paste(c(threatcodes.plants),collapse="|"),trade$ThreatCodes))] <- "1b.yes" # (iii)


#####################################################################################################################################################


# 2a. IF 1(a) THEN is there evidence that use and/or trade is a potential future threat? -----------
# (a) No == U
# (b) Yes == I

trade$DT2a[trade$DT1=="1a.no"]<-"2aa.no" # default = no

trade$DT2a[trade$DT1=="1a.no" & !!rowSums(sapply(trade[c(cn)], grepl, pattern = paste(FutureKeywords,collapse="|")))] <- "2ab.yes" #(ii)

trade$Score[which(trade$Score=="" & trade$DT1=="1a.no" & trade$DT2a=="2aa.no")]<-"U"
trade$Score[which(trade$Score=="" & trade$DT1=="1a.no" & trade$DT2a=="2ab.yes")]<-"I"

trade$BasisOfScore[which(trade$BasisOfScore=="" & trade$DT1=="1a.no" & trade$DT2a=="2aa.no")]<-"DecisionTree"
trade$BasisOfScore[which(trade$BasisOfScore=="" & trade$DT1=="1a.no" & trade$DT2a=="2ab.yes")]<-"DecisionTree"


#####################################################################################################################################################


# 2b. IF 1(b) THEN is there evidence that use and/or trade is NOT international -----------
# (a) No
# (b) Yes == U
# if 2b(b) ELSE 2b(a)
# 2b(b) = international is not ticked but national/subsistance is,  use/trade NOT international keywords

# evidence that use/trade is not international
# i. International_EndUses = NA AND subsistence|national > 0 
# ii. NoInt_TradeKeywords

trade$DT2b[trade$DT1=="1b.yes"]<-"2ba.no" # default = no

trade$DT2b[is.na(trade$International_EndUses) & (!is.na(trade$Subsistence)|!is.na(trade$National))]<-"2bb.yes" #(i)

trade$DT2b[!!rowSums(sapply(trade[c(cn)], grepl, pattern = paste(NoInt_TradeKeywords,collapse="|")))] <- "2bb.yes" #(ii)

trade$Score[which(trade$Score=="" & trade$DT2b=="2bb.yes")]<-"U"
trade$BasisOfScore[which(trade$BasisOfScore=="" & trade$DT2b=="2bb.yes")]<-"DecisionTree"


#####################################################################################################################################################


# 3. IF 2b(a) THEN is there evidence that use/trade is NOT a threat? -----------
# (a) No 
# (b) Yes == U
# if 3(b) ELSE 3(a)
# 3(b) = threat codes OTHER THAN use threat codes ticked (or only "Past, unlikely to return"), use/trade NO threat keywords

# evidence that use/trade is not a threat
# i. NoThreat_TradeKeywords

trade$DT3[trade$DT2b=="2ba.no"]<-"3a.no" # default = no

trade$DT3[!!rowSums(sapply(trade[c(cn)], grepl, pattern = paste(NoThreat_TradeKeywords,collapse="|")))] <- "3b.yes" #(i)


trade$Score[which(trade$Score=="" & trade$DT3=="3b.yes")]<-"U"
trade$BasisOfScore[which(trade$BasisOfScore=="" & trade$DT3=="3b.yes")]<-"DecisionTree"


#####################################################################################################################################################


# 4. IF 3(a) THEN is there evidence that use/trade IS a threat? -----------
# (a) No/past/future/potential/possible = I
# (b) Yes/probable
# if 4(b) ELSE 4(a)
# 4(b) = use-related threat code selected, use/trade as threat keywords

# evidence that use/trade is a threat
# i. use-related threat code selected (that is NOT coded up as "Future" or "Past unlikely to return")
# ii. Threat_TradeKeywords

trade$DT4[trade$DT3=="3a.no"]<-"4a.no" # default = no 

#trade$DT4[grepl(paste(threatcodes.universal,collapse="|"),trade$ThreatCodes) | (trade$Kingdom=="PLANTAE" & grepl(paste(c(threatcodes.plants),collapse="|"),trade$ThreatCodes))]<-"4b.yes"

# (b) yes/probable
trade$DT4[!!rowSums(sapply(trade[c(cn)], grepl, pattern = paste(c(Threat_TradeKeywordsPref,"over.?exploit","over.?harvest","over.?fish"),collapse="|")))] <- "4b.yes" #(ii)
trade$DT4[!!rowSums(sapply(trade[c(cn)], grepl, pattern = paste(Threat_TradeKeywordsSuff,collapse="|")))] <- "4b.yes" #(ii)

trade$DT4[grepl(paste(TradeKeywords,collapse="|"),trade$threats)] <- "4b.yes" #(ii) # test what happens when adding tradekeywords in threats without other threat text

# (a) potential/future threat
trade$DT4[!!rowSums(sapply(trade[c(cn)], grepl, pattern = paste(PotentialThreat_TradeKeywords,collapse="|")))] <- "4a.no" #(ii)

# (a) past threat
trade$DT4[is.na(trade$ThreatCodes) &
            ((trade$kingdomName %in% c("ANIMALIA","CHROMISTA","FUNGI") & grepl(paste(threatcodes.animals,collapse="|"),trade$ThreatCodes_PastUnlikely)) |
               (trade$kingdomName=="PLANTAE" & grepl(paste(c(threatcodes.plants),collapse="|"),trade$ThreatCodes_PastUnlikely)))] <- "4a.no" #(ii)


trade$Score[which(trade$Score=="" & trade$DT4=="4a.no")]<-"I"
trade$BasisOfScore[which(trade$BasisOfScore=="" & trade$DT4=="4a.no")]<-"DecisionTree"


#####################################################################################################################################################


# 5. IF 4(b) THEN is there evidence that use/trade is international? -----------
# (a) yes/probably = L
# (b) No/past/future/potential/possible = I

# evidence that the use/trade considered a threat is international
# i. If trade = threat in free text AND that kind of trade is in International_EndUses (e.g. "pet trade is a threat" and Pets/display animals is in International_EndUses)
# ii. Int_TradeKeywords

trade$DT5 <- NA
trade$DT5[trade$DT4=="4b.yes"]<-"5b.no" # default = no

# general international keywords
trade$DT5[trade$DT4=="4b.yes" & !!rowSums(sapply(trade[c(cn)], grepl, pattern = paste(c(Int_TradeKeywords),collapse="|")))] <- "5a.yes" #(ii) 

################################################################

# keywords for specific end uses #
# consider (a) end use AND (b) related text in use and trade, rationale, and threats. (trade keywords in threats, trade + threat keywords in other three)
# for use and trade, and rationale, also consider (c) threat code (text in threats is assumed to indicate threat)

# THREATS 

# pets #

# relevant threat trade keywords in any of the free text fields AND international end use
trade$DT5[!!rowSums(sapply(trade[c(cn)], function(x)
  grepl(paste(c(Threat_TradeKeywordsPref[grepl(paste(c("ornamental","pet","aquari","bird","falcon","fish"),collapse="|"),Threat_TradeKeywordsPref)],
                Threat_TradeKeywordsSuff[grepl(paste(c("ornamental","pet","aquari","bird","falcon","fish"),collapse="|"),Threat_TradeKeywordsSuff)], #"collect.?.?.? export","harvest trade","ornamental","aquari.?.?"),
                "over.?fish","over.?harvest"),collapse="|"),x) &
    grepl("Pets",trade$International_EndUses)))] <- "5a.yes"

# relevant trade keywords in threat AND international end use AND threatened by WT under DT4 (4b.yes)
trade$DT5[trade$DT4=="4b.yes" & c(grepl(paste(c(TradeKeywords[grepl(paste(c("ornamental","pet","aquari","bird","falcon","fish"),collapse="|"),TradeKeywords)],
                                                "over.?fish","over.?harvest","collect.?.?.? export","harvest trade","ornamental","aquari.?.?"),collapse="|"),trade$threats) & #
                                    grepl("Pets",trade$International_EndUses))] <- "5a.yes"


##############################

# food #

# relevant threat trade keywords in any of the free text fields AND international end use
trade$DT5[!!rowSums(sapply(trade[c(cn)], function(x){
  grepl(paste(c(Threat_TradeKeywordsPref[grepl(paste(c("fish","food"),collapse="|"),Threat_TradeKeywordsPref)],
                Threat_TradeKeywordsSuff[grepl(paste(c("fish","food"),collapse="|"),Threat_TradeKeywordsSuff)],
                "over.?fish","over.?harvest","exploit.?.?.?.?.? line.?fish"),collapse="|"),x) &
    grepl("Food",trade$International_EndUses)}))] <- "5a.yes"

# relevant trade keywords in threat AND international end use AND threatened by WT under DT4 (4b.yes)
trade$DT5[trade$DT4=="4b.yes" &  grepl(paste(c(TradeKeywords[grepl(paste(c("fish","food"),collapse="|"),TradeKeywords)],
                                               "over.?fish","over.?harvest","exploit.?.?.?.?.? line.?fish","trade food"),collapse="|"),trade$threats) &
            grepl("Food",trade$International_EndUses)] <- "5a.yes"


##############################

# Medicine #

# relevant threat trade keywords in any of the free text fields AND international end use
trade$DT5[!!rowSums(sapply(trade[c(cn)], function(x){
  grepl(paste(c(Threat_TradeKeywordsPref[grepl(paste(c("medicin.?.?"),collapse="|"),Threat_TradeKeywordsPref)],
                Threat_TradeKeywordsSuff[grepl(paste(c("medicin.?.?"),collapse="|"),Threat_TradeKeywordsSuff)]), #,"medicinal use"
              collapse="|"),x) & grepl("Medicine",trade$International_EndUses)}))] <- "5a.yes"

# relevant trade keywords in threat AND international end use AND threatened by WT under DT4 (4b.yes)
trade$DT5[trade$DT4=="4b.yes" &  grepl(paste(c(TradeKeywords[grepl(paste(c("medicin.?.?"),collapse="|"),TradeKeywords)]),collapse="|"),trade$threats) &
            grepl("Medicine",trade$International_EndUses)] <- "5a.yes"


##############################

# Handicrafts, jewellery, etc. #

# relevant threat trade keywords in any of the free text fields AND international end use
trade$DT5[!!rowSums(sapply(trade[c(cn)], function(x){
  grepl(paste(c(Threat_TradeKeywordsPref[grepl(paste(c("curio","cosmetic"),collapse="|"),Threat_TradeKeywordsPref)],
                Threat_TradeKeywordsSuff[grepl(paste(c("curio","cosmetic"),collapse="|"),Threat_TradeKeywordsSuff)]), #,"medicinal use"
              collapse="|"),x) & grepl("Handicraft",trade$International_EndUses)}))] <- "5a.yes"

# relevant trade keywords in threat AND international end use AND threatened by WT under DT4 (4b.yes)
trade$DT5[trade$DT4=="4b.yes" &  grepl(paste(c(TradeKeywords[grepl(paste(c("curio","cosmetic"),collapse="|"),TradeKeywords)]),collapse="|"),trade$threats) &
            grepl("Handicraft",trade$International_EndUses)] <- "5a.yes"


##############################

# Sport hunting/specimen collecting #

# relevant threat trade keywords in any of the free text fields AND international end use
trade$DT5[!!rowSums(sapply(trade[c(cn)], function(x){
  grepl(paste(c(Threat_TradeKeywordsPref[grepl(paste(c("sport"),collapse="|"),Threat_TradeKeywordsPref)],
                Threat_TradeKeywordsSuff[grepl(paste(c("sport"),collapse="|"),Threat_TradeKeywordsSuff)]), #,"medicinal use"
              collapse="|"),x) & grepl("hunting",trade$International_EndUses)}))] <- "5a.yes"

# relevant trade keywords in threat AND international end use AND threatened by WT under DT4 (4b.yes)
trade$DT5[trade$DT4=="4b.yes" &  grepl(paste(c(TradeKeywords[grepl(paste(c("sport"),collapse="|"),TradeKeywords)]),collapse="|"),trade$threats) &
            grepl("hunting",trade$International_EndUses)] <- "5a.yes"


################################################################

trade$Score[which(trade$DT5=="5a.yes" & !trade$BasisOfScore %in% c("Animal_PlantThreat", "TradeAndUse_keywords"))]<-"L"
trade$Score[which(trade$Score=="" & trade$DT5=="5b.no")]<-"I"

trade$BasisOfScore[which(trade$DT5=="5a.yes" & !trade$BasisOfScore %in% c("Animal_PlantThreat", "TradeAndUse_keywords"))]<-"DecisionTree"
trade$BasisOfScore[which(trade$BasisOfScore=="" & trade$DT5=="5b.no")]<-"DecisionTree"


#####################################################################################################################################################
#####################################################################################################################################################

# End