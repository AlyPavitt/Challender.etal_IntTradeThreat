################################################################
################################################################

##### RL THREAT FROM INTERNATIONAL TRADE AUTOMATED SCORING ##### 

################################################################
################################################################

rm(list = ls())

################################################################

# Packages ----------

library("here")
library("reshape2")

################################################################

# read in data -----------
trade <- read.csv(paste0(here(),"RedListAssessment_data.csv")) # read in pre-prepared data extracted and combined from the Red List database

# read in "threats"
# downloaded from https://www.iucnredlist.org/
# note that if multiple downloads are taken from the IUCN Red List website (e.g. for different taxonomic subsets), you will need to rbind() them before running the script below

# threats (species threats)
threats <- read.csv(paste0(here(),"/Input/Threats.csv"))

################################################################
################################################################

####################### CLEAN FREE TEXT ########################

# key column names ----------
cn <- c("Threats","Rationale","Use.and.trade")
usetrade <- c("Subsistence","National","International")

# standardise case and remove unwanted characters -----------
trade[,c(cn)] <- sapply(trade[,c(cn)],function(x) tolower(x)) # to lower
trade[,c(cn)] <- sapply(trade[,c(cn)],function(x) gsub("<[^>]+>","",x)) # exclude everything between < and > 
trade[,c(cn)] <- sapply(trade[,c(cn)],function(x) gsub("[^a-zA-Z|^ ]", "", x)) # remove all characters except alpha and space
trade[,c(cn)] <- sapply(trade[,c(cn)],function(x) gsub("   |  ", " ", x)) # remove double and triple spaces

# standardise abbreviations -----------
trade[,c(cn)] <- sapply(trade[,c(cn)],function(x) gsub("cannot|can not", "cant", x)) 
trade[,c(cn)] <- sapply(trade[,c(cn)],function(x) gsub("is not", "isnt", x)) 
trade[,c(cn)] <- sapply(trade[,c(cn)],function(x) gsub("will not", "wont", x))
trade[,c(cn)] <- sapply(trade[,c(cn)],function(x) gsub("does not", "doesnt", x))
trade[,c(cn)] <- sapply(trade[,c(cn)],function(x) gsub("has not", "hasnt", x))
trade[,c(cn)] <- sapply(trade[,c(cn)],function(x) gsub("have not", "havent", x))
trade[,c(cn)] <- sapply(trade[,c(cn)],function(x) gsub("it is", "its", x))

# standardise other keywords ----------
trade[,c(cn)] <- sapply(trade[,c(cn)],function(x) gsub("utiliz", "utilis", x)) 

# remove stopwords ----------
stop <- c("a","also","an","and","any","are","as","be","been","being","by","for","from","has","have","if","in","is","it","of","or","still","that","the","there","this","to","were","where","which","with","would", # standard stops
          "species","wild") # specific stops
stop <- paste0(" ",c(stop)," ")
trade[,c(cn)] <- sapply(trade[,c(cn)],function(x) gsub(paste(stop,collapse="|"), " ", x))

head(trade[,c(cn)])

################################################################
################################################################

# Defining parameters of automated scoring ----------

# Four basis for automated scoring
# #1: "Use and Trade" field contains keywords indicating that use and trade are (a) unlikely to be a threat (score "U") or (b) there is insufficient information (score "I") 
# #2: "International trade is a significant threat" field is selected (score "L")
# #3: An animal was selected in the initial IUCN RL query on the basis of plant-related threat keywords or threat scores only (score "U")

################################################################
################################################################

################### #1 Use and Trade keywords ##################

################################################################

# Set keyword vectors #
# all lower case
# NB each ".?" indicates up to one character wildcard

# trade-related keywords -----------
TradeKeywords <- c("trad.?.?.?","harvest.?.?.?","export.?.?.?.?.?.?","exploit.?.?.?.?.?.?","utilis.?.?.?.?.?.?","trapp.?.?.?.?","collect.?.?.?.?")

################################################################

###### CODE U ######

# keywords indicating that there is no trade/ international trade or it is not a threat (U) -----------

NoTrade_Prefix <- c("no","not","not known","not found","no.? current.?.?","current.?.? no.?","no evidence","no evidence current","no.? direct.?.?","no existing",
                    "neither","appears not","no information suggests","not subject", "not thought use","not thought","no.? document.?.?","no.? apparent",
                    "no known specific","no.? known current.?.?","no.? know.?.?.?.?.?.?","no documentation",
                    "no.? report.?.?","no.? record.?.?","no record.?.?.?.?","no report.?.?.?.?","no.? re.?.?.?.?ed.?.?.?","no report.?.?.?.? know.?","never re.?.?.?.?ed",
                    "not believed","no.? thought","doesnt appear","doesnt appear use.?","unlikely","no.? significant.?.? use.?")

NoTrade_Keywords <- c(paste(sapply(NoTrade_Prefix,paste,c(TradeKeywords,"international trad","use"))),"not targeted specifically","none known",
                      "no known human use","no use.? know","no.? know.? use","no report.?.?.?.? use","no know.? record.?.?.?.? use","no know.? report.?.?.?.? use",
                      "no known removal trade","not found commercial pet trade","no way use.? trad...?","trade isnt apparent",
                      "trade unrecorded") 

###### CODE I ######

# keywords indicating that  there is insufficient information (I) to know whether trade is a threat -----------

UnknownTrade_Keywords <- c("no use.? information","no trad.?.?.? information","no use.?trade information","no.? known use","no.? known whether use",
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
                           "information regarding trade.?use isnt known","information.?.? missing","trade unrecord.?.?",
                           "no report.?.? information","no report.?.? use.? information","no record.?.? information use","no.? re.?.?.?.?.?.? information trade")

################################################################

# score based on keywords in use and trade that ---------- 
# U = indicate it is not in trade, or 
# I = indicate that it is unknown whether it is in trade due to insufficient information
trade$Score <- ""
trade$Score[which(grepl(paste(c(NoTrade_Keywords),collapse="|"),trade$Use.and.trade))] <- "U"
trade$Score[which(grepl(paste(c(UnknownTrade_Keywords),collapse="|"),trade$Use.and.trade))] <- "I"

# indicate basis of the score
trade$BasisOfScore <- ""
trade$BasisOfScore[which(grepl(paste(c(NoTrade_Keywords),collapse="|"),trade$Use.and.trade))] <- "TradeAndUse_keywords"
trade$BasisOfScore[which(grepl(paste(c(UnknownTrade_Keywords),collapse="|"),trade$Use.and.trade))] <- "TradeAndUse_keywords"

################################################################
################################################################

######## #2 International trade is a significant threat ########

################################################################

# #2: "International trade is a significant threat" field is selected (score "L") ----------
trade$Score[which(trade$Is.international.trade.a.significant.driver.of.threat...1.Yes..2.No..3.Unknown..blank.don.t.have.data. == "1")] <- "L"
trade$BasisOfScore[which(trade$Is.international.trade.a.significant.driver.of.threat...1.Yes..2.No..3.Unknown..blank.don.t.have.data. == "1")] <- "InternationalTradeisaThreat"

################################################################
################################################################

############ #3 Animal selected under plant threats ############

################################################################

# #4: An animal was selected in the initial IUCN RL query on the basis of plant-related threat keywords or threat scores only (score "U")

# a. initial keyword query ----------
#NB exclude "timber" from querykeywords2.pt1 and include it separately JUST for plants
querykeywords1 <- c("trade","export","collect","enthusiast","overharvest","overexploit","demand","market","cites","zoo","botanic","herbari")
querykeywords2.1 <- c("nternational","ransboundary","ntercontinental","border","pet","aquari","horticultur","commercial","world","ornamental","cage bird",
                      "cagebird","curio","medicin","fisher","regional","pharmaceutical","unsustainabl")
querykeywords2.2 <- c("use","utilis","utiliz")                   

# b. identify taxa meeting one of these keyword queries ----------
# search for timber separately to exclude animals that are indirectly threatened by timber harvest
trade$meet.initial.query_keywords1 <- !!rowSums(sapply(trade[c(cn)], grepl, pattern = paste(querykeywords1,collapse="|")))
trade$meet.initial.query_keywords2.pt1 <- !!rowSums(sapply(trade[c(cn)], grepl, pattern = paste(querykeywords2.1,collapse="|")))
trade$meet.initial.query_keywords2.pt2 <- !!rowSums(sapply(trade[c(cn)], grepl, pattern = paste(querykeywords2.2,collapse="|")))
trade$meet.initial.query_keywords2[trade$meet.initial.query_keywords2.pt1=="TRUE" & trade$meet.initial.query_keywords2.pt2=="TRUE"] <- "TRUE"

# c. identify animal taxa meeting one of the threat criteria (see under #3) ----------
trade$meet.initial.query_threats[trade$Species %in% querythreatspecies.animals]<-"TRUE"

# d. identify species that meet either the initial query keywords OR the initial query threats ----------
# with additional caveat of "timber" as a keyword and five threat codes relating to plants/timber only refering to plants
trade$Score[trade$Kingdom %in% c("ANIMALIA","CHROMISTA","FUNGI") & !trade$meet.initial.query_keywords1=="TRUE" &
              is.na(trade$meet.initial.query_keywords2) & is.na(trade$meet.initial.query_threats)] <- "U"

trade$BasisOfScore[trade$Kingdom %in% c("ANIMALIA","CHROMISTA","FUNGI") & !trade$meet.initial.query_keywords1=="TRUE" &
                     is.na(trade$meet.initial.query_keywords2) & is.na(trade$meet.initial.query_threats)] <- "Animal_PlantThreat"

################################################################

