# Challender.etal_IntTradeThreat
The protection of species threatened by international trade


# Trade-threatened species SI 2..4_automated scoring---------------------------
Script used as one step in the main data analysis process - alongside Red List database query and manual coding. 

Using data extracted and formated by the initial query of the Red List database, this script scored species as "L" (likely threatened by international trade), "I" (insufficient information) or "U" (unlikely to be threatened by international trade) based on three rules summarised below (full detail provided in 2.4 of the supplementary information).
Species that could not be automatically scored based on these parameters were manually coded. 

Automated coding rules
 #1: "Use and Trade" field contains keywords indicating that use and trade are (a) unlikely to be a threat (score "U") or (b) there is insufficient information (score "I") 
 #2: "International trade is a significant threat" field is selected (score "L")
 #3: An animal was selected in the initial IUCN RL query on the basis of plant-related threat keywords or threat scores only (score "U")


# Full_AutoSelection---------------------------
Script developed to automate the entire process presented in the main body of the paper.

For this method to be repeatable into the future as new/updated Red List assessments are published, the entire approach was automated to score taxa as "L" (likely threatened by international trade), "I" (insufficient information) or "U" (unlikely to be threatened by international trade) based directly on publicly available raw downloads from the Red List website (code works on downloads of the following files from https://www.iucnredlist.org/ : "assessment", "taxonomy", "threats" and "usetrade").

This combines the coding rules detailed in Trade-threatened species SI 2..4_automated scoring and used in the main paper analysis, with an additional set of rules to recreate the manual assessment approach carried out by authors (detailed in 2.7 of the supplementary information).

Currently this has been fully developed for two classes (Actinopterygii and Amphibia), but will work for all Animalia (noting that there are not taxon-specific commands or keywords for any other classes, and so scoring for classes other than Actinopterygii and Amphibia will be less accurate - see 4.2 of supplementary information for more details). Work is ongoing (beyond the scope of this paper) to finalise this process for all taxa. 
