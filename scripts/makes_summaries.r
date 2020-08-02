# Author: Jake Harvey - jakekharvey@gmail.com
library(xlsx)


main <- function(){
  
  checklist <- read_template()
  dataset <- read.xlsx("odonate2020raw.xlsx", sheetIndex = "All_data")
  
  for(sp in checklist$Species){
    
    sp  <- checklist[sp,]$Species
    gen <- checklist[sp,]$Genus
    checklist[sp,]$X2020_paras.     <- count_parasitized(gen, sp, dataset = dataset)
    checklist[sp,]$X2020_non-paras. <- count_non_parasitized(gen, sp, dataset = dataset)
  }
  
  
  
  
  
  
  
  
}

count_parasitized <- function(genus, species, dataset){
  
  y <- dataset[which(dataset$Genus == genus & dataset$Species == species & dataset$Mite_Prevelence == 1),]
  return(length(y$Species)) # just to select any column
}

count_non_parasitized <- function(genus, species, dataset){
  
  y <- dataset[which(dataset$Genus == genus & dataset$Species == species),]
  return(length(y$Species)) # just to select any column
}

read_template <- function(){
x<-"Suborder |Family            |Genus         |Species         |prevalence |total_paras |total_abundance |2020_paras. |2020_non-paras. |2019_paras. |2019_non-paras. |2015_paras. |2015 non-paras. |Distribution
Anisoptera   |Aeshnidae         |Aeshna        |canadensis      |           |            |                |            |                |0           |3               |2           |8               |BMF
Anisoptera   |Aeshnidae         |Aeshna        |clepsydra       |           |            |                |            |                |0           |0               |0           |0               |mf
Anisoptera   |Aeshnidae         |Aeshna        |constricta      |           |            |                |            |                |0           |0               |0           |0               |mF
Anisoptera   |Aeshnidae         |Aeshna        |eremita         |           |            |                |            |                |0           |0               |0           |0               |SBMF
Anisoptera   |Aeshnidae         |Aeshna        |interrupta      |           |            |                |            |                |0           |0               |0           |2               |SBMF
Anisoptera   |Aeshnidae         |Aeshna        |juncea          |           |            |                |            |                |0           |0               |0           |0               |SBm
Anisoptera   |Aeshnidae         |Aeshna        |septentrionalis |           |            |                |            |                |0           |0               |0           |0               |Sb
Anisoptera   |Aeshnidae         |Aeshna        |sitchensis      |           |            |                |            |                |0           |0               |0           |0               |SBMf
Anisoptera   |Aeshnidae         |Aeshna        |subarctica      |           |            |                |            |                |0           |0               |0           |0               |SBMf
Anisoptera   |Aeshnidae         |Aeshna        |tuberculifera   |           |            |                |            |                |0           |0               |0           |0               |SBmF
Anisoptera   |Aeshnidae         |Aeshna        |umbrosa         |           |            |                |            |                |0           |0               |0           |0               |SBMF
Anisoptera   |Aeshnidae         |Aeshna        |verticalis      |           |            |                |            |                |0           |0               |0           |0               |f
Anisoptera   |Aeshnidae         |Anax          |junius          |           |            |                |            |                |0           |0               |0           |0               |bmF
Anisoptera   |Aeshnidae         |Anax          |longipes        |           |            |                |            |                |0           |0               |0           |0               |f
Anisoptera   |Aeshnidae         |Basiaeschna   |janata          |           |            |                |            |                |0           |0               |0           |0               |bMF
Anisoptera   |Aeshnidae         |Boyeria       |grafiana        |           |            |                |            |                |0           |0               |0           |0               |BMF
Anisoptera   |Aeshnidae         |Boyeria       |vinosa          |           |            |                |            |                |0           |0               |0           |1               |MF
Anisoptera   |Aeshnidae         |Epiaeschna    |heros           |           |            |                |            |                |0           |0               |0           |0               |f
Anisoptera   |Aeshnidae         |Gomphaeschna  |furcillata      |           |            |                |            |                |0           |0               |0           |0               |f
Anisoptera   |Aeshnidae         |Nasiaeschna   |pentacantha     |           |            |                |            |                |0           |0               |0           |0               |f
Anisoptera   |Aeshnidae         |Rhionaeschna  |mutata          |           |            |                |            |                |0           |0               |0           |0               |f
Anisoptera   |Cordulegastridae  |Cordulegaster |diastatops      |           |            |                |            |                |0           |0               |0           |1               |BMF
Anisoptera   |Cordulegastridae  |Cordulegaster |maculata        |           |            |                |            |                |0           |0               |0           |1               |SBMF
Anisoptera   |Cordulegastridae  |Cordulegaster |obliqua         |           |            |                |            |                |0           |0               |0           |0               |F
Anisoptera   |Corduliidae       |Cordulia      |shurtleffii     |           |            |                |            |                |0           |3               |1           |15              |SBMF
Anisoptera   |Corduliidae       |Dorocordulia  |libera          |           |            |                |            |                |1           |15              |1           |7               |SBMF
Anisoptera   |Corduliidae       |Epitheca      |canis           |           |            |                |            |                |0           |0               |0           |4               |BMF
Anisoptera   |Corduliidae       |Epitheca      |cynosura        |           |            |                |            |                |0           |7               |0           |0               |MF
Anisoptera   |Corduliidae       |Epitheca      |princeps        |           |            |                |            |                |0           |0               |0           |0               |mF
Anisoptera   |Corduliidae       |Epitheca      |spinigera       |           |            |                |            |                |0           |16              |0           |0               |BMF
Anisoptera   |Corduliidae       |Helocordulia  |uhleri          |           |            |                |            |                |0           |0               |0           |2               |bMF
Anisoptera   |Corduliidae       |Neurocordulia |michaeli        |           |            |                |            |                |0           |0               |0           |0               |f
Anisoptera   |Corduliidae       |Neurocordulia |yamaskanensis   |           |            |                |            |                |0           |0               |0           |0               |mF
Anisoptera   |Corduliidae       |Somatochlora  |albicincta      |           |            |                |            |                |0           |0               |0           |17              |SBMf
Anisoptera   |Corduliidae       |Somatochlora  |brevicincta     |           |            |                |            |                |0           |0               |0           |0               |SB
Anisoptera   |Corduliidae       |Somatochlora  |cingulata       |           |            |                |            |                |0           |0               |0           |0               |SBMf
Anisoptera   |Corduliidae       |Somatochlora  |elongata        |           |            |                |            |                |0           |0               |0           |0               |BMf
Anisoptera   |Corduliidae       |Somatochlora  |filosa          |           |            |                |            |                |0           |0               |0           |0               |f
Anisoptera   |Corduliidae       |Somatochlora  |forcipata       |           |            |                |            |                |0           |0               |0           |1               |SBMF
Anisoptera   |Corduliidae       |Somatochlora  |franklini       |           |            |                |            |                |0           |0               |0           |0               |SBmf
Anisoptera   |Corduliidae       |Somatochlora  |incurvata       |           |            |                |            |                |0           |0               |0           |0               |mf
Anisoptera   |Corduliidae       |Somatochlora  |kennedyi        |           |            |                |            |                |0           |0               |0           |2               |SBMF
Anisoptera   |Corduliidae       |Somatochlora  |linearis        |           |            |                |            |                |0           |0               |0           |0               |f
Anisoptera   |Corduliidae       |Somatochlora  |minor           |           |            |                |            |                |0           |0               |0           |0               |SBMf
Anisoptera   |Corduliidae       |Somatochlora  |septentrionalis |           |            |                |            |                |0           |0               |0           |0               |SBm
Anisoptera   |Corduliidae       |Somatochlora  |tenebrosa       |           |            |                |            |                |0           |0               |0           |0               |mf
Anisoptera   |Corduliidae       |Somatochlora  |walshii         |           |            |                |            |                |0           |0               |0           |0               |SBMF
Anisoptera   |Corduliidae       |Somatochlora  |whitehousei     |           |            |                |            |                |0           |0               |0           |0               |SB
Anisoptera   |Corduliidae       |Somatochlora  |williamsoni     |           |            |                |            |                |0           |0               |0           |1               |BMF
Anisoptera   |Corduliidae       |Williamsonia  |fletcheri       |           |            |                |            |                |0           |0               |0           |0               |f
Anisoptera   |Gomphidae         |Arigomphus    |cornutus        |           |            |                |            |                |0           |0               |0           |0               |F
Anisoptera   |Gomphidae         |Arigomphus    |furcifer        |           |            |                |            |                |0           |0               |0           |0               |F
Anisoptera   |Gomphidae         |Dromogomphus  |spinosus        |           |            |                |            |                |0           |0               |0           |0               |MF
Anisoptera   |Gomphidae         |Gomphurus     |fraternus       |           |            |                |            |                |0           |0               |0           |0               |F
Anisoptera   |Gomphidae         |Gomphurus     |vastus          |           |            |                |            |                |0           |0               |0           |0               |F
Anisoptera   |Gomphidae         |Gomphurus     |ventricosus     |           |            |                |            |                |0           |0               |0           |0               |f
Anisoptera   |Gomphidae         |Hagenius      |brevistylus     |           |            |                |            |                |0           |0               |0           |1               |MF
Anisoptera   |Gomphidae         |Hylogomphus   |adelphus        |           |            |                |            |                |0           |0               |0           |0               |BMF
Anisoptera   |Gomphidae         |Lanthus       |parvulus        |           |            |                |            |                |0           |0               |0           |0               |bMf
Anisoptera   |Gomphidae         |Ophiogomphus  |anomalus        |           |            |                |            |                |0           |0               |0           |0               |mF
Anisoptera   |Gomphidae         |Ophiogomphus  |aspersus        |           |            |                |            |                |0           |0               |0           |0               |bMF
Anisoptera   |Gomphidae         |Ophiogomphus  |carolus         |           |            |                |            |                |0           |0               |0           |0               |bmf
Anisoptera   |Gomphidae         |Ophiogomphus  |colubrinus      |           |            |                |            |                |0           |0               |0           |0               |SBMF
Anisoptera   |Gomphidae         |Ophiogomphus  |mainensis       |           |            |                |            |                |0           |0               |0           |0               |Mf
Anisoptera   |Gomphidae         |Ophiogomphus  |rupinsulensis   |           |            |                |            |                |0           |0               |0           |0               |MF
Anisoptera   |Gomphidae         |Phanogomphus  |borealis        |           |            |                |            |                |0           |0               |0           |0               |MF
Anisoptera   |Gomphidae         |Phanogomphus  |descriptus      |           |            |                |            |                |0           |0               |0           |0               |bmF
Anisoptera   |Gomphidae         |Phanogomphus  |exilis          |           |            |                |            |                |0           |0               |0           |18              |BMF
Anisoptera   |Gomphidae         |Phanogomphus  |lividus         |           |            |                |            |                |0           |0               |0           |0               |bf
Anisoptera   |Gomphidae         |Phanogomphus  |spicatus        |           |            |                |            |                |0           |0               |5           |21              |BMF
Anisoptera   |Gomphidae         |Progomphus    |obscurus        |           |            |                |            |                |0           |0               |0           |0               |f
Anisoptera   |Gomphidae         |Stylogomphus  |albistylus      |           |            |                |            |                |0           |0               |0           |0               |bMF
Anisoptera   |Gomphidae         |Stylurus      |amnicola        |           |            |                |            |                |0           |0               |0           |0               |mf
Anisoptera   |Gomphidae         |Stylurus      |notatus         |           |            |                |            |                |0           |0               |0           |0               |F
Anisoptera   |Gomphidae         |Stylurus      |scudderi        |           |            |                |            |                |0           |0               |0           |0               |BMF
Anisoptera   |Gomphidae         |Stylurus      |spiniceps       |           |            |                |            |                |0           |0               |0           |0               |mF
Anisoptera   |Libellulidae      |Celithemis    |elisa           |           |            |                |            |                |0           |4               |1           |2               |F
Anisoptera   |Libellulidae      |Celithemis    |eponina         |           |            |                |            |                |0           |0               |0           |0               |f
Anisoptera   |Libellulidae      |Erythemis     |simplicicollis  |           |            |                |            |                |0           |0               |0           |0               |F
Anisoptera   |Libellulidae      |Erythrodiplax |berenice        |           |            |                |            |                |0           |0               |0           |0               |f
Anisoptera   |Libellulidae      |Ladona        |julia           |           |            |                |            |                |7           |63              |9           |11              |BMF
Anisoptera   |Libellulidae      |Leucorrhinia  |frigida         |           |            |                |            |                |26          |58              |12          |10              |MF
Anisoptera   |Libellulidae      |Leucorrhinia  |glacialis       |           |            |                |            |                |0           |1               |1           |4               |SBMF
Anisoptera   |Libellulidae      |Leucorrhinia  |hudsonica       |           |            |                |            |                |0           |1               |3           |15              |SBMF
Anisoptera   |Libellulidae      |Leucorrhinia  |intacta         |           |            |                |            |                |5           |1               |0           |5               |MF
Anisoptera   |Libellulidae      |Leucorrhinia  |patricia        |           |            |                |            |                |0           |0               |0           |7               |SBmf
Anisoptera   |Libellulidae      |Leucorrhinia  |proxima         |           |            |                |            |                |5           |1               |1           |4               |SBMF
Anisoptera   |Libellulidae      |Libellula     |incesta         |           |            |                |            |                |4           |26              |0           |2               |F
Anisoptera   |Libellulidae      |Libellula     |luctuosa        |           |            |                |            |                |0           |0               |0           |0               |F
Anisoptera   |Libellulidae      |Libellula     |pulchella       |           |            |                |            |                |0           |0               |0           |0               |bmF
Anisoptera   |Libellulidae      |Libellula     |quadrimaculata  |           |            |                |            |                |1           |32              |0           |8               |SBMF
Anisoptera   |Libellulidae      |Libellula     |semifasciata    |           |            |                |            |                |0           |0               |0           |0               |f
Anisoptera   |Libellulidae      |Nannothemis   |bella           |           |            |                |            |                |0           |0               |0           |0               |mf
Anisoptera   |Libellulidae      |Pachydiplax   |longipennis     |           |            |                |            |                |0           |0               |0           |0               |f
Anisoptera   |Libellulidae      |Pantala       |flavescens      |           |            |                |            |                |0           |0               |0           |0               |BMF
Anisoptera   |Libellulidae      |Pantala       |hymenaea        |           |            |                |            |                |0           |0               |0           |0               |bmf
Anisoptera   |Libellulidae      |Perithemis    |tenera          |           |            |                |            |                |0           |0               |0           |0               |f
Anisoptera   |Libellulidae      |Plathemis     |lydia           |           |            |                |            |                |0           |0               |0           |0               |MF
Anisoptera   |Libellulidae      |Sympetrum     |corruptum       |           |            |                |            |                |0           |0               |0           |0               |f
Anisoptera   |Libellulidae      |Sympetrum     |costiferum      |           |            |                |            |                |0           |0               |10          |3               |BMF
Anisoptera   |Libellulidae      |Sympetrum     |danae           |           |            |                |            |                |0           |0               |0           |1               |SBMF
Anisoptera   |Libellulidae      |Sympetrum     |internum        |           |            |                |            |                |0           |0               |0           |1               |SBMF
Anisoptera   |Libellulidae      |Sympetrum     |janeae          |           |            |                |            |                |0           |0               |0           |0               |F
Anisoptera   |Libellulidae      |Sympetrum     |obtrusum        |           |            |                |            |                |0           |0               |1           |53              |BMF
Anisoptera   |Libellulidae      |Sympetrum     |rubicundulum    |           |            |                |            |                |0           |0               |0           |0               |F
Anisoptera   |Libellulidae      |Sympetrum     |semicinctum     |           |            |                |            |                |0           |0               |0           |2               |MF
Anisoptera   |Libellulidae      |Sympetrum     |vicinum         |           |            |                |            |                |0           |0               |2           |17              |MF
Anisoptera   |Libellulidae      |Tramea        |lacerata        |           |            |                |            |                |0           |0               |0           |0               |f
Anisoptera   |Macromiidae       |Didymops      |transversa      |           |            |                |            |                |0           |0               |0           |5               |bMF
Anisoptera   |Macromiidae       |Macromia      |illinoiensis    |           |            |                |            |                |0           |0               |0           |0               |MF
Zygoptera    |Calopterygidae    |Calopteryx    |aequabilis      |           |            |                |            |                |0           |0               |0           |0               |BMF
Zygoptera    |Calopterygidae    |Calopteryx    |amata           |           |            |                |            |                |0           |0               |0           |0               |bMF
Zygoptera    |Calopterygidae    |Calopteryx    |maculata        |           |            |                |            |                |0           |7               |0           |1               |BMF
Zygoptera    |Calopterygidae    |Hetaerina     |americana       |           |            |                |            |                |0           |0               |0           |0               |f
Zygoptera    |Coenagrionidae    |Amphiagrion   |saucium         |           |            |                |            |                |0           |0               |0           |0               |MF
Zygoptera    |Coenagrionidae    |Argia         |apicalis        |           |            |                |            |                |0           |0               |0           |0               |f
Zygoptera    |Coenagrionidae    |Argia         |fumipennis      |           |            |                |            |                |0           |34              |25          |18              |bMF
Zygoptera    |Coenagrionidae    |Argia         |moesta          |           |            |                |            |                |0           |0               |0           |1               |MF
Zygoptera    |Coenagrionidae    |Chromagrion   |conditum        |           |            |                |            |                |2           |41              |8           |27              |BMF
Zygoptera    |Coenagrionidae    |Coenagrion    |resolutum       |           |            |                |            |                |0           |0               |0           |0               |SBMF
Zygoptera    |Coenagrionidae    |Coenagrion    |interrogatum    |           |            |                |            |                |0           |0               |1           |13              |SBMF
Zygoptera    |Coenagrionidae    |Enallagma     |anna            |           |            |                |            |                |0           |0               |0           |0               |f
Zygoptera    |Coenagrionidae    |Enallagma     |annexum         |           |            |                |            |                |0           |0               |1           |1               |SBMF
Zygoptera    |Coenagrionidae    |Enallagma     |antennatum      |           |            |                |            |                |0           |0               |0           |0               |F
Zygoptera    |Coenagrionidae    |Enallagma     |aspersum        |           |            |                |            |                |0           |0               |0           |0               |mF
Zygoptera    |Coenagrionidae    |Enallagma     |boreale         |           |            |                |            |                |1           |1               |16          |43              |SBMF
Zygoptera    |Coenagrionidae    |Enallagma     |carunculatum    |           |            |                |            |                |0           |2               |0           |0               |mF
Zygoptera    |Coenagrionidae    |Enallagma     |civile          |           |            |                |            |                |0           |0               |2           |4               |mf
Zygoptera    |Coenagrionidae    |Enallagma     |clausum         |           |            |                |            |                |0           |1               |0           |0               |mf
Zygoptera    |Coenagrionidae    |Enallagma     |ebrium          |           |            |                |            |                |64          |116             |34          |110             |BMF
Zygoptera    |Coenagrionidae    |Enallagma     |exsulans        |           |            |                |            |                |1           |0               |4           |6               |mF
Zygoptera    |Coenagrionidae    |Enallagma     |geminatum       |           |            |                |            |                |0           |6               |0           |1               |mF
Zygoptera    |Coenagrionidae    |Enallagma     |hageni          |           |            |                |            |                |32          |101             |39          |77              |SBMF
Zygoptera    |Coenagrionidae    |Enallagma     |signatum        |           |            |                |            |                |0           |0               |0           |0               |F
Zygoptera    |Coenagrionidae    |Enallagma     |traviatum       |           |            |                |            |                |0           |0               |1           |2               |f
Zygoptera    |Coenagrionidae    |Enallagma     |vernale         |           |            |                |            |                |1           |0               |8           |35              |bMF
Zygoptera    |Coenagrionidae    |Enallagma     |vesperum        |           |            |                |            |                |0           |0               |0           |0               |mF
Zygoptera    |Coenagrionidae    |Ischnura      |hastata         |           |            |                |            |                |0           |0               |0           |0               |f
Zygoptera    |Coenagrionidae    |Ischnura      |posita          |           |            |                |            |                |7           |5               |0           |0               |F
Zygoptera    |Coenagrionidae    |Ischnura      |verticalis      |           |            |                |            |                |35          |149             |63          |92              |BMF
Zygoptera    |Coenagrionidae    |Nehalennia    |gracilis        |           |            |                |            |                |0           |0               |0           |0               |MF
Zygoptera    |Coenagrionidae    |Nehalennia    |irene           |           |            |                |            |                |0           |0               |6           |44              |BMF
Zygoptera    |Lestidae          |Lestes        |congener        |           |            |                |            |                |0           |5               |1           |6               |BMF
Zygoptera    |Lestidae          |Lestes        |disjunctus      |           |            |                |            |                |1           |42              |5           |39              |SBMF
Zygoptera    |Lestidae          |Lestes        |dryas           |           |            |                |            |                |0           |0               |0           |0               |SBMF
Zygoptera    |Lestidae          |Lestes        |eurinus         |           |            |                |            |                |0           |0               |0           |0               |MF
Zygoptera    |Lestidae          |Lestes        |forcipatus      |           |            |                |            |                |0           |0               |0           |0               |BMF
Zygoptera    |Lestidae          |Lestes        |inaequalis      |           |            |                |            |                |3           |3               |0           |0               |mF
Zygoptera    |Lestidae          |Lestes        |rectangularis   |           |            |                |            |                |0           |0               |0           |0               |MF
Zygoptera    |Lestidae          |Lestes        |unguiculatus    |           |            |                |            |                |0           |0               |0           |0               |bMF
Zygoptera    |Lestidae          |Lestes        |vigilax         |           |            |                |            |                |0           |0               |6           |8               |f" 

x <- read.table(text = x, header = TRUE, sep = "|", stringsAsFactors = FALSE) 
return(x)
}

main()
