# library
library(tidyverse)
library(wesanderson)
library(ggsci)
library(scales)

# set theme
theme_set(theme_bw()+ 
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                strip.background = element_blank(),
                panel.border = element_rect(size=0.5,colour="black"),
                axis.text = element_text(color="black")
                )
          )
pal <- wes_palette("FantasticFox1")
heat_pal <- c(pal[3],pal[2],pal[4],pal[5])
pal <- wes_palette("Rushmore1")

pal_cont_polar <- wes_palette("Zissou1", 21, type = "continuous")

# functions
plotter <- function(print_obj,pre = "temp",out_dir = "~/Dropbox/Data/tmp",h = 4, w = 6){
  # if not interactive, just print it
  if( ! interactive()) {
    print(print_obj)
  } else {
    # if interactive, I want to export this to a pdf file
    outpath <- paste0(out_dir,"/",pre,".pdf")
    message(paste("outputting plot to", outpath))
    pdf(outpath,height = h, width = w, useDingbats = F)
    print(print_obj)
    dev.off()
  }
}


repeatColors =c("DNA"="#C19936",
                "DNA?"="#C19935",
                "LINE"="#FFA500",
                "Low_complexity"="#75B043",
                "LTR"="#51B756",
                "RC"="#53BB74",
                "Retroposon"="#55BE9D",
                "RNA"="#ff4000",
                "rRNA"="#52BEBB",
                "scRNA"="#6B9AD3",
                "Simple_repeat"="#8992C9",
                "SINE"="#9A8AC1",
                "snRNA"="#A886BC",
                "srpRNA"="#B67EB6",
                "Unspecified"="#C378B2",
                "tRNA"="#006400",
                "Unknown"="#BA55D3",
                "Satellite"="#53B0E4")

defaultColor = "#000080"

censatColors =c("TE" = "#E87C71",
                "MER"="#E28455",
                "HOR"="#D78C32",
                "BSAT"="#E370AB",
                "CER" = "#CE9334",
                "HSAT2"="#C19935",
                "HSAT1"="#A2A638",
                "HSAT3"="#8CAC3E",
                "L1"="#75B042",
                "LSAU"="#54B346",
                "LTR"="#51B756",
                "MST"="#53BB73",
                "GSAT"="#55BE8D",
                "GSATII"="#54C0A5",
                "rRNA"="#52BEBB",
                "SAR"="#51BDCE",
                "novel"="#9400D3",
                "HSAT4"="#53B0E4",
                "SATR"="#5AA5DA",
                "CT"="#6B9AD2",
                "HERV"="#8992C8",
                "MSAT"="#9A8AC2",
                "MON"="#A885BC",
                "SST"="#C378B2",
                "HSAT5"="#ED72A5",
                "HSAT6"="#EF768C", 
                "gap-rDNA"="#ff4000",
                "L1" = "#ffbf00", 
                "TAR"= "#0080ff",
                "ACRO"="#9400D4",
                "Alu"="#9A8AC3")
