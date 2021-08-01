
library(tidyverse)

odonates <- read_csv("datasets_derived/odonate_summaries_predicted.csv")

p <- ggplot(aes(x = abundance, y = ses_median - ses_5, colour = suborder),
            data = odonates) +
            geom_point()
ggsave("figures/abundance_v_specialization.svg", p)