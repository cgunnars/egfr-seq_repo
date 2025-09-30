source('./src/utilities.R')

prepEnv()
library(readxl)
lux <- read_excel('./data/lux_data/fig4/clean-data.xlsx')

lux <- lux %>% group_by(condition, day, donor) %>% summarize(mean_mfi = mean(`norm mfi`))

ggplot(data = lux %>% filter(day == 1), aes(x=condition, y=mean_mfi, fill=condition)) + geom_boxplot(show.legend = FALSE) + geom_point(show.legend=FALSE) +
    scale_fill_manual(values=getDrugColors()) + theme_classic()

ggsave('./fig/growth/RNAseq_lux.pdf', width = 4, height=2, create.dir=T)
