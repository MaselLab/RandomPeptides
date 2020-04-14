SeqAny.filt <- SeqAny %>%
  mutate(FIT.neme = ifelse(is.na(Up.any), 0, 1)) %>% 
  select(PEP.ID = X, FIT.neme)
  
Pool.FIT <- inner_join(SeqAny.filt, FitTableMLE, by = "PEP.ID")

Pool.FIT <- Pool.FIT %>% 
  mutate(FIT.BIN = ifelse(FIT.EST > 1, 1, 0)) %>% 
  select(PEP.ID, FIT.neme, FIT.BIN, FIT.EST)

sum(Pool.FIT$FIT.BIN)
# This is the number of overlapping PEP.ID's that WE are convinced are beneficial
# [1] 280

Pool.UpFIT <- Pool.FIT %>% 
  filter(FIT.neme == 1)

sum(Pool.UpFIT$FIT.neme)
# This is the number of overlapping PEP.ID's that NEME is convinced are beneficial
# [1] 277

sum(Pool.UpFIT$FIT.neme == Pool.UpFIT$FIT.BIN)
# The number of PEP.ID fitnesses we agree on
# [1] 274

# there are 3 we therefore do not agree on

Excess.UpFIT <- Pool.UpFIT %>% 
  filter(FIT.BIN == 0)

#             PEP.ID FIT.neme FIT.BIN   FIT.EST
# 1 PEPNR00000000665        1       0 0.9817836
# 2 PEPNR00000000933        1       0 0.9675737
# 3 PEPNR00000001049        1       0 0.9935439

# Very close, barely fall out

Pool.DownFIT <- Pool.UpFIT <- Pool.FIT %>% 
  filter(FIT.neme == 0)

dim(Pool.DownFIT)
# The number of PEP.ID NEME says have a negative selection coefficient
# [1] 436   4

sum(Pool.DownFIT$FIT.neme == Pool.UpFIT$FIT.BIN)
# The number of PEP.ID's we agree have a negative selection coefficient
# [1] 430

# There are therefore 6 we do not agree on
Excess.DownFIT <- Pool.DownFIT %>% 
  filter(FIT.BIN == 1)

#            PEP.ID FIT.neme FIT.BIN  FIT.EST
# 1 PEPNR00000000051        0       1 1.004408
# 2 PEPNR00000000650        0       1 1.076314
# 3 PEPNR00000000668        0       1 1.002184
# 4 PEPNR00000000775        0       1 1.089601
# 5 PEPNR00000000807        0       1 1.036029
# 6 PEPNR00000000943        0       1 1.228099

# Again very close, except for the last one


# However the kicker is the number of PEP.ID's we find beneficial
sum(FitTableMLE$FIT.EST > 1)

# [1] 519

# Which identifies 242 possibly beneficial PEP.ID's that Neme didnt consider at all

FitTableMLE.UpFIT <- FitTableMLE %>% filter(FIT.EST > 1)

ggplot() +
  geom_histogram(mapping = aes(Pool.UpFIT$FIT.EST, fill = "red"), alpha = 0.6) + 
  geom_histogram(mapping = aes(FitTableMLE.UpFIT$FIT.EST, fill = "blue"), alpha = 0.6) +
  scale_fill_discrete(name = "Method", labels = c("MLE", "NEME")) +
  labs(x = "Fitness Estimate", y = "Count")
                 