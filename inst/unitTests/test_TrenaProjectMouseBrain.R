library(TrenaProjectMouseBrain)
library(RUnit)
library(trenaSGM)
library(org.Hs.eg.db)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tp")) {
   message(sprintf("--- creating instance of TrenaProjectMouseBrain"))
   tp <- TrenaProjectMouseBrain();
}
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   # test_supportedGenes()
   # test_variants()
   # test_footprintDatabases()
   # test_expressionMatrices()
   # test_setTargetGene()
   # test_buildSingleGeneModel()
   # test_buildSingleGeneModel_slowGenes()
   # test_buildSingleGeneModel_footprintsAndWithout_MEF2C()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   message(sprintf("--- test_constructor"))

   checkTrue(all(c("TrenaProjectMouseBrain", "TrenaProjectMM10") %in% is(tp)))
   checkEquals(getFootprintDatabasePort(tp), 5432)

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
