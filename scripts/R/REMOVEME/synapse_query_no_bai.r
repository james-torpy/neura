library(rjson)
library(RCurl)
library(digest)
library(RUnit)
library(synapseClient)

projectId = "syn4645130"
cohort = "MSSM-Penn-Pitt"
queryString = sprintf("select id, name, cohort, dataType from file where projectId=='%s' and cohort=='%s' and fileType != '.bai'", projectId, cohort)
queryResults = synQuery(queryString, blockSize=200)
queryResults = queryResults$collectAll()