.headers on

SELECT ga.geneId, ga.geneName, d.diseaseId, lower(d.diseaseName) FROM geneDiseaseNetwork gdn, geneAttributes ga, diseaseAttributes d WHERE ga.geneNID = gdn.geneNID AND d.diseaseNID = gdn.diseaseNID;



